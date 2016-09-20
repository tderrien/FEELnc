#!/bin/bash

### Script to run the FEELnc pipeline in one command line using shuffle mRNA
### Usage: FEELnc_pipeline.sh --candidate=<TRANSCRIPT_MODEL_GTF> --reference=<REFERENCE_GTF> --refSequence=<REFERENCE_SEQUENCE_FASTA> --name=<ANALYSIS_NAME>
###        These options are mandatories


### Chack if the three FEELnc modules are in the PATH
command -v FEELnc_filter.pl     >/dev/null 2>&1 || { echo -e >&2 "FEELnc_filter.pl is not in accessible, please modify your \$PATH variable accordingly\nExit";     exit 1;}
command -v FEELnc_codpot.pl     >/dev/null 2>&1 || { echo -e >&2 "FEELnc_codpot.pl is not in accessible, please modify your \$PATH variable accordingly\nExit";     exit 1;}
command -v FEELnc_classifier.pl >/dev/null 2>&1 || { echo -e >&2 "FEELnc_classifier.pl is not in accessible, please modify your \$PATH variable accordingly\nExit"; exit 1;}


usage="Usage: FEELnc_pipeline.sh --candidate=<TRANSCRIPT_MODEL_GTF> --reference=<REFERENCE_GTF> --refSequence=<REFERENCE_SEQUENCE_FASTA> --outname=<OUTPUT_NAME> --outdir=<OUTPUT_DIRECTORY>\nThese options are mandatories"

for opt in "$@"
do
    case $opt in
	--candidate=*)
	    candidate="${opt#*=}"
	    shift # past argument=value
	    ;;
	--reference=*)
	    reference="${opt#*=}"
	    shift # past argument=value
	    ;;
	--refSequence=*)
	    refSequence="${opt#*=}"
	    shift # past argument=value
	    ;;
	--outname=*)
	    outname="${opt#*=}"
	    shift # past argument=value
	    ;;
	--outdir=*)
	    outdir="${opt#*=}"
	    shift # past argument=value
	    ;;
	*)
            echo -e "Unknown option $opt\nExit\n$usage" # unknown option
	    exit 1
	    ;;
    esac
done


### Test if candidate, reference and refSequence options are provided
if [[ $candidate == "" ]]
then
    echo -e "Option --candidate is empty, it is mandatory\nExit\n$usage"
    exit 1
fi
if [[ $reference == "" ]]
then
    echo -e "Option --reference is empty, it is mandatory\nExit\n$usage"
    exit 1
fi
if [[ $refSequence == "" ]]
then
    echo -e "Option --refSequence is empty, it is mandatory\nExit\n$usage"
    exit 1
fi
if [[ $outname == "" ]]
then
    echo -e "Option --outname is empty, it is mandatory\nExit\n$usage"
    exit 1
fi
if [[ $outdir == "" ]]
then
    echo -e "Option --outdir is empty, it is mandatory\nExit\n$usage"
    exit 1
fi


### Test if candidate, reference and refSequence files are exist/are provided
if [[ ! (-s $candidate) ]]
then
    echo -e "The candidate file: '$candidate' didn't exists or is empty\nExit\n$usage"
    exit 1
fi
if [[ ! (-s $reference) ]]
then
    echo -e "The reference file: '$reference' didn't exists or is empty\nExit\n$usage"
    exit 1
fi
if [[ ! (-s $refSequence) ]]
then
    echo -e "The refSequence file: '$refSequence' didn't exists or is empty\nExit\n$usage"
    exit 1
fi
   

### Create output directory
echo -e "####################\nCreate output directory $outdir"
mkdir -p $outdir
echo



### Run the FILTER module
echo -e "####################\nStep 1/3: Running the FEELnc filter module:"

mkdir -p $outdir"/filter/"

afterFilter="$outdir/filter/$outname.filter.gtf"
afterFilterLog="$outdir/filter/$outname.filter.log"

cmd="FEELnc_filter.pl --infile=$candidate --mRNAfile=$reference --biotype=transcript_biotype=protein_coding --monoex=-1 --outlog=$afterFilterLog > $afterFilter"
echo -e "\ncommande line:\n\t$cmd\n"

FEELnc_filter.pl --infile=$candidate --mRNAfile=$reference --biotype=transcript_biotype=protein_coding --monoex=-1 --outlog=$afterFilterLog > $afterFilter
echo



### Run the CODPOT module
echo -e "####################\nStep 2/3: Running the FEELnc codpot module:"

codpotdir="$outdir/codpot/"
mkdir -p $codpotdir

afterCodpot="$outname.codpot"

cmd="FEELnc_codpot.pl --infile=$afterFilter --mRNAfile=$reference --genome=$refSequence -b transcript_biotype=protein_coding --mode=shuffle --outname=$afterCodpot --outdir=$codpotdir"
echo -e "\ncommande line:\n\t$cmd\n"

FEELnc_codpot.pl --infile=$afterFilter --mRNAfile=$reference --genome=$refSequence -b transcript_biotype=protein_coding --mode=shuffle --outname=$afterCodpot --outdir=$codpotdir
echo



### Run the CLASSIFIER module
echo -e "####################\nStep 3/3: Running the FEELnc classifier module:"

mkdir -p $outdir"/classifier/"

afterClass="$outdir/classifier/$outname.classifier.txt"
afterClassLog="$outdir/classifier/$outname.classifier.log"

cmd="FEELnc_classifier.pl --lncrna=$codpotdir/$afterCodpot.lncRNA.gtf --mrna=$reference --log=$afterClassLog > $afterClass"
echo -e "\ncommande line:\n\t$cmd\n"

FEELnc_classifier.pl --lncrna=$codpotdir"/"$afterCodpot".lncRNA.gtf" --mrna=$reference --log=$afterClassLog > $afterClass
echo


### Done
echo -e "\n####################\nDone: results in $outdir/{filter;codpot;classifier}"
echo
