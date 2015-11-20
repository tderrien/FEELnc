# FEELnc
 FlExible Extraction of Long non-coding RNAs

*Version (10/11/2015): Major Update. Correction of bug  during ORF sequences extraction from a GTF containing CDS information (the sequence was containing the RNA sequence in addition to the ORF). If you launched FEELnc analysis with a GTF including CDS information, please re-run your analysis with this new version. Sorry for the inconvenience.*

*Version (05/11/2015): Major update, replace ORF length feature by ORF coverage feature (ORF length divided by RNA length)*

*Version (26/05/2015): Major update includes the use of RandomForest for FEELnc_codpot*


Please note that FEELnc project is still in working progress. But feel free to use it and send any comments/bug/suggestions. Thanks!


## Introduction

This document is intended to give a (minimal) description of the FEELnc pipeline in order to annotate long non-coding RNAs (lncRNAs).

Currently, FEELnc is composed of 3 modules (See *Launch FEELnc 3-step pipeline* for more details):

	* FEELnc_filter.pl	: Extract, filter candidate transcripts
	* FEELnc_codpot.pl	: Compute the coding potential of candidate transcripts
	* FEELnc_classifier.pl: Classify lncRNAs based on their genomic localization wrt others transcripts


To get help on each module, you can type :

	FEELnc_filter.pl --help
	# Or
    FEELnc_filter.pl --man


## Input files

The formats used to describe genes, transcripts, exon is **.GTF** and **.FASTA** for genome file.

Basically, FEELnc users should have the following minimal input files:

	- Infile.GTF          (-i,--infile)   : input GTF file (e.g cufflinks transcripts.GTF)
	- ref_annotation.GTF  (-a,--mRNAfile) : GTF annotation file*
	- ref_genome.FASTA    (-g,--genome)   : genome FASTA file or directory with individual chrom FASTA files


\* *Note: It is recommended to only extract protein_coding transcripts (mRNAs) from the reference annotation file (ref_annotation.GTF) when this information is available, either manually or better by using the option :*
**--biotype transcript_biotype=protein_coding**.
In doing so, you will not remove lncRNAs overlapping other non-coding RNAs or pseudogenes....


-------------------------
## Installation and requirements


### Requirements

The following software and libraries must be installed on your machine:

- [Perl5+](https://www.perl.org/) : tested with version 5.18.2
 * [Bioperl](http://www.bioperl.org/wiki/Main_Page)  : tested with version BioPerl-1.6.924
 * [Parralell::ForkManager](http://search.cpan.org/~szabgab/Parallel-ForkManager-1.07/lib/Parallel/ForkManager.pm) : tested with version 1.07
- R [Rscript](http://cran.r-project.org): tested with version 3.1.0.
 * [ROCR](https://rocr.bioinf.mpi-sb.mpg.de/) test with version 1.0-5
 * [randomForest](http://cran.r-project.org/web/packages/randomForest/index.html) tested with version 4.6-10

* Note: R librairies should be installed automatically when running FEELnc. In case it does not work, please type in a R session:
	install.packages('ROCR')
	install.packages('randomForest')

### Installation

Clone the FEELnc git:

	git clone https://github.com/tderrien/FEELnc.git

Go to FEELnc directory

	cd FEELnc

Export PERL5LIB and FEELNCPATH variables

	export FEELNCPATH=${PWD}
	export PERL5LIB=$PERL5LIB:${FEELNCPATH}/lib/

Add FEELnc scripts to your PATH and add the distribution-specific binary of KmerInShort to your PATH or copy it to your bin directory

	export PATH=$PATH:${FEELNCPATH}/scripts/

	# for MAC
	export PATH=$PATH:${FEELNCPATH}/bin/MAC/
	# or
	cp ${FEELNCPATH}/bin/MAC/ ~/bin/

	# for LINUX
	export PATH=$PATH:${FEELNCPATH}/bin/LINUX/
	# or
	cp ${FEELNCPATH}/bin/LINUX/ ~/bin/

### Test with toy example:

	cd test/

	# Filter
	FEELnc_filter.pl -i transcript_chr38.gtf -a annotation_chr38.gtf \
    -b transcript_biotype=protein_coding > candidate_lncRNA.gtf

	# Coding_Potential
	# Note1 :  as a test, the training is only done on  1000 tx (-n 1000 option)
	FEELnc_codpot.pl -i candidate_lncRNA.gtf -a annotation_chr38.gtf -g genome_chr38.fa -n 1000

	# Classifier
	FEELnc_classifier.pl -i feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf -a annotation_chr38.gtf > candidate_lncRNA_classes.txt


### Head of the principal output files on the toy exemple:
**Filter**

*candidate_lncRNA.gtf*: a GTF file containing each exon of all transcripts that passed the filter.

	head -n2 candidate_lncRNA.gtf
	 38	Cufflinks	exon	516103	518188	.	+	.	gene_id "XLOC_090599"; transcript_id "TCONS_00231416"; class_code "u"; exon_number "1"; oId "CUFF.138034.16"; tss_id "TSS139476";
	 38	Cufflinks	exon	519376	519454	.	+	.	gene_id "XLOC_090599"; transcript_id "TCONS_00231416"; class_code "u"; exon_number "2"; oId "CUFF.138034.16"; tss_id "TSS139476";

*transcript_chr38.feelncfilter.log*: the log of the filter, write each transcripts that have been removed and the reason.

	head -n4 transcript_chr38.feelncfilter.log
	 COMMAND
	 /usr/bin/perl -w ${FEELNCPATH}/FEELnc_filter.pl -i transcript_chr38.gtf -a annotation_chr38.gtf -b transcript_biotype=protein_coding
	 Filter monoexonic (option 0): TCONS_00234417 =  1 exon (with strand .)...
	 Filter monoexonic (option 0): TCONS_00234233 =  1 exon (with strand .)...


**Coding Potential**

*candidate_lncRNA.gtf.lncRNA.gtf*: extraction of the original GTF file containing the transcripts predicted as lncRNAs.

	head -n2 candidate_lncRNA.gtf.lncRNA.gtf
	 38	Cufflinks	exon	516103	518188	.	+	.	gene_id "XLOC_090599"; transcript_id "TCONS_00231416"; class_code "u"; exon_number "1"; oId "CUFF.138034.16"; tss_id "TSS139476";
	 38	Cufflinks	exon	519376	519454	.	+	.	gene_id "XLOC_090599"; transcript_id "TCONS_00231416"; class_code "u"; exon_number "2"; oId "CUFF.138034.16"; tss_id "TSS139476";

*candidate_lncRNA.gtf.mRNA.gtf*: same as the previous file but for predicted mRNA transcripts.

	head -n2 candidate_lncRNA.gtf.mRNA.gtf
	 38	Cufflinks	exon	23052109	23052291	.	+	.	gene_id "XLOC_090839"; transcript_id "TCONS_00232417"; class_code "u"; exon_number "1"; oId "CUFF.139754.1"; tss_id "TSS139994";
	 38	Cufflinks	exon	23053697	23053936	.	+	.	gene_id "XLOC_090839"; transcript_id "TCONS_00232417"; class_code "u"; exon_number "2"; oId "CUFF.139754.1"; tss_id "TSS139994";

*candidate_lncRNA.gtf_RF_summary.txt*: file containing the threshold use for the prediction and the number of predicted lncRNAs and mRNAs.

	cat candidate_lncRNA.gtf_RF_summary.txt
	 # Summary file:
	 -With_cutoff:	0.526
	 -Nb_lncRNAs:	282
	 -Nb_mRNAs:	59

*candidate_lncRNA.gtf_RF.txt*: the result file of the random forest prediction. The kmerScore_* columns are the kmer score for the transcript for each kmer size, higer values representing a high occurence of common mRNA kmer. The ORF_cover is the ratio between the ORF size and the RNA size. The coding_potential column is the ratio between the number of trees who have vote for the transcript to be a mRNAand the total number of trees. The label column is the biotype prediction of the transcript (0: non coding, 1: coding) depending on is coding_potential and the threshold (see *candidate_lncRNA.gtf_RF_summary.txt*).

	head -n3 candidate_lncRNA.gtf_RF.txt | column -t
	 name            kmerScore_1mer  kmerScore_2mer  kmerScore_3mer  kmerScore_6mer  kmerScore_9mer  kmerScore_12mer  ORF_cover           RNA_size  coding_potential  label
	 TCONS_00231401  0.489996        0.478795        0.447059        0.353341        0.355575        0.490291         0.0234859675036928  13540     0.136             0
	 TCONS_00231402  0.485523        0.467847        0.435232        0.333433        0.203866        0.360465         0.0155639755173419  17155     0.114             0


**Classifier**

*candidate_lncRNA_classes.txt*: the file report statistics on the interaction and every interactions between lncRNA and RNA contained in the input annotation. Lines begining with a '#' are comments and lines begining with a '*' are the best interactions for a specific lncRNA.

	head -n13 candidate_lncRNA_classes.txt
	 #FEELnc Classification
	 #lncRNA file :  lncrna : feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf
	 #mRNA file : annotation_chr38.gtf
	 #Minimal window size : 10000
	 #Maximal window size : 100000
	 #Number of lncRNA : 282
	 #Number of mRNA : 293
	 #Number of interaction : 345
	 #Number of lncRNA without interaction : 30
	 #List of lncRNA without interaction : TCONS_00231718 TCONS_00232885 TCONS_00231717 TCONS_00231711 TCONS_00232954 TCONS_00231758 TCONS_00232925 TCONS_00231770 TCONS_00231712 TCONS_00232969 TCONS_00232893 TCONS_00231948 TCONS_00231773 TCONS_00232894 TCONS_00232805 TCONS_00231757 TCONS_00232958 TCONS_00232876 TCONS_00232871 TCONS_00231772 TCONS_00232968 TCONS_00232872 TCONS_00232806 TCONS_00231962 TCONS_00232926 TCONS_00231942 TCONS_00232804 TCONS_00232808 TCONS_00231769 TCONS_00231774
	 #INTERACTIONS
	 *lncRNA	XLOC_090640	TCONS_00231611	RNA_partner	ENSCAFG00000010206	ENSCAFT00000016203	antisense	intergenic	1555	 Status=convergent	 Subtype=downstream
	 *lncRNA	XLOC_090615	TCONS_00231491	RNA_partner	ENSCAFG00000009901	ENSCAFT00000046245	sense	intergenic	44590	 Status=same_strand	 Subtype=upstream


For more details, see the specific parts on each modules


-------------------------
## Launch the 3-step pipeline

### 1- FEELnc_filter.pl

The first step of the pipeline (FEELnc_filter) consists in filtering out unwanted/spurious transcripts and/or transcripts overlapping (in sense) exons of the reference annotation
and especially protein_coding exons as they more probably correspond to new mRNA isoforms (see -b,--biotype option).

	# Usage:
    FEELnc_filter.pl -i infile.gtf -a annotation_mRNA.gtf > candidate_lncRNA.gtf


If your reference annotation ("*ref_annotation.GTF*") contains transcript_biotype information (e.g protein_coding, pseudogene, miRNA...), you can subselect a specific transcript biotype to make the overlap with.

    FEELnc_filter.pl -i infile.gtf \
	-a ref_annotation.GTF \
	-b transcript_biotype=protein_coding \
	> candidate_lncRNA.gtf

This option is highly recommended if you don't want to remove transcripts
overlapping with other transcripts than mRNAs (e.g lincRNA, miRNA, pseudogene...).
For stranded RNASeq protocol, it is also possible  to include monoexonic lncRNAs that are antisense to mRNAs e.g

	FEELnc_filter.pl -i infile.gtf \
	-a ref_annotation.GTF \
	-b transcript_biotype=protein_coding \
	--monoex=-1
	> candidate_lncRNA.gtf


**- FULL OPTIONS (FEELnc_filter.pl --help) :**
```
  * General:
      --help                Print this help
      --man                 Open man page
      --verbosity           Level of verbosity

  * Mandatory arguments:
      -i,--infile=file.gtf          Specify the GTF file to be filtered (such as a cufflinks transcripts/merged .GTF file)
      -a,--mRNAfile=file.gtf        Specify the annotation GTF file to be filtered on based on sense exon overlap (file of protein coding annotation or whole reference annotation 'ref_annotation.GTF')

  * Filtering arguments:
      -s,--size=200                 Keep transcript with a minimal size (default 200)
      -b,--biotype                  Only consider transcript(s) from the reference annotation having this(these) biotype(s) (e.g : -b transcript_biotype=protein_coding,pseudogene) [default undef i.e all transcripts]
      -l,--linconly                 Keep only long intergenic/interveaning ncRNAs [default FALSE].
      --monoex=-1|0|1               Keep monoexonic transcript(s): mode to be selected from : -1 keep monoexonic antisense (for RNASeq stranded protocol), 1 keep all monoexonic, 0 remove all monoexonic   [default 0]
      --biex=25                     Discard biexonic transcripts having one exon size lower to this value (default 25)

  * Overlapping specification:
      -f,--minfrac_over=0           minimal fraction out of the candidate lncRNA size to be considered for overlap [default 0 i.e 1nt]
      -p,--proc=4                   number of thread for computing overlap [default 4]

  * Log output:
      -o,--outlog=file.log          Specify the log file of output which [default infile.log]

```


### 2- FEELnc_codpot.pl

The second step of the pipeline (FEELnc_codpot) aims at computing the CPS i.e the coding potential score (between [0-1]) foreach of the candidate transcripts in the candidate_lncRNA.gtf file.

**- INPUT :**

It makes use of the intrinsic properties of input sequences (ORF coverage and mRNA sizes, k-mer frequencies...) based on 2 training files:


	- known_mRNA.gtf (or .fa)   : a set of known protein_coding transcripts
	- known_lncRNA.gtf  (or .fa): a set of known lncRNA transcripts

If you have a set of known lncRNAs, you could run the module like:

	FEELnc_codpot.pl -i candidate_lncRNA.gtf -a known_mRNA.gtf -l known_lncRNA.gtf

However, for most organisms, the set of known_lncRNA transcripts is not known and thus
a set of genomic intergenic regions are automatically extracted as the lncRNA training set.
In this case, the reference genome file is required (ref_genome.FA)

    FEELnc_codpot.pl -i candidate_lncRNA.gtf -a known_mRNA.gtf -g ref_genome.FA


As in the previous module, if your reference annotation file  ("*ref_annotation.GTF*") contains additionnal fields such **transcript_biotype** and/or **transcript_status** in the [GENCODE annotation](http://www.gencodegenes.org/gencodeformat.html) or [ENSEMBL](http://www.ensembl.org), you can extract them manually or by using the **-b option** (as  to get the best training set of known mRNAs.

    FEELnc_codpot.pl -i candidate_lncRNA.gtf -a ref_annotation.GTF \
    -g ref_genome.FA \
    -b transcript_biotype=protein_coding -b transcript_status=KNOWN



To calculate the CPS cutoff separating coding (mRNAs) versus long non-coding RNAs (lncRNAs),
FEELnc_codpot uses a R script that will make a 10 fold cross-validation on the input training files and finally, extracts the CPS that maximizes sensitivity (Sn) and Specificity (Sp) (thanks to the ROCR library)


**- OUTPUT :**

If your input file is called **INPUT**, this second module will create these output files:

	 - {INPUT}_RF_learningData.txt: FEELnc metrics scores (ORF coverage and mRNA sizes, k-mer frequencies and labels) for the training files.
	 - {INPUT}_RF_statsLearn_CrossValidation.txt: statistics for n cross validation on the training files.
	 - {INPUT}_RF_TGROC.png: TwoGraph ROC curve plot to select the best coding potential cutoff.
	 - {INPUT}_RF.txt: FEELnc metrics scores (ORF coverage and mRNA sizes, k-mer frequencies and labels) for the testing file.
	 - {INPUT}_RF_varImpPlot.png: Dotchart plot of variable importance as measured by a Random Forest.

	 - {INPUT}.lncRNA.gtf || {INPUT}.lncRNA.fa: a .GTF/.FA file of the transcripts below the CPS (i.e the final set of lncRNAs).
	 - {INPUT}.mRNA.gtf || {INPUT}.mRNA.fa: a .GTF/.FA file of the transcripts above the coding potential cutoff (i.e the final set of mRNAs).
	 - [Possibly] {INPUT}.noORF.gtf || {INPUT}.noORF.fa: a .GTF/.FA file of the transcripts without any ORF found by FEELnc using the specified --testorftype option (see FEELnc_codpot.pl options description for more details). Transcripts contained in this file most probably correspond to lncRNAs.

An example of an {INPUT}_RF_TGROC.png graphic obtained using automatic threshold, i.e. sensibility equal to specificity on 10-fold cross-validation:
![ScreenShot](./image/FEELnc_codpot_performance.png)

**- FULL OPTIONS (FEELnc_codpot.pl --help) :**

```

Usage:
    FEELnc_codpot.pl -i transcripts.GTF -a known_mRNA.GTF -g genome.FA -l
    known_lnc.GTF [options...]

Options:
  General:
      --help                Print this help
      --man                 Open man page
      --verbosity           Level of verbosity

  Mandatory arguments:
      -i,--infile=file.gtf/.fasta           Specify the .GTF or .FASTA file  (such as a cufflinks transcripts/merged .GTF or .FASTA file)
      -a,--mRNAfile=file.gtf/.fasta         Specify the annotation .GTF or .FASTA file  (file of protein coding transcripts .GTF or .FASTA file)

  Optional arguments:
      -g,--genome=genome.fa                 Genome file or directory with chr files (mandatory if input is .GTF) [ default undef ]
      -l,--lncRNAfile=file.gtf/.fasta       Specify a known set of lncRNA for training .GTF or .FASTA  [ default undef ]
      -b,--biotype                          Only consider transcripts having this(these) biotype(s) from the reference annotation (e.g : -b transcript_biotype=protein_coding,pseudogene) [default undef i.e all transcripts]
      -n,--numtx=undef                      Number of transcripts required for the training and for all transcripts in the annotation use undef [ default undef ]
      -r,--rfcut=[0-1]                      Random forest voting cutoff [ default undef i.e will compute best cutoff ]
      --spethres=undef                      Two specificity threshold based on the 10-fold cross-validation, first one for mRNA and the second for lncRNA, need to be in ]0,1[ on separated by a ','
      -k,--kmer="3,6,9"                     Kmer size list with sizes separated by ',' as string [ default "3,6,9" ], the maximum value for one size is '15'
      -o,--outname="./"                     Output filename [ default infile_name ]
      --outdir="./"                         Output directory [ default current directory ]
      -s,--sizeinter=0.75                   Ratio between mRNA sequence and non-coding intergenic extracted region sizes [default 0.75 ]
      --learnorftype=1                      Integer [0,1,2,3,4] to specify the type of longest ORF computation [ default: 1 ] for learning data set. If the CDS is annotated in the .GTF, then the CDS is considered as the longest ORF, whatever the --orftype value.
                                                    '0': ORF with start and stop codon;
                                                    '1': same as '0' and ORF with only a start codon, take the longest;
                                                    '2': same as '0' and ORF with only a stop codon,  take the longest;
                                                    '3': same as '0' and ORF with a start or a stop,  take the longest (see '1' and '2');
                                                    '4': same as '3' but if no ORF is found, take the input sequence as ORF.
      --testorftype=1                       Integer [0,1,2,3,4] to specify the type of longest ORF calculate [ default: 1 ] for test data set. See --learnortype description for more informations.
      --ntree                               Number of trees used in random forest [ default 500 ]

  Debug arguments:
      --keeptmp                           To keep the temporary files in a 'tmp' directory the outdir, by default don't keep it (0 value). Any other value than 0 will keep the temporary files
      -v,--verbosity=0                         Which level of information that need to be print [ default 0 ]
      --seed=1234                           Used to fixe the seed value for the extraction of intergenic DNA region to get lncRNA like sequences and for the random forest [ default 1234 ]

  Intergenic lncrna extraction:
            -to be added

```

### 3- FEELnc_classifier.pl

The last step of the pipeline consists in classifying new lncRNAs w.r.t to the localisation and the direction of transcription of proximal RNA transcripts

Indeed, classifying lncRNAs with mRNA could help to predict functions for lncRNAs.
For all newly identified lncRNAs transcripts, a sliding window strategy is used to check for possible overlap with transcripts in the reference annotation.

If an overlap is found, the lncRNAs is considered as **GENIC** otherwise it is **INTERGENIC** (lincRNA) .

Then, subclasses are defined according the features and direction of overlap (See OUTPUT for full details).

Foreach lncRNA interaction, a best lncRNA:mRNA interaction is identified in the output file by a line starting with '*' and defined as followed:

 - for **INTERGNIC** : the best RNA partner is the closest to the lincRNA
 - for **GENIC**	: the best RNA partner is by rule of priority (exonic then the fraction of exonicof overlap) then intronic then containing.


```
	FEELnc_classifier.pl -i lncRNA.gtf -a  ref_annotation.GTF > lncRNA_classes.txt
```


**- OUTPUT :**


A summary of the number of the input parameters and numbers of interactions is given in the beginning of the file.

The classes are defined as in Derrien et al, Genome Research. 2012, and can be prioritized according to :

- **Intergenic lncRNAs** (i.e lincRNAs)
 - *divergent*  : when the lincRNA is transcribed in an opposite direction (head to head) w.r.t to the closest RNA partner (depending on distance, they could share a bi-directional promoter).
 - *convergent*: when the lincRNA is transcribed in a convergent direction w.r.t to the closest RNA partner.
 - *same_strand*: when the lincRNA is transcribed in a same starnd w.r.t to the closest RNA partner

- **Genic lncRNAs** (lncRNAs overlapping an annotated RNA)
 - *Exonic* :
    - antisense : at least one lncRNA exon overlaps in antisense an RNA exon
    - sense :  should correspond to lncRNAs overlapping non protein-coding transcripts (depending on the Filter option) and most probably lncRNAs host transcripts for small ncRNAs (snoRNAs, snRNAs, miRNAs...)
 - *Intronic* :
    - antisense : lncRNA exon overlaps in antisense RNA introns (but none exons)
    - sense : lncRNA exon overlaps in sense RNA introns (but none exons)
 - *containing*:
    - antisense : lncRNA intron overlaps antisense RNA
    - sense : lncRNA intron overlaps sense RNA exons

Illustration of the classification:
![ScreenShot](./image/FEELnc_lncRNA_classification_intergenic.png)
![ScreenShot](./image/FEELnc_lncRNA_classification_genic.png)

Example:

```
#FEELnc Classification
#lncRNA file :  lncrna : lncRNA_34.gtf
#mRNA file : 34.gtf
#Minimal window size : 1000
#Maximal window size : 1000
#Number of lncRNA : 462
#Number of mRNA : 374
#Number of interaction : 180
#Number of lncRNA without interaction : 369
#List of lncRNA without interaction : TCONS_00162422 TCONS_00163593 TCONS_00160948 TCONS_00162636 TCONS_00161917 TCONS_00161884 TCONS_00162654 TCONS_00164153 TCONS_00162630 TCONS_00161117 TCONS_00163030 TCONS_00162985 TCONS_00162306 TCONS_00162438
#INTERACTIONS
*lncRNA  XLOC_053965  TCONS_00160885  RNA_partner  ENSCAFG00000025966  ENSCAFT00000040249  antisense  genic       0    Status=containing   Subtype=intronic
lncRNA   XLOC_053965  TCONS_00160885  RNA_partner  ENSCAFG00000008990  ENSCAFT00000014276  antisense  intergenic  217  Status=divergent    Subtype=upstream
lncRNA   XLOC_053965  TCONS_00160885  RNA_partner  ENSCAFG00000008997  ENSCAFT00000014283  antisense  intergenic  979  Status=convergent   Subtype=downstream
lncRNA   XLOC_053965  TCONS_00160885  RNA_partner  ENSCAFG00000008997  ENSCAFT00000014282  antisense  intergenic  979  Status=convergent   Subtype=downstream
*lncRNA  XLOC_055197  TCONS_00163966  RNA_partner  ENSCAFG00000014726  ENSCAFT00000046107  antisense  genic       0    Status=nested       Subtype=exonic
lncRNA   XLOC_055197  TCONS_00163966  RNA_partner  ENSCAFG00000014726  ENSCAFT00000023388  antisense  genic       0    Status=overlapping  Subtype=exonic
*lncRNA  XLOC_054505  TCONS_00162285  RNA_partner  ENSCAFG00000022834  ENSCAFT00000034941  sense      genic       0    Status=containing   Subtype=intronic
lncRNA   XLOC_054505  TCONS_00162285  RNA_partner  ENSCAFG00000014605  ENSCAFT00000023189  antisense  genic       0    Status=containing   Subtype=intronic
*lncRNA  XLOC_054288  TCONS_00161746  RNA_partner  ENSCAFG00000013503  ENSCAFT00000043054  antisense  genic       0    Status=overlapping  Subtype=intronic
lncRNA   XLOC_054288  TCONS_00161746  RNA_partner  ENSCAFG00000032221  ENSCAFT00000046514  sense      genic       0    Status=containing   Subtype=intronic
```

Here is showed 4 interactions concerning 3 lncRNAs (TCONS_00160885, TCONS_00161746, TCONS_00162285 and TCONS_00163966) where one lncRNA (TTCONS_00160885) has 4 interactions with a window size of 1,000 nt.
(The best interactions are marked as ***lncRNA**)

\* **Note1**: At the moment, the interactions are computed with the reference file (-a option).
Therefore, the possibly newly identified mRNAs in the previous step are not included by default (but you could include them by (cuff)merging with you reference annotation).

\* **Note2**:  you may see a warning message like this:

	lncRNA_ID and ENSXXXX are overlapping in the same strand
Depending on your filtering options, this may correspond to a non-protein-coding transcript (pseudogene, miRNA) which overlaps the lncRNA

**- FULL OPTIONS (FEELnc_classifier.pl --help) :**

```
    General:
      --help                Print this help
      --man                 Open man page
      --verbosity           Level of verbosity

    Mandatory arguments:
      -i,--lncrna=file.(gtf/gff) Specify the lncRNA GTF/GFF file
      -a,--mrna=file.(gtf/gff)    Specify the annotation GTF/GFF file (file of protein coding annotaion)

    Filtering arguments:
      -w,--window=200               Size of the window around the lncRNA to compute interactins/classification [default 10000]
      -m, --maxwindow=10000 Maximal size of the window during the expansion process [default 10000]
```

## Warnings

 - In the installation (and/or) the classifier step, you may see a warning like
```
Can't call method "close" on an undefined value Bio/DB/SeqFeature/Store/berkeleydb.pm
```

## Authors

 - Valentin Wucher
 - Fabrice Legeai
 - Thomas Derrien

## Acknowledgments

