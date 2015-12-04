#!/usr/bin/perl -w

#
# Modification by V.Wucher april 16 2015:
# 	Modification of the predicting method: use now random forest (package R randomForest)
#

# Perl libs
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Data::Dumper;
use List::Util qw( min max );

# lib directory: {FEELnc github directory}/lib/
use Parser;
use ExtractFromHash;
use ExtractFromFeature;
use Intersect;
use Utils;
use Orf;
use RandomForest;
use ExtractCdnaOrf;

# Program name
my $progname = basename($0);

# Variables
my $infile     = '';
my $mRNAfile   = '';
my $genome     = undef;
my $lncRNAfile = undef;
my %biotype;
my $man        = 0;
my $help       = 0;
my $verbosity  = 1;
# my $outputlog;
my $numtx    = undef; # number of mRNAs and lncRNAs tx for training separate by a ','. undef or '.'  for all transcripts
my $minnumtx = 100;   # Min number of tx for training (a too small value will result in a bad learning)


# VW Add a variable to get the kmer size which are used to get the kmer scores
my $kmerList = '1,2,3,6,9,12';

# VW Add a variable to keep tmp file, default don't keep
my $keepTmp = 0;

# VW If random forest (rf/RF) cutoff is defined, no need to compute it on TP lncRNA and mRNA
#    and a $speThres for the thresolds on mRNA and lncRNA specificity
my $rfcut        = undef;
my $speThres     = undef;
my @speThresList = undef;

# VW Add option to select the calculate orf for learning and test data sets
my $orfTypeLearn = 1;
my $orfTypeTest  = 1;

# VW Add an option to specify the output directory, default current directoryand an out name
my $outDir  = "./feelnc_codpot_out";
my $outName = "";

# VW Add an option to change the number of trees use in random forest
my $nTree = 500;

# VW Add an option to fixe the seed
my $seed = 1234;

# VW Add nbr proc
my $proc = 1;

# VW Add a percentage to get two learning file
my $perc = 0.1;

# Intergenic extraction:
my $maxTries   = 10;
my $maxN       = 5;
my $sizecorrec = 0.75; # a float value between 0 and 1

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions(
    'i|infile=s'     => \$infile,
    'a|mRNAfile=s'   => \$mRNAfile,
    'l|lncRNAfile=s' => \$lncRNAfile,
    'g|genome=s'     => \$genome,
    'n|numtx=s'      => \$numtx,
    'b|biotype=s'    => \%biotype,
    'r|rfcut=f'      => \$rfcut,
    'spethres=s'     => \$speThres,
    'k|kmer=s'       => \$kmerList,
    's|sizeinter=f'  => \$sizecorrec,
    'learnorftype=i' => \$orfTypeLearn,
    'testorftype=i'  => \$orfTypeTest,
    'ntree=i'        => \$nTree,
    'outdir=s'       => \$outDir,
    'o|outname=s'    => \$outName,
    'percentage=f'   => \$perc,
    'keeptmp'        => \$keepTmp,
    'v|verbosity=i'  => \$verbosity,
    'p|processor=i'  => \$proc,
    'seed=i'         => \$seed,
    'help|?'         => \$help,
    'man'            => \$man
    ) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;


# Test parameters
pod2usage("- Error: Cannot read your input GTF file '$infile'...\nFor help, see:\n$progname --help\n") unless( -r $infile);
pod2usage("- Error: Cannot read your input annotation file '$mRNAfile'...\nFor help, see:\n$progname --help\n") unless( -r $mRNAfile);
if (defined $rfcut){
    pod2usage ("- Error: --rfcut option '$rfcut' should be a float between 0 and 1 [0-1] \n") unless ($rfcut >= 0 and $rfcut <= 1);
}
pod2usage ("- Error: --sizecorrec option (ratio between mRNAs sequence lenghts and intergenic non coding sequence lenghts) '$sizecorrec' should be a float between 0 and 1 [0-1]\n") unless ($sizecorrec >= 0 and $sizecorrec <= 1);
pod2usage ("- Error: --orfTypeLearn option '$orfTypeLearn' should be equal to 0, 1, 2, 3 or 4 (see 'FEELnc_codpot.pl --help' for more information)\n") unless ($orfTypeLearn==0 || $orfTypeLearn==1 || $orfTypeLearn==2 || $orfTypeLearn==3 || $orfTypeLearn==4);
pod2usage ("- Error: --orfTypeTest option '$orfTypeTest' should be equal to 0, 1, 2, 3 or 4 (see 'FEELnc_codpot.pl --help' for more information)\n") unless ($orfTypeTest==0 || $orfTypeTest==1 || $orfTypeTest==2 || $orfTypeTest==3 || $orfTypeTest==4);
# pod2usage ("- Error: --outDir option '$outDir' is not a directory or it does not exist \n") unless (-d $outDir);
pod2usage ("- Error: --nTree option '$nTree' should be strictly positive\n") unless ($nTree > 0);
pod2usage ("- Error: --rfcut and --spethres specified, only one of the two options can be used (default one threshold defined on a 10-fold cross-validation)\n") if((defined $rfcut) && (defined $speThres));
pod2usage ("- Error: -p/--processor option '$proc' should be a positive integer\n") unless ($proc >= 1);
pod2usage ("- Error: --percentage option '$perc' should be a number in ]0;1[\n") unless ($perc>0 && $perc<1);


# Check the max kmersize
my @kmerTable = split(/,/,$kmerList);
my $kmerMax   = max @kmerTable;
pod2usage ("- Error: \$kmerList option '$kmerList' is not valid. One of the size is stricly greater than '15' (see --help for more details)\n") unless ($kmerMax <= 15);

# Check threshold values for the mRNAs and lncRNAs specificity
if((defined $speThres))
{
    @speThresList = split(/,/, $speThres);
    if(@speThresList!=2)
    {
	pod2usage ("- Error: --speThres option '$speThres' should be a list of two value separated by a ',' (see --help for more details)\n");
    }
    if(($speThresList[0]<=0) || ($speThresList[0]>=1) || ($speThresList[1]<=0) || ($speThresList[1]>=1))
    {
	pod2usage ("- Error: one value of --speThres option '$speThres' is equal or greater than 1 or equal or lesser than 0, should be in ]0,1[ (see --help for more details)\n");
    }
}



# Check the number of mRNAs ($numtxCod) and lncRNAs ($numtxNon) transcripts use for learning
my $numtxCod = undef;
my $numtxNon = undef;

if(defined $numtx)
{
    if(($numtx =~ tr/,//) != 1)
    {
	pod2usage ("- Error: --numtx option '$numtx' should be a list of two value separated by a ',' (see --help for more details)\n");
    }
    else
    {
	($numtxCod,$numtxNon) = split(/,/, $numtx);
    }

    if($numtxCod =~ /^\d+$/ &&  $numtxCod < $minnumtx)
    {
	print $numtxCod."     ".$minnumtx."\n";
	pod2usage("- Error: number of mRNAs transcripts for training in --numtx option '$numtx' is not valid. Should be greater than $minnumtx or void to keep all the annotation (see --help for more details)\n");
    }
    elsif(($numtxCod !~ /^\d+$/) && $numtxCod ne "")
    {
	pod2usage("- Error: number of mRNAs transcripts for training in --numtx option '$numtx' is not valid. Should be greater than $minnumtx or void to keep all the annotation (see --help for more details)\n");
    }
    elsif($numtxCod eq "")
    {
	$numtxCod = undef;
    }

    if($numtxNon =~ /^\d+$/ &&  $numtxNon < $minnumtx)
    {
	print $numtxNon."     ".$minnumtx."\n";
	pod2usage("- Error: number of mRNAs transcripts for training in --numtx option '$numtx' is not valid. Should be greater than $minnumtx or void to keep all the annotation (see --help for more details)\n");
    }
    elsif(($numtxNon !~ /^\d+$/) && $numtxNon ne "")
    {
	pod2usage("- Error: number of mRNAs transcripts for training in --numtx option '$numtx' is not valid. Should be greater than $minnumtx or void to keep all the annotation (see --help for more details)\n");
    }
    elsif($numtxNon eq "")
    {
	$numtxNon = undef;
    }
}


# Default option for $outName
if($outName eq "")
{
    $outName=basename($infile);
}

# For $outDiradd a '/' at the end of the path
if (-d $outDir ){
	print STDERR "Warning: Output directory '$outDir' already exists... files might be overwritten!\n";
}elsif (-r $outDir){
	die "Error: Output directory '$outDir' is a file... aborting!\n";
} else{
	my $cmdline="mkdir -p $outDir";
	system($cmdline);
}
$outDir = $outDir."/"; # add "/" at the end in case it is forgotten in pasting outdir and outname


# Create the directory for temporary files with the job id ($$) and the temporary name
my $outTmp = "/tmp/";
if($keepTmp!=0)
{
    $outTmp = $outDir."/tmp/";
    mkdir $outTmp;
}
my $nameTmp = $outTmp."/".$$."_".$outName;


# Initialize the seed
srand($seed);

# If $numtx is undef, then learning on all transcripts, can be long so print a warning...
print STDERR "You do not have specified a maximum number mRNAs transcripts for the training. Use all the annotation, can be long...\n"  if(!defined $numtxCod);
print STDERR "You do not have specified a maximum number lncRNA transcripts for the training. Use all the annotation, can be long...\n" if(!defined $numtxNon);


#############################################################

# test path
die "Error: You should set the environnment variable FEELNCPATH to the dir of installation\nexport FEELNCPATH=my_dir_of_install/\n(See README)\n" unless (defined $ENV{'FEELNCPATH'});
# Rscript
my $rprogpath = $ENV{'FEELNCPATH'}."/utils/codpot_randomforest.r";
die "Error: The environnment variable FEELNCPATH does not reach the 'utils/codpot_randomforest.r' script\n" unless (-r $rprogpath);
# KIS path
my $kisPath = Utils::pathProg("KmerInShort");


# Die if lnc training file is not set and mRNA file is in FASTA: no possibility of intergenic extraction
my $mRNAfileformat = Utils::guess_format($mRNAfile);
pod2usage ("- Error: Cannot train the program if lncRNA training file (-l option) is not defined and mRNA file (-a option) is in FASTA format!\nPlease, provide the mRNA/annotation file in .GTF format so that I could extract intergenic sequences for training...\n") if (!defined $lncRNAfile && $mRNAfileformat eq "fasta");


# Define fasta file names
my $codFile    = $nameTmp.".coding_rna.fa";
my $codOrfFile = $nameTmp.".coding_orf.fa";
my $nonFile    = $nameTmp.".noncoding_rna.fa";
my $nonOrfFile = $nameTmp.".noncoding_orf.fa";
my $testFile    = $nameTmp.".test_rna.fa";
my $testOrfFile = $nameTmp.".test_orf.fa";

# Define output name
my $rfout = $outDir.$outName."_RF.txt";


##########################################################
# mRNA file
#######
# add a refhash that will contain the mRNA ID that passed cDNA and ORF steps
# Will be used to checked for randomization
my $ref_cDNA_passed;
my $refmrna;

# Get cDNA and ORF for coding training file
if($mRNAfileformat eq "gtf")      # -- if GTF
{

    print STDERR "> Extract ORFs/cDNAs for mRNAs from a GTF file\n";
    ($ref_cDNA_passed, $refmrna) = ExtractCdnaOrf::CreateORFcDNAFromGTF($mRNAfile, $codFile, $codOrfFile, $numtxCod, $minnumtx, $genome, 'exon,CDS,stop_codon,start_codon', \%biotype, $orfTypeLearn, $verbosity, $kmerMax);
}
elsif($mRNAfileformat eq "fasta") # -- if FASTA
{
    print STDERR "> Extract ORFs/cDNAs for mRNAs from a FASTA file\n";
    ExtractCdnaOrf::CreateORFcDNAFromFASTA($mRNAfile, $codFile, $codOrfFile, $numtxCod, $minnumtx, $orfTypeLearn, $verbosity, $kmerMax);
}
else
{
    die "Error : Unrecognized format for annotation file '$mRNAfile'\n";
}


##########################################################
# lncRNA file
#######
# Get cDNA and ORF for non coding training file
if(defined $lncRNAfile) # -- if file is defined, it means that we do not have to extract from intergenic
{
    my $lncRNAfileformat = Utils::guess_format($lncRNAfile);

    if ($lncRNAfileformat eq "gtf") # -- if GTF
    {
	print STDERR "> Extract ORFs/cDNAs for lncRNAs from a GTF file\n";
	($ref_cDNA_passed, $refmrna) = ExtractCdnaOrf::CreateORFcDNAFromGTF($lncRNAfile, $nonFile, $nonOrfFile, $numtxNon, $minnumtx, $genome, 'exon', undef, $orfTypeLearn, $verbosity, $kmerMax);
    }
    elsif($lncRNAfileformat eq "fasta")
    {
	print STDERR "> Extract ORFs/cDNAs for lncRNAs from a FASTA file\n";
	ExtractCdnaOrf::CreateORFcDNAFromFASTA($lncRNAfile, $nonFile, $nonOrfFile, $numtxNon, $minnumtx, $orfTypeLearn, $verbosity, $kmerMax);
    }
    else
    {
	die "Error: Unrecognized format for lncRNA training file '$lncRNAfile'\n";
    }
}
else # -- if lncRNA training file not defined
{
    # To get mRNA annotation
    # Relocated mRNA sequence in intergenic regions to be used as a training lncRNA file
    print STDERR "> The lncRNA training file is not set. Extract ORFs/cDNAs for lncRNAs from intergenic regions (can take a while)\n";
    ExtractCdnaOrf::randomizedGTFtoFASTA($refmrna, $ref_cDNA_passed, $nonFile, $nonOrfFile, $genome, $numtxNon, $minnumtx, $sizecorrec, $orfTypeLearn, $maxTries, $maxN, $verbosity, $kmerMax);
}


##########################################################
# test file
#######
# Get cDNA and ORF for test file
if(Utils::guess_format($infile) eq "gtf")      # -- if GTF
{
    print STDERR "> Extract ORFs/cDNAs for candidates RNAs from a GTF file\n";
    # Use undef for $nbtx/$numtx to get all sequences and ORFs
    ExtractCdnaOrf::CreateORFcDNAFromGTF($infile, $testFile, $testOrfFile, undef, undef, $genome, 'exon,CDS,stop_codon,start_codon', undef, $orfTypeTest, $verbosity, $kmerMax);
}
elsif(Utils::guess_format($infile) eq "fasta") # -- if FASTA
{
    print STDERR "> Extract ORFs/cDNAs for candidates RNAs from a FASTA file\n";
    # Use undef for $nbtx/$numtx to get all sequences and ORFs
    ExtractCdnaOrf::CreateORFcDNAFromFASTA($infile, $testFile, $testOrfFile, undef, undef, $orfTypeTest, $verbosity, $kmerMax);
}
else
{
    die "Error: Unrecognized format for input file '$infile'...\n";
}


#################################
# Divide each learning file into two file: one for the kmermodel and one for the random forest

my @codFileKmRf    = ($codFile.".forKmerModel.fa",    $codFile.".forRandomForest.fa");
my @codOrfFileKmRf = ($codOrfFile.".forKmerModel.fa", $codOrfFile.".forRandomForest.fa");
my @nonFileKmRf    = ($nonFile.".forKmerModel.fa",    $nonFile.".forRandomForest.fa");
my @nonOrfFileKmRf = ($nonOrfFile.".forKmerModel.fa", $nonOrfFile.".forRandomForest.fa");

Utils::divFasta($codFile,    $codFileKmRf[0],    $codFileKmRf[1],    $perc, $verbosity);
Utils::divFasta($codOrfFile, $codOrfFileKmRf[0], $codOrfFileKmRf[1], $perc, $verbosity);
Utils::divFasta($nonFile,    $nonFileKmRf[0],    $nonFileKmRf[1],    $perc, $verbosity);
Utils::divFasta($nonOrfFile, $nonOrfFileKmRf[0], $nonOrfFileKmRf[1], $perc, $verbosity);


#################################
# Launch RF on $infile in fasta

print STDERR "> Run random Forest on '$testFile'\n";
if(! defined $speThres)
{
    RandomForest::runRF(\@codFileKmRf, \@codOrfFileKmRf, \@nonFileKmRf, \@nonOrfFileKmRf, $testFile, $testOrfFile, $rfout, $kmerList, $rfcut, $nTree, $outDir, $verbosity, $nameTmp, $keepTmp, $seed, $proc);
}
else
{
    RandomForest::runRF(\@codFileKmRf, \@codOrfFileKmRf, \@nonFileKmRf, \@nonOrfFileKmRf, $testFile, $testOrfFile, $rfout, $kmerList, $speThres, $nTree, $outDir, $verbosity, $nameTmp, $keepTmp, $seed, $proc);
}

# Parse RF result
RandomForest::rfPredToOut($infile, $rfout, $outDir, $outName);


# Check if a TUCp file exists, if there is one so write this message
if( -e $outDir.$outName.".TUCp.gtf" )
{
    print STDERR "\nTranscripts of Unknown Coding Potential (TUCps) found: check file '".$outDir.$outName.".TUCp.gtf'.\n";
    print STDERR "\tThese transcripts correspond to transcripts with a coding potential between the mRNA and the lncRNA specificity threshold.\n";
}
if( -e $outDir.$outName.".TUCp.fa" )
{
    print STDERR "\nTranscripts of Unknown Coding Potential (TUCps) found: check file '".$outDir.$outName.".TUCp.fa'.\n";
    print STDERR "\tThese transcripts correspond to transcripts with a coding potential between the mRNA and the lncRNA specificity threshold.\n";
}

# Check if a noORF file exists, if there is one so write this message
if( -e $outDir.$outName.".noORF.gtf" )
{
    print STDERR "\nTranscripts without ORF found: check file '".$outDir.$outName.".noORF.gtf'.\n";
    print STDERR "\tThese transcripts correspond probably to lncRNAs since no ORF was found. You might also want to change the --testorftype option.\n";
}
if( -e $outDir.$outName.".noORF.fa" )
{
    print STDERR "\nTranscripts without ORF found: check file '".$outDir.$outName.".noORF.fa'.\n";
    print STDERR "\tThese transcripts correspond probably to lncRNAs since no ORF was found. You might also want to change the --testorftype option.\n";
}


# Cleaning temporary files
if($keepTmp==0)
{
    unlink $codFile;
    unlink $codOrfFile;
    unlink $nonFile;
    unlink $nonOrfFile;
    unlink $testFile;
    unlink $testOrfFile;
    unlink $codFileKmRf[0];
    unlink $codFileKmRf[1];
    unlink $codOrfFileKmRf[0];
    unlink $codOrfFileKmRf[1];
    unlink $nonFileKmRf[0];
    unlink $nonFileKmRf[1];
    unlink $nonOrfFileKmRf[0];
    unlink $nonOrfFileKmRf[1];
}


__END__

=pod

=encoding UTF-8

=head1 NAME

FEELnc_codpot.pl - Compute the coding potential of an candidate transcripts

=head1 VERSION

version 0.01

=head1 SYNOPSIS

FEELnc_codpot.pl -i transcripts.GTF -a known_mRNA.GTF -g genome.FA -l known_lnc.GTF  [options...]

=head1 DESCRIPTION

FEELnc (Fast and Effective Extraction of Long non-coding RNAs) is dedicated to the annotation of lncRNAs
based on a set of transcripts as input (basically a cufflink transcripts.gtf file)
The second step if the pipeline (FEELnc_codpot) aims at computing coding potential of the input transcripts.

=head1 OPTIONS

=head2 General

  --help				Print this help
  --man					Open man page
  --verbosity				Level of verbosity


=head2 Mandatory arguments

  -i,--infile=file.gtf/.fasta		Specify the .GTF or .FASTA file  (such as a cufflinks transcripts/merged .GTF or .FASTA file)
  -a,--mRNAfile=file.gtf/.fasta		Specify the annotation .GTF or .FASTA file  (file of protein coding transcripts .GTF or .FASTA file)


=head2 Optional arguments

  -g,--genome=genome.fa			Genome file or directory with chr files (mandatory if input is .GTF) [ default undef ]
  -l,--lncRNAfile=file.gtf/.fasta	Specify a known set of lncRNA for training .GTF or .FASTA  [ default undef ]
  -b,--biotype				Only consider transcripts having this(these) biotype(s) from the reference annotation (e.g : -b transcript_biotype=protein_coding,pseudogene) [default undef i.e all transcripts]
  -n,--numtx=undef			Number of mRNA and lncRNA transcripts required for the training. mRNAs and lncRNAs numbers need to be separate by a ',': i.e. 1500,1000 for 1500 mRNAs and 1000 lncRNAs. For all the annotation, let it blank [ default undef, all the two annotations ]
  -r,--rfcut=[0-1]			Random forest voting cutoff [ default undef i.e will compute best cutoff ]
  --spethres=undef			Two specificity threshold based on the 10-fold cross-validation, first one for mRNA and the second for lncRNA, need to be in ]0,1[ on separated by a ','
  -k,--kmer=1,2,3,6,9,12		Kmer size list with size separate by ',' as string [ default "1,2,3,6,9,12" ], the maximum value for one size is '15'
  -o,--outname={INFILENAME}		Output filename [ default infile_name ]
  --outdir="feelnc_codpot_out/"		Output directory [ default "./feelnc_codpot_out/" ]
  -s,--sizeinter=0.75			Ratio between mRNA sequence lengths and non coding intergenic region sequence lengths as, by default, ncInter = mRNA * 0.75
  --learnorftype=1			Integer [0,1,2,3,4] to specify the type of longest ORF calculate [ default: 1 ] for learning data set.
					If the CDS is annotated in the .GTF, then the CDS is considered as the longest ORF, whatever the --orftype value.
						'0': ORF with start and stop codon;
						'1': same as '0' and ORF with only a start codon, take the longest;
						'2': same as '1' but with a stop codon;
						'3': same as '0' and ORF with a start or a stop, take the longest (see '1' and '2');
						'4': same as '3' but if no ORF is found, take the input sequence as ORF.
  --testorftype=1			Integer [0,1,2,3,4] to specify the type of longest ORF calculate [ default: 1 ] for test data set. See --learnortype description for more informations.
  --ntree				Number of trees used in random forest [ default 500 ]
  --percentage=0.1			Percentage of the training file use for the training of the kmer model. What remains will be used to train the random forest

=head2 Debug arguments

  --keeptmp=0				To keep the temporary files in a 'tmp' directory the outdir, by default don't keep it (0 value). Any other value than 0 will keep the temporary files
  --verbosity=1				Integer [0,1,2]: which level of information that need to be print [ default 1 ]. Note that that printing is made on STDERR
  --seed=1234				Use to fixe the seed value for the extraction of intergenic DNA region to get lncRNA like sequences and for the random forest [ default 1234 ]


=head2 Intergenic lncRNA extraction

	-to be added



=head1 AUTHORS

=over 4


=item *

Valentin WUCHER <vwucher@univ-rennes1.fr>


=item *

Thomas DERRIEN <tderrien@univ-rennes1.fr>

=item *

Fabrice LEGEAI <fabrice.legeai@inria.fr>


=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2014 by IGDR - CNRS

=cut
