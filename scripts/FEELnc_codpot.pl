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
my $verbosity  = 0;
# my $outputlog;
my $numtx    = 3000;	# number of tx for training
my $minnumtx = 100;	# Min number of tx for training (a too small value will result in a bad regression)


# VW Add a variable to get the kmer size which are used to calculat the kmer scores
my $kmerList = '1,2,3,4,5,6';

# VW Add a variable to keep tmp file, default don't keep
my $keepTmp = 0;

# VW If random forest (rf/RF) cutoff is defined, no need to compute it on TP lncRNA and mRNA
my $rfcut = undef;

# VW Add option to select the calculate orf for learning and test data sets
my $orfTypeLearn = 0;
my $orfTypeTest  = 3;

# VW Add an option to specify the output directory, default current directory
my $outDir = "./";

# Intergenic extraction:
my $maxTries   = 10;
my $maxN       = 5;
my $sizecorrec = 1; # a float value between 0 and 1

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions(
    'i|infile=s'     => \$infile,
    'a|mRNAfile=s'   => \$mRNAfile,
    'l|lncRNAfile=s' => \$lncRNAfile,
    'g|genome=s'     => \$genome,
    'n|numtx=i'      => \$numtx,
    'b|biotype=s'    => \%biotype,
    'r|rfcut=f'      => \$rfcut,
    'k|kmer=s'       => \$kmerList,
    's|sizeinter=f'  => \$sizecorrec,
    'learnorftype=i' => \$orfTypeLearn,
    'testorftype=i'  => \$orfTypeTest,
    'o|outdir=s'     => \$outDir,
    'keeptmp'        => \$keepTmp,
    'v|verbosity=i'  => \$verbosity,
    'help|?'         => \$help,
    'man'            => \$man
    # 	"o|outlog=s"     => \$outputlog,
    ) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;


# Test parameters
pod2usage("Error: Cannot read your input GTF file '$infile'...\nFor help, see:\n$progname --help\n") unless( -r $infile);
pod2usage("Error: Cannot read your input annotation file '$mRNAfile'...\nFor help, see:\n$progname --help\n") unless( -r $mRNAfile);
pod2usage ("- Error: \$numtx option (number of transcripts for training) '$numtx' should be greater than $minnumtx  \n") unless ($numtx >= $minnumtx);
if (defined $rfcut){
    pod2usage ("- Error: \$rfcut option '$rfcut' should be a float between 0 and 1 [0-1] \n") unless ($rfcut >= 0 and $rfcut <= 1);
}
pod2usage ("- Error: \$sizecorrec option (ratio between mRNAs sequence lenghts and intergenic non coding sequence lenghts) '$sizecorrec' should be a float between 0 and 1 [0-1] \n") unless ($sizecorrec >= 0 and $sizecorrec <= 1);
pod2usage ("- Error: \$orfTypeLearn option '$orfTypeLearn' should be equal to 0, 1, 2, 3 or 4 (see 'FEELnc_codpot.pl --help' for more information) \n") unless ($orfTypeLearn==0 || $orfTypeLearn==1 || $orfTypeLearn==2 || $orfTypeLearn==3 || $orfTypeLearn==4);
pod2usage ("- Error: \$orfTypeTest option '$orfTypeTest' should be equal to 0, 1, 2, 3 or 4 (see 'FEELnc_codpot.pl --help' for more information) \n") unless ($orfTypeTest==0 || $orfTypeTest==1 || $orfTypeTest==2 || $orfTypeTest==3 || $orfTypeTest==4);
pod2usage ("- Error: \$outDir option '$outDir' is not a directory or it does not exist \n") unless (-d $outDir);


# Check the max kmersize
my @kmerTable = split(/,/,$kmerList);
my $kmerMax   = max @kmerTable;
pod2usage ("- Error: \$kmerList option '$kmerList' is not valid. One of the size is stricly greater than '15' (see --help for more details). \n") unless ($kmerMax <= 15);



# For $outDiradd a '/' at the end of the path
$outDir = $outDir."/";

# Create the directory for temporary files
my $outTmp = "/tmp/";
if($keepTmp!=0)
{
    $outTmp = $outDir."tmp/";
    mkdir $outTmp;
}




#############################################################

# test path
die "Error: You should set the environnment variable FEELNCPATH to the dir of installation\nexport FEELNCPATH=my_dir_of_install/\n(See README)\n" unless (defined $ENV{'FEELNCPATH'});
# my $rprogpath   = $ENV{'FEELNCPATH'}."/bin/".$rprog;
# pod2usage("Error: Cannot access FEELnc bin dir with path '$rprogpath'...\nCheck the environnment variable FEELNCPATH\n") unless( -r $rprogpath);
# my $pathRscript = Utils::pathProg("Rscript");
#VW: don't need cpat anymore my $pathlogit   = Utils::pathProg("cpat.py");
# test PYTHONPATH from CPAT : http://dldcc-web.brc.bcm.edu/lilab/liguow/CGI/cpat/_build/html/index.html#installation
# die "Error: You should set the PYTHONPATH env. variable to CPAT installation
# export PYTHONPATH=/home/user/CPAT/usr/local/lib/python2.7/site-packages:\$PYTHONPATH. #setup PYTHONPATH
# (See http://dldcc-web.brc.bcm.edu/lilab/liguow/CGI/cpat/_build/html/index.html#installation)\n" unless (defined $ENV{'PYTHONPATH'});

# Log File
##########
# my $commandline = qx/ps -o args $$/;
# if (!defined $outputlog){
# 	$outputlog	=	Utils::renamefile($infile, ".feelnccodpot.log");
# }
# open(LOG,">$outputlog") or die("Cannot open '$outputlog'");
#
#
# print LOG $commandline;
# print STDERR "> Results will be available in file: '$outputlog'\n";


# Die if lnc training file is not set and mRNA file is in FASTA: no possibility of intergenic extraction
my $mRNAfileformat = Utils::guess_format($mRNAfile);
pod2usage ("- Error: Cannot train the program if lncRNA training file (-l option) is not defined and mRNA file (-a option) is in FASTA format!\nPlease, provide the mRNA/annotation file in .GTF format so that I could extract intergenic sequences for training...\n") if (!defined $lncRNAfile && $mRNAfileformat eq "fasta");


warn "> Preparing files for random forest...\n";

# Define training file names
my $codFile    = $outTmp.Utils::renamefile($mRNAfile, ".codingTrain.fa");
my $codOrfFile = $outTmp.Utils::renamefile($mRNAfile, ".codingOrfTrain.fa");
my $nonFile;
my $nonOrfFile;
if(defined $lncRNAfile)
{
    $nonFile    = $outTmp.Utils::renamefile($lncRNAfile, ".nonCodingTrain.fa");
    $nonOrfFile = $outTmp.Utils::renamefile($lncRNAfile, ".nonCodingOrfTrain.fa");
}
else
{
    $nonFile    = $outTmp.Utils::renamefile($mRNAfile, ".nonCodingTrain.fa");
    $nonOrfFile = $outTmp.Utils::renamefile($mRNAfile, ".nonCodingOrfTrain.fa");
}

# Define test file names
my $testFile    = $outTmp.Utils::renamefile($infile, ".test.fa");
my $testOrfFile = $outTmp.Utils::renamefile($infile, ".testOrf.fa");

# Define output name
my $rfout = $outDir.basename($infile)."_RF.txt";


##########################################################
# mRNA file
#######
# add a refhash that will contain the mRNA ID that passed cDNA and ORF steps
# Will be used to checked for randomization
my $ref_cDNA_passed;

# Get cDNA and ORF for coding training file
if($mRNAfileformat eq "gtf")      # -- if GTF
{
    $ref_cDNA_passed = ExtractCdnaOrf::CreateORFcDNAFromGTF($mRNAfile, $codFile, $codOrfFile, $numtx, $minnumtx, $genome, 'exon,CDS,stop_codon,start_codon', \%biotype, $orfTypeLearn, $verbosity, $kmerMax);
}
elsif($mRNAfileformat eq "fasta") # -- if FASTA
{
    ExtractCdnaOrf::CreateORFcDNAFromFASTA($mRNAfile, $codFile, $codOrfFile, $numtx, $minnumtx, $orfTypeLearn, $verbosity, $kmerMax);
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
	$ref_cDNA_passed = ExtractCdnaOrf::CreateORFcDNAFromGTF($lncRNAfile, $nonFile, $nonOrfFile, $numtx, $minnumtx, $genome, 'exon', undef, $orfTypeLearn, $verbosity, $kmerMax);
    }
    elsif($lncRNAfileformat eq "fasta")
    {
	ExtractCdnaOrf::CreateORFcDNAFromFASTA($lncRNAfile, $nonFile, $nonOrfFile, $numtx, $minnumtx, $orfTypeLearn, $verbosity, $kmerMax);
    }
    else
    {
	die "Error: Unrecognized format for lncRNA training file '$lncRNAfile'\n";
    }
}
else                    # -- if lncRNA training file not defined
{
    # To get mRNA annotation
    my $refmrna = Parser::parseGTF($mRNAfile, 'exon,CDS,stop_codon,start_codon', undef , \%biotype , $verbosity);
    # Relocated mRNA sequence in intergenic regions to be used as a training lncRNA file
    print STDERR "> The lncRNA training file is not set...will extract intergenic region for training (can take a while...)\n";
    ExtractCdnaOrf::randomizedGTFtoFASTA($refmrna, $ref_cDNA_passed, $nonFile, $nonOrfFile, $genome, $numtx, $minnumtx, $sizecorrec, $orfTypeLearn, $maxTries, $maxN, $verbosity, $kmerMax);
}


##########################################################
# test file
#######
# Get cDNA and ORF for test file
if(Utils::guess_format($infile) eq "gtf")      # -- if GTF
{
    # Use undef for $nbtx to get all sequences and ORFs
    ExtractCdnaOrf::CreateORFcDNAFromGTF($infile, $testFile, $testOrfFile, undef, $minnumtx, $genome, 'exon', undef, $orfTypeTest, $verbosity, $kmerMax);
}
elsif(Utils::guess_format($infile) eq "fasta") # -- if FASTA
{
    ExtractCdnaOrf::CreateORFcDNAFromFASTA($infile, $testFile, $testOrfFile, undef, $minnumtx, $orfTypeTest, $verbosity, $kmerMax);
}
else
{
    die "Error: Unrecognized format for input file '$infile'...\n";
}


# # VW modif crade !
# # besoin des ORF pour lnc et test
# my $lncOrfFile  = $outTmp."lncRNA_ORF.fa";
# my $testOrfFile = $outTmp."test_ORF.fa";


# # VW : Récupère les ORF du jeu lncRNA et test, crade !!!!
# if (Utils::guess_format($lncfile) eq "gtf")
# {
#     &CreateORFcDNAFromGTF($lncfile, "/tmp/poubelle1", $lncOrfFile, $numtx, $genome, $orfTypeLearn, $verbosity);

# }
# else
# {
#     &CreateORFcDNAFromFASTA($lncfile, "/tmp/poubelle1", $lncOrfFile, $numtx, $orfTypeLearn, $verbosity);
# }

# # VW : utilise undef pour avoir l'ensemble des ORF
# if (Utils::guess_format($infile) eq "gtf")
# {
#     &CreateORFcDNAFromGTF($refin, "/tmp/poubelle2",  $testOrfFile, undef, $genome, $orfTypeTest, $verbosity);
# }
# elsif (Utils::guess_format($infile) eq "fasta")
# {
#     &CreateORFcDNAFromFASTA($infile, "/tmp/poubelle2",  $testOrfFile, undef, $orfTypeTest, $verbosity);
# }




#################################
# Launch RF on $infile in fasta

print STDERR "> Run random Forest on '$testFile':\n";
RandomForest::runRF($codFile, $codOrfFile, $nonFile, $nonOrfFile, $testFile, $testOrfFile, $rfout, $kmerList, $rfcut, $outDir, $verbosity, $keepTmp);

# Parse RF result
RandomForest::rfPredToOut($infile, $rfout, $outDir);


# Cleaning temporary files
if($keepTmp==0)
{
    unlink $codFile;
    unlink $codOrfFile;
    unlink $nonFile;
    unlink $nonOrfFile;
    unlink $testFile;
    unlink $testOrfFile;
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

  --help                Print this help
  --man                 Open man page
  --verbosity		Level of verbosity


=head2 Mandatory arguments

  -i,--infile=file.gtf/.fasta		Specify the .GTF or .FASTA file  (such as a cufflinks transcripts/merged .GTF or .FASTA file)
  -a,--mRNAfile=file.gtf/.fasta		Specify the annotation .GTF or .FASTA file  (file of protein coding transcripts .GTF or .FASTA file)


=head2 Optional arguments

  -g,--genome=genome.fa			Genome file or directory with chr files (mandatory if input is .GTF) [ default undef ]
  -l,--lncRNAfile=file.gtf/.fasta	Specify a known set of lncRNA for training .GTF or .FASTA  [ default undef ]
  -b,--biotype				Only consider transcripts having this(these) biotype(s) from the reference annotation (e.g : -b transcript_biotype=protein_coding,pseudogene) [default undef i.e all transcripts]
  -n,--numtx=2000			Number of transcripts required for the training [ default 2000 ]
  -r,--rfcut=[0-1]			Random forest voting cutoff [ default undef i.e will compute best cutoff ]
  -k,--kmer="2,3,4,5,6"			Kmer size list with size separate by ',' as string [ default "2,3,4,5,6" ], the maximum value for the size is '15'
  -o,--outdir="./"			Output directory [ default current directory ]
  -s,--sizeinter=0.75			Ratio between mRNA sequence lengths and non coding intergenic region sequence lengths as, by default, ncInter = mRNA * 0.75
  --learnorftype=0			Integer [0,1,2,3,4] to specify the type of longest ORF calculate (default: 0) for learning data set.
					If the CDS is annotated in the .GTF, then the CDS is considered as the longest ORF, whatever the --orftype value.
						'0': only ORF with start and stop codon;
						'1': same as '0' and if no ORF found, take the longest with a start codon;
						'2': same as '1' but with a stop codon;
						'3': same as '0' and if no ORF found, take the longest between ORF with start or stop codon (see '1' and '2');
						'4': same as '3' but if no ORF is found, take the input sequence as ORF.
  --testorftype=3			Integer [0,1,2,3,4] to specify the type of longest ORF calculate (default: 3) for test data set. See --learnortype description for more informations.
  --keeptmp=0				To keep the temporary files in a 'tmp' directory the outdir, by default don't keep it (0 value). Any other value than 0 will keep the temporary files


=head2 Intergenic lncRNA extraction

	-to be added


=head2 Log output

  -o,--outlog=file.log		Specify the log file of output which [default infile.log]



=head1 AUTHORS

=over 4

=item *

Thomas DERRIEN <tderrien@univ-rennes1.fr>

=item *

Fabrice LEGEAI <fabrice.legeai@inria.fr>

=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2014 by IGDR - CNRS

=cut
