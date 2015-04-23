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



# lib directory : ~tderrien/bin/perl/lib/
use Parser;
use ExtractFromHash;
use ExtractFromFeature;
use Intersect;
use Utils;
use Orf;
#VW use Cpat;
use RandomForest;

# my $pathRcrossvalidation = "~tderrien/bin/perl/script/FEELnc/bin/crossValidation_cutoff.r";
my $rprog    = "10crossValidation_cutoff.r";
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


# VW Adding a variable to get the kmer size which are used to calculat the kmer scores
my $kmerList = '2,3,4,5,6';

# VW Adding a variable to keep tmp file, default don't keep
my $keepTmp = 0;

# If random forest (rf/RF) cutoff is defined, no need to compute it on TP lncRNA and mRNA
my $rfcut = undef;
my $sn_sp = undef;



# Intergenic extraction:
my $maxTries   = 10;
my $maxN       = 5;
my $sizecorrec = 1; # a float value between 0 and 1
#my $sizecorrec = 0.1; # VW

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
#############################################################

# test path
die "Error: You should set the environnment variable FEELNCPATH to the dir of installation\nexport FEELNCPATH=my_dir_of_install/\n(See README)\n" unless (defined $ENV{'FEELNCPATH'});
my $rprogpath   = $ENV{'FEELNCPATH'}."/bin/".$rprog;
pod2usage("Error: Cannot access FEELnc bin dir with path '$rprogpath'...\nCheck the environnment variable FEELNCPATH\n") unless( -r $rprogpath);
my $pathRscript = Utils::pathProg("Rscript");
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


# store mRNA annotation = hashref
my $refmrna;



warn "> Preparing files for random forest...\n";

##########################################################
# mRNA file
#######
# Training file
my $cdnafile = Utils::renamefile($mRNAfile, ".cdnatrain.fa");
my $orffile  = Utils::renamefile($mRNAfile, ".orftrain.fa");
my $lncfile;
if (defined $lncRNAfile){
	$lncfile = Utils::renamefile($lncRNAfile, ".lnctrain.fa");
} else {
	$lncfile = Utils::renamefile($mRNAfile, ".mRNAlinctrain.fa");
}


# add a refhash that will contain the mRNA ID that passed cDNA and ORF steps
# Will be used to checked for randomization
my $ref_cDNA_passed;

# if GTF
# ------
if ($mRNAfileformat eq "gtf"){

	# die if genome not specified
	pod2usage("Error: Cannot read your genome file '$genome' (-g option)...\nFor help, see:\n$progname --help\n") if (! -r $genome && !-d $genome);

	$refmrna  = Parser::parseGTF($mRNAfile, 'exon,CDS,stop_codon,start_codon', undef , \%biotype , $verbosity);
	my $sizeh = keys(%{$refmrna});


	die "Your input mRNA file '", basename($mRNAfile),"' contains only *$sizeh* transcripts.\nNot enough to train the program with the '--numtx|n $numtx' option (default option == 3000)\n" if ($sizeh < $numtx);
	print STDERR "\tYour input mRNA training file '", basename($mRNAfile),"' contains *$sizeh* transcripts\n" if ($verbosity > 0 );

	# Create cDNA and ORF 2 files for training and testing CPAT
	$ref_cDNA_passed = &CreateORFcDNAFromGTF($refmrna, $cdnafile, $orffile, $numtx, $genome, $verbosity); # for reproducibility

# if FASTA
# ------
} elsif ($mRNAfileformat eq "fasta") {

	# Create cDNA and ORF 2 files for training and testing CPAT
	&CreateORFcDNAFromFASTA($mRNAfile, $cdnafile, $orffile, $numtx, $verbosity);

} else {
	die "Error : Unrecognized format for annotation file '$mRNAfile'\n";
}

##########################################################
# lncRNA file
#######
# if file is defined, it means that we do not have to extract from intergenic
if (defined $lncRNAfile){

	my $computeORF       = undef; # we do not have to compute/extract ORF
	my $lncRNAfileformat = Utils::guess_format($lncRNAfile);

	# if GTF
	# ------
	if ($lncRNAfileformat eq "gtf"){

		my $reflnc = Parser::parseGTF($lncRNAfile, 'exon' , undef, undef, $verbosity);
		my $sizeh  = scalar keys(%{$reflnc});

		die "Your input lncRNA training file '", basename($lncRNAfile),"' contains only *$sizeh* transcripts.\nNot enough to train the program with the '--numtx|n $numtx' option (default option == 3000)\n" if ($sizeh < $numtx);
		print STDERR "\tYour lncRNA training file '", basename($lncRNAfile),"' contains *$sizeh* transcripts\n" if ($verbosity > 0 );

		# Create cDNA and ORF 2 files for training and testing CPAT
		$ref_cDNA_passed = &CreateORFcDNAFromGTF($reflnc, $lncfile, $computeORF, $numtx, $genome, $verbosity);

	# if FASTA
	# ------
	}elsif ($lncRNAfileformat eq "fasta") {

		# Create cDNA and ORF 2 files for training and testing CPAT
		&CreateORFcDNAFromFASTA($lncRNAfile, $lncfile, $computeORF, $numtx, $verbosity);

	} else {
		die "Error: Unrecognized format for lncRNA training file '$lncRNAfile'\n";
	}

} else { # lncRNA training file not defined

	print STDERR "> The lncRNA training file is not set...will extract intergenic region for training (can take a while...)\n";

	# Relocated mRNA sequence in intergenic regions to be used as a training lncRNA file

	&randomizedGTFtoFASTA ($refmrna, $ref_cDNA_passed, $lncfile, $genome, $numtx, $maxTries, $maxN, $verbosity);

}

#################################
# Launch RF on $infile in fasta
my $infile_outfa;
my $refin;
if (Utils::guess_format($infile) eq "gtf"){

	$refin        = Parser::parseGTF($infile, 'exon', undef , undef , $verbosity);
	$infile_outfa = $infile.".fa";
	ExtractFromHash::hash2fasta($refin, $genome, $infile_outfa,  $verbosity);

} elsif (Utils::guess_format($infile) eq "fasta"){
	$infile_outfa = $infile;

} else {
	die "Error: Unrecognized format for input file '$infile'...\n";
}


# VW modif crade !
# besoin des ORF pour lnc et test
my $lncOrfFile  = "./lncRNA_ORF.fa";
my $testOrfFile = "./test_ORF.fa";

# VW : Récupère les ORF du jeu lncRNA et test, crade !!!!
&CreateORFcDNAFromFASTA($lncfile, "/tmp/poubelle1", $lncOrfFile, $numtx, $verbosity);
# VW : en dur une valeur de 10000 pour avoir l'ensemble des ORF
#&CreateORFcDNAFromGTF($refin, "/tmp/poubelle2",  $testOrfFile, 10000, $genome, $verbosity);
&CreateORFcDNAFromGTF($refin, "/tmp/poubelle2",  $testOrfFile, undef, $genome, $verbosity);

print STDERR "> Run random Forest on '$infile_outfa':\n";
my $rfout = basename($infile)."_RF.out";


# VW: Run de façon crade !
RandomForest::runRF($cdnafile, $orffile, $lncfile, $lncOrfFile, $infile_outfa, $testOrfFile, $rfout, $kmerList, $rfcut, $verbosity, $keepTmp);
# Parsing of the random forest output
RandomForest::rfPredToGTF($infile,$rfout);

################################
##### END OF THE MAIN CODE #####
################################



###################################
##### DEFINITION OF FUNCTIONS #####
###################################

sub CreateORFcDNAFromGTF{

	my ($h, $cdnafile, $orffile, $nbtx, $genome, $verbosity) = @_;

	# Note if $orffile is not defined, we just extract cDNA
	# Note if $nbtx is undefined, we extract all ORF and cDNA

	my $orfob;
	my %h_orf;              # for storing and printing ORF sequence
	my %h_cdna;             # for storing and printing cDNA sequence
	my $allow_no_start = 0; # do not allow for CDS start not found
	my $allow_no_stop  = 0; # do not allow for CDS stop not found
	my $countseqok     = 0; # counter on good ORF (start and end found)
	my $filterforCDS   = 0; # get only line with CDS level


	for my $tr (sort keys(%{$h})){ # for reproducibility

		# shortcut for feature2seq sub
		my $chr    = $h->{$tr}->{'chr'};
		my $strand = $h->{$tr}->{'strand'};

		# Check Biotype
		my $biotype = $h->{$tr}->{'feature'}[0]->{'transcript_biotype'} if (defined $h->{$tr}->{'feature'}[0]->{'transcript_biotype'});


		# get cDNA sequence for transcript tr
		$filterforCDS = 0; # do we filter seq for CDS
		my $cdnaseq   = ExtractFromFeature::feature2seq($h->{$tr}->{'feature'}, $genome, $chr , $strand, $filterforCDS, $verbosity);
		die "ERROR: Tx '$tr' returns an empty sequence...\n" if (!defined $cdnaseq);
		#######################################
		# ORF
		if (defined $orffile){
			my $containCDS = ExtractFromFeature::checkCDS($h->{$tr}->{'feature'});
			if (! $containCDS ){
				warn "\tYour input GTF file does not contain CDS information... the program will extract the longest one for each transcript...\n" if ($countseqok < 1 && $verbosity > 5);
				# we create an ORF hash based on extraction of longest ORF
				$orfob = Orf::longestORF2($cdnaseq,$strand, $allow_no_start, $allow_no_stop, undef, 1);

			} else {
				warn "\tYour input GTF file does contain CDS information...\n" if ($countseqok < 1 && $verbosity > 5);
				$filterforCDS = 1; # we activate filter to get only CDS and stop codon DNA sequence
				my $orfseq    = ExtractFromFeature::feature2seq($h->{$tr}->{'feature'}, $genome, $chr , $strand, $filterforCDS, $verbosity);
				# we create an ORF hash
				$orfob        = Orf::orfSeq2orfOb($orfseq, $strand, $verbosity);

			}

			# Add ORF to a hash %h_orf only if the ORF is complete
			if ($orfob->{'check_start'} && $orfob->{'check_stop'}){
				$h_orf{$tr} = $orfob->{'cds_seq'};
				print STDERR "\tExtracting ORFs&cDNAs ", $countseqok++,"/$nbtx...\r" if( defined $nbtx);
				print STDERR "\tExtracting ORFs&cDNAs ", $countseqok++,"...\r"       if(!defined $nbtx);
			} else {
				warn "Tx: $tr ('$biotype') with CDS features: $containCDS is not complete...skipping for training\n" if ($verbosity > 10);
				next; # next if ORF is not OK
			}
		}

		#######################################
		# ADD cDNA only (if ORF is OK) : see next in above block
		# store cDNA seq
		if (!defined $orffile){
			print STDERR "\tExtracting cDNAs ", $countseqok++,"/$nbtx...\r" if( defined $nbtx);
			print STDERR "\tExtracting cDNAs ", $countseqok++,"...\r"       if(!defined $nbtx);
		}
		$h_cdna{$tr} = $cdnaseq;


		if (defined $nbtx && $countseqok == $nbtx){
			print STDERR "\tMax ORF/cDNAs sequences '$nbtx' reached..ending!\n";
			last;
		}


	}
	# if dedfined ORFfile, we write ORF and cDNA file
	if (defined $orffile){
		# Final Check if the number of complete ORF is ok
		my $sizehorf = keys(%h_orf);
		die "The number of complete ORF found with computeORF mode is *$sizehorf* transcripts... That's not enough to train the program\n" if (defined $nbtx && $sizehorf < $minnumtx);

		&writefastafile(\%h_orf,  $orffile, $verbosity);
		&writefastafile(\%h_cdna, $cdnafile, $verbosity);



	# we write only  cDNA file
	} else {

		my $sizeh = keys(%h_cdna);
		die "The number of cDNA sequences is *$sizeh* transcripts... That's not enough to train the program\n" if ($sizeh < $minnumtx);
		&writefastafile(\%h_cdna, $cdnafile, $verbosity);
	}

	return \%h_cdna;

}


sub CreateORFcDNAFromFASTA{

	my  ($fastafile, $cdnafile, $orffile, $nbtx, $verbosity)	=	@_;


	print STDERR "Extract ORF/cDNA from fasta file '$fastafile'..\n";

	my %h_orf;              # for storing and printing ORF sequence
	my %h_cdna;             # for storing and printing cDNA sequence

	my $allow_no_start = 0; # do not allow for CDS start not found
	my $allow_no_stop  = 0; # do not allow for CDS stop not found

	# counter for seq with ORF ok
	my $countseqok = 0;
	my $strand     = ".";

	# Create SeqIO objects
	my $seqin = Bio::SeqIO->new(-file => $fastafile, -format => "fasta");

	# count the nb of sequences
	my $nbseq = 0;
	$nbseq++ while( my $seq = $seqin->next_seq());
	die "Your input FASTA '$fastafile' contains only *$nbseq* sequences.\nNot enough to train the program (default option --ntx|-n)\n" if ($nbseq < $minnumtx);

	# weird have to recreate a seqio object
	$seqin = Bio::SeqIO->new(-file => $fastafile, -format => "fasta");

	# Go through each sequences
	while(my $seq = $seqin->next_seq()) {

		my $tr = $seq->id();

		# if not orf
		if (!defined $orffile){
			# store cDNA sequence
# 			my $new_seq = Bio::Seq->new(-id => $tr, -seq => $seq->seq(), -alphabet => 'dna');
			$h_cdna{$tr} = $seq->seq();
            print STDERR "\tExtracting cDNAs from FASTA ", $countseqok++,"/$nbtx complete cDNA(s)...\r";


		} else {# get also ORF

			my $orfob = Orf::longestORF2($seq->seq(),$strand, $allow_no_start, $allow_no_stop) if (defined $orffile);
			# Add ORF to a hash %h_orf only if the ORF is complete
			if ($orfob->{'check_start'} && $orfob->{'check_stop'}){
				$h_orf{$tr} = $orfob->{'cds_seq'};
				print STDERR "\tExtracting ORFs&cDNAs from FASTA ", $countseqok++,"/$nbtx complete ORF(s)...\r";

				# store cDNA sequence
# 				my $new_seq = Bio::Seq->new(-id => $tr, -seq => $seq, -alphabet => 'dna');
				$h_cdna{$tr} = $seq->seq();

			} else {
				warn "Tx: $tr : ORF is not complete...skipping for training\n" if ($verbosity > 5);
				next; # next if ORF is not OK
			}
		}
		# Check if numtx is reached
		if ($countseqok == $nbtx){
			print STDERR "Max cDNAs/ORF sequences '$nbtx' reached..ending!\n";
			last;
		}

	}
	# if dedfined ORFfile, we write ORF and cDNA file
	if (defined $orffile){

		# Final Check if the number of complete ORF is ok
		my $sizehorf = keys(%h_orf);
		die "The number of complete ORF found with computeORF mode is *$sizehorf* ... That's not enough to train the program\n" if ($sizehorf < $minnumtx);

		&writefastafile(\%h_orf,  $orffile, $verbosity);
		&writefastafile(\%h_cdna, $cdnafile, $verbosity);


	# we write only  cDNA file
	} else {

		my $sizeh = keys(%h_cdna);
		die "The number of cDNA sequences is *$sizeh* transcripts... That's not enough to train the program\n" if ($sizeh < $minnumtx);
		&writefastafile(\%h_cdna, $cdnafile, $verbosity);
	}
}


sub writefastafile{

	my ($h, $filename, $verbosity) = @_;

	print STDERR "\tWriting FASTA file '$filename'\n" if ($verbosity > 5);

	# cDNA
	my $seq = Bio::SeqIO ->new(-format => 'fasta', -file => '>'.$filename, -alphabet =>'dna');
	foreach my $id (keys %{$h}){
		my $new_seq = Bio::Seq->new(-id => $id, -seq => $h->{$id});
		$seq->write_seq($new_seq);
	}

}

sub randomizedGTFtoFASTA{

	my ($h, $ref_cDNA_passed, $cdnafile, $genome, $nbtx, $maxTries, $maxN, $verbosity) = @_;

	$nbtx      ||= 1000; # number of random tx required
	$maxTries  ||= 10;   # max tries to for computing both overlap and N
	$maxN      ||= 5;    # Proportion (in 100%) of N's authorized in new random sequence
	$verbosity ||= 0;


	my $split         = 1;
	my $hlightforover = Parser::GTF2GTFgnlight ($h, $split, $verbosity);

	# Get genome sequences size
	print STDERR "- Get chromosome sizes \n" if ($verbosity > 0);
	my $db = Bio::DB::Fasta->new($genome);
	my $refgenomesize;
	foreach my $id ( $db->ids){
		next if ($id =~ /^AAEX|^JH/ ); # for dog chromosome
		next if ($id =~ /^KI|^GL/ ); # for human chromosome GRCh38 : Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
		$refgenomesize->{$id} = $db->length($id); # populate hash with id => seq_length
	}
# 	print Dumper $refgenomesize;

	#  hashref tx annotation sizes
	my $refannotsize = ExtractFromHash::getCumulSizeFromGtfHash ($h,$verbosity, 0);

# 	print Dumper $refannotsize;

	print STDERR "- Relocate Transcripts \n" if ($verbosity > 0);
	my $i = 0;
	my $h_transcript_size = keys(%{$h});

	my %h_cdna_rdm; # to store correclty relocated sequences

	TX:
	foreach my $tx (sort keys %{$refannotsize}){ # sort for reproducibility

		next if ( ! exists $ref_cDNA_passed->{$tx}); # only keep mRNA tx that are in the cDNA fasat file for sorting CPAT :  for reproducibility


		my $overlap    = 1; # Initialize variable for iterative search for selfoverlap
		my $includeN   = 1; # Initialize variable for iterative search for N
		my $countTries = 0; # Number of tries

		# data for new sequence
		my ($chrrdm, $beg, $end, $seq);
		$seq = ""; # new fasta sequence ==> initialize in case pb with bio::db index

		if (defined $nbtx && $i == $nbtx){
			print STDERR "- Max number of transcripts (--nbtx == $nbtx) reached... ending!\n";
			last;
		}

		# while there is an overlap with known annotation
		while ($overlap || $includeN){

			# maxTries
			$countTries++;
			if ( $countTries ==  $maxTries){
				print  STDERR "MaxTries reached ($maxTries) for $tx...skipping it\n";
				next TX;
			}

			my $seed = $i+$countTries; # the seed is initiated accroding to the $i (tx) and the nb of try... if only $i, the same chr:pos will be returned...
			# Initialize srand foreach tx
			srand($seed);
			# define a rand indice for all chr hash
			my $randindex = int( rand(scalar keys %{$refgenomesize}) );
			my @chrrdm    = sort keys(%{$refgenomesize}); # sort hash for reproducibility
			$chrrdm       = $chrrdm[$randindex];

			# define a random start/begin position on the random chr (and thus the end)
			$beg = int(rand($refgenomesize->{$chrrdm}));
			$end = $beg + int( $refannotsize->{$tx}->{size} * $sizecorrec);

			#VW print "chr: $chrrdm; begin: $beg; end: $end\n";

			# Self - Overlap
			$overlap = overlapwithH($chrrdm,$beg,$end, $hlightforover, $countTries, $verbosity);
			if ($overlap){
				next;
			} else{
				$overlap =0;
			}

			# Test for Ns
			#############
			my $propN;
			($propN,$seq) = getPropN($chrrdm,$beg,$end, $db, 'N');
			if ($propN == -1){
				warn "Try: $countTries -> Extract sequences for $tx ($chrrdm:$beg-$end) returns an undefined sequence... skipping it\n" if ($verbosity > 10);
			} elsif ($propN > $maxN){
				warn "Try: $countTries -> Extract sequences for $tx ($chrrdm:$beg-$end) returns a $propN % with N!... skipping it\n" if ($verbosity > 10);
			}else {
				$includeN = 0;
			}
		}
		# Write New random sequence
		my $id           = $tx."_random_($chrrdm:$beg-$end)";
		$h_cdna_rdm{$id} = $seq;

		# verbosity
		$i++;
		if ($verbosity > 0){
			Utils::showProgress($nbtx, $i, "Print ".$tx.": ");
		}
	}

	my $sizeh = keys(%h_cdna_rdm);
	die "The number of RANDOMLY relocated cDNA sequences =  *$sizeh* transcripts... That's not enough to train the program\n" if ($sizeh < $minnumtx);
	&writefastafile(\%h_cdna_rdm, $cdnafile, $verbosity);

}

# test for overlap between a chr:start-end and a refh splited by chr
sub overlapwithH{

	my ($chr,$start,$end, $rehchr, $count, $verbosity)	= @_;

	my $overlap = 0;
	if (exists $rehchr->{$chr}){ # for the chromosome in the annotation test overlap

		my $refhchr = $rehchr->{$chr};

		# Test for overlap with annotation $h
		foreach my $locus (ExtractFromHash::sortGnsStartg($refhchr)){

			my $annbeg = $rehchr->{$chr}->{$locus}->{"startg"};
			my $annend = $rehchr->{$chr}->{$locus}->{"endg"};
			my $strand = $rehchr->{$chr}->{$locus}->{"strand"};

			# trick to speed  loop
			next if ($annend < $start);
			if      ($annbeg > $end){
				$overlap = 0;
				last;
			}

			# test overlap
			$overlap = Utils::foverlap($start,$end,$annbeg,$annend, $strand, ".", 0);

			if ($overlap){
				print STDERR "Try: $count -> Overlap $chr:$start-$end -- $chr:$annbeg-$annend ($strand) $locus \n" if ($verbosity > 10);
				last;
			}
		}
	} else { # if new chromosome is not in the exclusion file (chr without feature)
		$overlap = 0;
	}

	return $overlap;

}

# get proportion of N ($nucleotide) in a sequence defined by
# -$chr,$start,$end,
# -$db a bio::db::fasta object
sub getPropN{

	my ($chr,$start,$end, $db, $nucleotide) = @_;

	my $propN = -1; # default values
	my $seq   = "";

	# Get sequence
	$seq = $db->seq($chr, $start => $end);
	# test if good sequence
	if ($seq eq ""){
		warn "getPropN:: Sequence ($chr:$start-$end) returns an empty string!...skipping it\n";
	} else {
		my $numberofN = () = $seq  =~ /$nucleotide/gi;
		$propN        = int( $numberofN *100 / ($end-$start) );
	}

	return ($propN, $seq);

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
  -k,--kmer="2,3,4,5,6"			Kmer size list with size separate by ',' as string [ default "2,3,4,5,6" ]
  --keeptmp=0				To keep the temporary files, by default don't keep it (0 value). Any other value than 0 will keep the temporary files

=head2 Intergenic lncRNA extraction

	-to be added


=head2 Log output

  -o,--outlog=file.log		Specify the log file of output which [default infile.log]



=head1 AUTHORS

=over 4

=item *

Thomas DERRIEN <tderrien@univ-rennes1.fr>
Fabrice LEGEAI <fabrice.legeai@inria.fr>

=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2014 by IGDR - CNRS

=cut
