package  ExtractCdnaOrf;

$version = v0.0.1;

# Perl libs
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Data::Dumper;

# lib directory: {FEELnc github directory}/lib/
use Parser;
use ExtractFromHash;
use ExtractFromFeature;
use Intersect;
use Utils;
use Orf;
use List::Util 'shuffle';


# Compare length of two or three orf objects and return the longest
sub compOrfLen
{
    my ($obj1, $ty1, $obj2, $ty2, $obj3, $ty3) = @_;
    $obj3 //= undef;
    $ty3  //= undef;

    die "Error compOrfLen: need at least two ORF object. Exit.\n" if( (! defined $obj1) || (! defined $ty1) || (! defined $obj2) || (! defined $ty2) );

    if( (! defined $obj3) )
    {
	if($obj1->{'orflength'} >= $obj2->{'orflength'})
	{
	    return( ($obj1, $ty1) );
	}
	else
	{
	    return( ($obj2, $ty2) );
	}
    }
    else
    {
	if($obj1->{'orflength'} >= $obj2->{'orflength'} && $obj1->{'orflength'} >= $obj3->{'orflength'})
	{
	    return( ($obj1, $ty1) );
	}
	elsif($obj2->{'orflength'} >= $obj1->{'orflength'} && $obj2->{'orflength'} >= $obj3->{'orflength'})
	{
	    return( ($obj2, $ty2) );
	}
	else
	{
	    return( ($obj3, $ty3) );
	}
    }

    return( (undef, undef) );
}

# Return $typeOrf if there is an ORF extracted, -1 if no ORF found (regarding the parameters)
sub getTypeOrf
{
    my ($name, $seq, $str, $refOrf, $type, $kmerMax) = @_;
    my $orfob0;
    my $orfob1;
    my $orfob2;
    my $orfob4;
    my $best;
    my $bestTy;
    my $flag0 = 0;
    my $flag1 = 0;
    my $flag2 = 0;
    my $flag4 = 0;

    return(-1) if($seq eq "");
    # -- if the sequence is empty, return -1

    # Type 0
    if($type==0)
    {
	$orfob0 = Orf::longestORF2($seq,$str, 0, 0, undef, 1);

	if(defined $orfob0 && $orfob0->{'orflength'} >= 2*$kmerMax)
	{
	    $flag0 = 1;
	}

	if($flag0==1)
	{

	    $refOrf->{$name} = $orfob0->{'cds_seq'};
	    return(0);
	}
	else
	{
	    return(-1);
	}
    }

    # Type 1
    if($type==1)
    {
	$orfob0 = Orf::longestORF2($seq,$str, 0, 0, undef, 1);
	$orfob1 = Orf::longestORF2($seq,$str, 0, 1, undef, 1);

	if(defined $orfob0 && $orfob0->{'orflength'} >= 2*$kmerMax)
	{
	    $flag0 = 1;
	}
	if(defined $orfob1 && $orfob1->{'orflength'} >= 2*$kmerMax)
	{
	    $flag1 = 1;
	}

	if($flag0==1 && $flag1==1)
	{
	    ($best, $bestTy) = &compOrfLen($orfob0, 0, $orfob1, 1);
	    $refOrf->{$name} = $best->{'cds_seq'};
	    return($bestTy);
	}
	elsif($flag0==1 && $flag1==0)
	{
	    $refOrf->{$name} = $orfob0->{'cds_seq'};
	    return(0);
	}
	elsif($flag0==0 && $flag1==1)
	{
	    $refOrf->{$name} = $orfob1->{'cds_seq'};
	    return(1);
	}
	else
	{
	    return(-1);
	}
    }

    # Type 2
    if($type==2)
    {
	$orfob0 = Orf::longestORF2($seq,$str, 0, 0, undef, 1);
	$orfob2 = Orf::longestORF2($seq,$str, 1, 0, undef, 1);
	if(defined $orfob0 && $orfob0->{'orflength'} >= 2*$kmerMax)
	{
	    $flag0 = 1;
	}
	if(defined $orfob2 && $orfob2->{'orflength'} >= 2*$kmerMax)
	{
	    $flag2 = 1;
	}

	if($flag0==1 && $flag2==1)
	{
	    ($best, $bestTy) = &compOrfLen($orfob0, 0, $orfob2, 2);
	    $refOrf->{$name} = $best->{'cds_seq'};
	    return($bestTy);
	}
	elsif($flag0==1 && $flag2==0)
	{
	    $refOrf->{$name} = $orfob0->{'cds_seq'};
	    return(0);
	}
	elsif($flag0==0 && $flag2==1)
	{
	    $refOrf->{$name} = $orfob2->{'cds_seq'};
	    return(2);
	}
	else
	{
	    return(-1);
	}
    }

    # Type 3 or 4
    if($type==3 || $type==4)
    {
	$orfob0 = Orf::longestORF2($seq,$str, 0, 0, undef, 1);
	$orfob1 = Orf::longestORF2($seq,$str, 0, 1, undef, 1);
	$orfob2 = Orf::longestORF2($seq,$str, 1, 0, undef, 1);
	if(defined $orfob0 && $orfob0->{'orflength'} >= 2*$kmerMax)
	{
	    $flag0 = 1;
	}
	if(defined $orfob1 && $orfob1->{'orflength'} >= 2*$kmerMax)
	{
	    $flag1 = 1;
	}
	if(defined $orfob2 && $orfob2->{'orflength'} >= 2*$kmerMax)
	{
	    $flag2 = 1;
	}

	if($flag0==1 && $flag1==1 && $flag2==1)
	{
	    ($best, $bestTy) = &compOrfLen($orfob0, 0, $orfob1, 1, $orfob2, 2);
	    $refOrf->{$name} = $best->{'cds_seq'};
	    return($bestTy);
	}
	elsif($flag0==1 && $flag1==1 && $flag2==0)
	{
	    ($best, $bestTy) = &compOrfLen($orfob0, 0, $orfob1, 1);
	    $refOrf->{$name} = $best->{'cds_seq'};
	    return($bestTy);
	}
	elsif($flag0==1 && $flag1==0 && $flag2==1)
	{
	    ($best, $bestTy) = &compOrfLen($orfob0, 0, $orfob2, 2);
	    $refOrf->{$name} = $best->{'cds_seq'};
	    return($bestTy);
	}
	elsif($flag0==0 && $flag1==1 && $flag2==1)
	{
	    ($best, $bestTy) = &compOrfLen($orfob1, 1, $orfob2, 2);
	    $refOrf->{$name} = $best->{'cds_seq'};
	    return($bestTy);
	}
	elsif($flag0==1 && $flag1==0 && $flag2==0)
	{
	    $refOrf->{$name} = $orfob0->{'cds_seq'};
	    return(0);
	}
	elsif($flag0==0 && $flag1==1 && $flag2==0)
	{
	    $refOrf->{$name} = $orfob1->{'cds_seq'};
	    return(1);
	}
	elsif($flag0==0 && $flag1==0 && $flag2==1)
	{
	    $refOrf->{$name} = $orfob2->{'cds_seq'};
	    return(2);
	}
	elsif($type==3 && $flag0==0 && $flag1==0 && $flag2==0)
	{
	    return(-1);
	}
    }

    # Type 4
    if($type==4)
    {
	$orfob4 = Orf::longestORF2($seq,$str, 1, 1, undef, 1);
	if(defined $orfob4 && $orfob4->{'orflength'} >= 2*$kmerMax)
	{
	    $flag4 = 1;
	}

	if($flag4==1)
	{
	    $refOrf->{$name} = $orfob4->{'cds_seq'};
	    return(4);
	}
	else
	{
	    return(-1);
	}
    }

    return(-1);
}


sub CreateORFcDNAFromGTF
{
    my ($gtfFile, $cdnaFile, $orfFile, $nbtx, $minnumtx, $genome, $lineType, $refBiotype, $orfType, $verbosity, $kmerMax) = @_;
    # Note if $nbtx is undefined, we extract all ORF and cDNA

    # add default value for minnumtx otherwise error when launchin without option numtx
    $minnumtx  //= 100;
    $verbosity //= 1;

    # die if genome not specified
    pod2usage("Error: Cannot read your genome file '$genome' (-g option)...\nFor help, run with --help option\n") if (! -r $genome && !-d $genome);

    # Parse the GTF file
    my $refmrna  = Parser::parseGTF($gtfFile, $lineType, undef , $refBiotype , $verbosity);
    my $sizeh    = keys(%{$refmrna});

    # Die if not enough transcript for training
    die "Not enough to train the program with the '--nbtx|n $nbtx' option (minimum == 100)...\n" if ( defined $nbtx && ($nbtx <  $minnumtx) );
    die "Your input GTF file '", basename($gtfFile),"' contains only *$sizeh* transcripts !\n Not enough to train the program (minimum == 100)...\n " if ( $sizeh < $minnumtx);

    print STDERR "\tYour input GTF file '", basename($gtfFile),"' contains *$sizeh* transcripts\n" if ($verbosity > 0 );

    my $orfob;
    my %h_orf;              # for storing and printing ORF sequence
    my %h_cdna;             # for storing and printing cDNA sequence
    my $countseqok     = 0; # counter on good ORF (start and end found)
    my $filterforCDS   = 0; # get only line with CDS level
    my $orfFlag        = 0; # get the result of getTypeOrf

    # for my $tr (sort keys(%{$refmrna}))  # for reproducibility
    for my $tr ( shuffle( sort( keys(%{$refmrna}) ) ) ) # sort before shuffle to get the same shuffle at each run
    {
	# shortcut for feature2seq sub
	my $chr    = $refmrna->{$tr}->{'chr'};
	my $strand = $refmrna->{$tr}->{'strand'};

	# get cDNA sequence for transcript tr
	$filterforCDS = 0; # do we filter seq for CDS
	my $cdnaseq   = ExtractFromFeature::feature2seq($refmrna->{$tr}->{'feature'}, $genome, $chr , $strand, $filterforCDS, $verbosity);
	if (!defined $cdnaseq)
	{
	    print STDERR "Warning: Transcript '$tr' returns an empty sequence... Skipping this transcripts\n" if ($verbosity > 1);
	    next;
	}

	### Get ORF
	my $containCDS = ExtractFromFeature::checkCDS($refmrna->{$tr}->{'feature'});
	if (! $containCDS )
	{
	    $orfFlag = &getTypeOrf($tr, $cdnaseq, $strand, \%h_orf, $orfType, $kmerMax);

	    # Print accordingly to getTypeOrf result
	    if ($orfFlag != -1)
	    {
		# Get the sequence only if an ORF is found
		$h_cdna{$tr} = $cdnaseq;
		$countseqok++;
		print STDERR "\tExtracting ORFs/cDNAs ", $countseqok,"/$nbtx...\r"                 if( defined $nbtx);
		print STDERR "\tExtracting ORFs/cDNAs ", $countseqok,"/".keys(%{$refmrna})."...\r" if(!defined $nbtx);
	    }
	    else
	    {
		print STDERR "No ORF found for transcript $tr... Skipping this transcript\n" if ($verbosity > 1);
		next; # next if ORF is not OK
	    }
	}
	else
	{
	    $filterforCDS = 1; # we activate filter to get only CDS and stop codon DNA sequence
	    my $orfseq    = ExtractFromFeature::feature2seq($refmrna->{$tr}->{'feature'}, $genome, $chr , $strand, $filterforCDS, $verbosity);
	    # Check if the length of the orf is greater than 2*kmerMax, if not then skip it
	    if(length($orfseq) < 2*$kmerMax)
	    {
		print STDERR "ORF found on transcript $tr is smaller than 2*$kmerMax... Skipping this transcript\n" if ($verbosity > 1);
		next; # next if ORF is not OK
	    }
	    # we create an ORF hash
	    $orfob        = Orf::orfSeq2orfOb($orfseq, $strand, $verbosity);
	    $h_orf{$tr}   = $orfob->{'cds_seq'};
	    $h_cdna{$tr}  = $cdnaseq;
	    $countseqok++;
	    print STDERR "\tExtracting ORFs/cDNAs ", $countseqok,"/$nbtx...\r"                 if( defined $nbtx );
	    print STDERR "\tExtracting ORFs/cDNAs ", $countseqok,"/".keys(%{$refmrna})."...\r" if(!defined $nbtx );
	}

	if (defined $nbtx && $countseqok == $nbtx) # Check for transcript extraction limit
	{
	    print STDERR "\n\tMax ORF/cDNAs sequences '$nbtx' reached.\n";
	    last;
	}
    }

    my $sizehorf = keys(%h_orf);
    die "The number of complete ORF found with computeORF mode is *$sizehorf* transcripts... That's not enough to train the program\n" if ($sizehorf < $minnumtx);

    if(!defined $nbtx)
    {
	print STDERR "\n\tExtracted '$countseqok' ORF/cDNAs sequences on '".keys(%{$refmrna})."'.\n";
    }

    # Write output FASTA files
    &write2fastafile(\%h_cdna, $cdnaFile, \%h_orf,  $orfFile, $verbosity);

    return (\%h_cdna, $refmrna);
}


sub CreateORFcDNAFromFASTA
{
    my ($fastaFile, $cdnaFile, $orfFile, $nbtx, $minnumtx, $orfType, $verbosity, $kmerMax) = @_;
    # If $nbtx is undef, extract all sequences
    $verbosity //= 1;

    print STDERR "Extract ORF/cDNA from fasta file '$fastaFile'..\n";

    my %h_orf;       # for storing and printing ORF sequence
    my %h_cdna;      # for storing and printing cDNA sequence
    my $orfFlag = 0; # get getTypeOrf result

    # counter for seq with ORF ok
    my $countseqok = 0;
    my $strand     = ".";

    # Create SeqIO objects
    my $seqin = Bio::SeqIO->new(-file => $fastaFile, -format => "fasta");

    # count the nb of sequences
    my $nbseq = 0;
    $nbseq++ while( my $seq = $seqin->next_seq());
    die "Your input FASTA '$fastaFile' contains only *$nbseq* sequences.\nNot enough to train the program (default option --ntx|-n)\n" if (defined $minnumtx && $nbseq < $minnumtx);

    # weird have to recreate a seqio object
    $seqin = Bio::SeqIO->new(-file => $fastaFile, -format => "fasta");

    # Go through each sequences
    while(my $seq = $seqin->next_seq()) {

	my $tr = $seq->id();

	$orfFlag = &getTypeOrf($tr, $seq->seq(), $strand, \%h_orf, $orfType, $kmerMax);

	# Print according to getTypeOrf result
	if ($orfFlag != -1)
	{
	    # Add cDNA only if an ORF is found
	    $h_cdna{$tr} = $seq->seq();
	    $countseqok++;
	    print STDERR "\tExtracting ORFs/cDNAs ", $countseqok,"/$nbtx...\r"  if( defined $nbtx);
	    print STDERR "\tExtracting ORFs/cDNAs ", $countseqok,"/$nbseq...\r" if(!defined $nbtx);
	}
	else
	{
	    print STDERR "No ORF found for transcript $tr... Skipping this transcript\n" if ($verbosity > 1);
	    next; # next if ORF is not OK
	}

	# Check if nbtx is reached
	if (defined $nbtx && $countseqok == $nbtx)
	{
	    print STDERR "\n\tMax cDNAs/ORF sequences '$nbtx' reached.\n";
	    last;
	}
    }

    # Final Check if the number of complete ORF is ok
    my $sizehorf = keys(%h_orf);
    die "The number of complete ORF found with computeORF mode is *$sizehorf*... That's not enough to train the program\n" if (defined $nbtx && $sizehorf < $minnumtx);

    if(!defined $nbtx)
    {
	print STDERR "\n\tExtracted '$countseqok' ORF/cDNAs sequences on '$nbseq'.\n";
    }

    # Write output FASTA files
    &write2fastafile(\%h_cdna, $cdnaFile, \%h_orf,  $orfFile, $verbosity);
}


sub writefastafile{
    my ($h, $filename, $verbosity) = @_;
    $verbosity //= 1;

    print STDERR "\tWriting FASTA file '$filename'\n" if ($verbosity > 1);

    # cDNA
    my $seq = Bio::SeqIO ->new(-format => 'fasta', -file => '>'.$filename, -alphabet =>'dna');
    # foreach my $id (sort keys %{$h}){
    foreach my $id ( sort( keys(%{$h}) ) ){
	my $new_seq = Bio::Seq->new(-id => $id, -seq => $h->{$id});
	$seq->write_seq($new_seq);
    }

}

sub write2fastafile{
    my ($h1, $filename1, $h2, $filename2, $verbosity) = @_;
    $verbosity //= 1;

    print STDERR "\tWriting FASTA file '$filename1'\n" if ($verbosity > 1);
    print STDERR "\tWriting FASTA file '$filename2'\n" if ($verbosity > 1);

    # cDNA
    my $seq1 = Bio::SeqIO ->new(-format => 'fasta', -file => '>'.$filename1, -alphabet =>'dna');
    # ORF
    my $seq2 = Bio::SeqIO ->new(-format => 'fasta', -file => '>'.$filename2, -alphabet =>'dna');
    # foreach my $id (sort keys %{$h}){
    foreach my $id ( shuffle( sort( keys(%{$h1}) ) ) )
    {
	my $new_seq1 = Bio::Seq->new(-id => $id, -seq => $h1->{$id});
	$seq1->write_seq($new_seq1);

	my $new_seq2 = Bio::Seq->new(-id => $id, -seq => $h2->{$id});
	$seq2->write_seq($new_seq2);
    }

}

sub randomizedGTFtoFASTA{
    my ($h, $ref_cDNA_passed, $cdnafile, $orffile, $genome, $nbtx, $minnumtx, $sizecorrec, $orfType, $maxTries, $maxN, $verbosity, $kmerMax) = @_;
    $maxTries  //= 10;   # max tries to for computing both overlap and N
    $maxN      //= 5;    # Proportion (in 100%) of N's authorized in new random sequence
    $verbosity //= 1;

    my $split         = 1;
    my $hlightforover = Parser::GTF2GTFgnlight ($h, $split, $verbosity);

    # Get genome sequences size
    print STDERR "\t\t- Get chromosome sizes \n" if ($verbosity > 1);
    my $db = Bio::DB::Fasta->new($genome);
    my $refgenomesize;
    foreach my $id ( $db->ids){
	$refgenomesize->{$id} = $db->length($id); # populate hash with id => seq_length
    }

    #  hashref tx annotation sizes
    my $refannotsize = ExtractFromHash::getCumulSizeFromGtfHash($h,$verbosity, 0);

    print STDERR "\t\t-Relocate Transcripts \n" if ($verbosity > 1);
    my $i                 = 0;
    my $h_transcript_size = keys(%{$h});

    my %h_cdna_rdm;  # to store correclty relocated sequences
    my %h_orf;       # to store the orf
    my %h_orf_tmp;   # to temporary store the orf
    my $orfFlag = 0; # to get the type of ORF
    my $cptok   = 0; # to get the number of extracted sequence


  TX:
    foreach my $tx ( shuffle( sort( keys(%{$refannotsize}) ) ) ) # sort before shuffle to get the same shuffle at each run
    {
	next if ( ! exists $ref_cDNA_passed->{$tx}); # only keep mRNA tx that are in the cDNA fasta file:  for reproducibility

	my $overlap    =  1; # Initialize variable for iterative search for selfoverlap
	my $includeN   =  1; # Initialize variable for iterative search for N
	my $countTries = -1; # Number of tries

	# data for new sequence
	my ($chrrdm, $beg, $end, $seq);
	$seq       = ""; # new fasta sequence ==> initialize in case pb with bio::db index
	my $seqORF = ""; # to temporary stock ORF sequence
	my $id;          # transcript name

	if (defined $nbtx && $i == $nbtx){
	    print STDERR "\n\tMax number of transcripts '$nbtx' reached.\n";
	    last;
	}

	# While there is an overlap with known annotation, N in seq or no ORF
	while ($overlap || $includeN || $orfFlag==-1){

	    # maxTries
	    $countTries++;
	    if ( $countTries ==  $maxTries){
		print STDERR "\n\tMaxTries reached ($maxTries) for $tx. Skipping it\n" if($verbosity > 1);
		next TX;
	    }

	    # define a rand indice for all chr hash
	    my $randindex = int( rand(scalar keys %{$refgenomesize}) );
	    my @chrrdm    = sort keys(%{$refgenomesize}); # sort hash for reproducibility
	    $chrrdm       = $chrrdm[$randindex];

	    # define a random start/begin position on the random chr (and thus the end)
	    $beg = int(rand($refgenomesize->{$chrrdm}));
	    $end = $beg + int( $refannotsize->{$tx}->{size} * $sizecorrec);

	    # Check if the size is almost equal to the original transcript
	    if ($end >  $refgenomesize->{$chrrdm})
	    {
		print STDERR "Try: $countTries -> Extract sequences for $tx ($chrrdm:$beg-$end) got an end greater than chromosome end... skipping it\n" if ($verbosity > 1);
		next;
	    }

	    # If the final length is < 200bp (smaller than lncRNA definition), then the size is set to 200
	    if($end - $beg < 200)
	    {
		$end = $beg + 200;
	    }

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
	    ($propN,$seq) = getPropN($chrrdm,$beg,$end, $db, 'N', $verbosity);
	    if ($propN == -1){
		print STDERR "Try: $countTries -> Extract sequences for $tx ($chrrdm:$beg-$end) returns an undefined sequence... Skipping it\n" if ($verbosity > 1);
	    } elsif ($propN > $maxN){
		print STDERR "Try: $countTries -> Extract sequences for $tx ($chrrdm:$beg-$end) returns a $propN % with N!... Skipping it\n" if ($verbosity > 1);
	    }else {
		$includeN = 0;
	    }

	    $id = $tx."_random_($chrrdm:$beg-$end)";
	    # Get ORF sequence
	    $orfFlag = &getTypeOrf($id, $seq, "+", \%h_orf_tmp, $orfType, $kmerMax);
	    $seqORF = $h_orf_tmp{$id};
	}
	# Write New random sequence
	$h_cdna_rdm{$id} = $seq;
	$h_orf{$id}      = $seqORF;
	$cptok++;
	print STDERR "\tExtracting ORFs/cDNAs with random coordinate ", $cptok,"/$nbtx...\r"                      if( defined  $nbtx);
	print STDERR "\tExtracting ORFs/cDNAs with random coordinate ", $cptok,"/".keys(%{$refannotsize})."...\r" if( !defined $nbtx);

	# verbosity
	$i++;
	if ($verbosity > 0)
	{
	    # Utils::showProgress($nbtx,                          $i, "Print ".$tx.": ") if(defined $nbtx);
	    # Utils::showProgress(scalar(keys(%{$refannotsize})), $i, "Print ".$tx.": ") if(!defined $nbtx);
	    Utils::showProgress($nbtx,                          $i, "Extracting intergenic sequences: ") if(defined $nbtx);
	    Utils::showProgress(scalar(keys(%{$refannotsize})), $i, "Extracting intergenic sequences: ") if(!defined $nbtx);
	}
    }

    my $sizeh = keys(%h_cdna_rdm);
    die "The number of RANDOMLY relocated cDNA sequences =  *$sizeh* transcripts... That's not enough to train the program\n" if ($sizeh < $minnumtx);

    if(defined $nbtx)
    {
	print STDERR "\n\tExtracted '$cptok' intergenic sequences on '$nbtx'.\n";
    }
    else
    {
	print STDERR "\n\tExtracted '$cptok' intergenic sequences on '".keys(%{$refannotsize})."'.\n";
    }

    # Write output FASTA files
    &write2fastafile(\%h_orf, $orffile, \%h_cdna_rdm, $cdnafile, $verbosity);
}


# Test for overlap between a chr:start-end and a refh splited by chr
sub overlapwithH{
    my ($chr,$start,$end, $rehchr, $count, $verbosity)	= @_;
    $verbosity //= 1;

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
		print STDERR "Try: $count -> Overlap $chr:$start-$end -- $chr:$annbeg-$annend ($strand) $locus \n" if ($verbosity > 1);
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
    my ($chr,$start,$end, $db, $nucleotide, $verbosity) = @_;
    $verbosity //= 1;

    my $propN = -1; # default values
    my $seq   = "";

    # Get sequence
    $seq = $db->seq($chr, $start => $end);
    # test if good sequence
    if ($seq eq ""){
	print STDERR "getPropN:: Sequence ($chr:$start-$end) returns an empty string... Skipping it\n" if($verbosity > 1);
    } else {
	my $numberofN = () = $seq  =~ /$nucleotide/gi;
	$propN        = int( $numberofN *100 / ($end-$start) );
    }

    return ($propN, $seq);
}

1;
