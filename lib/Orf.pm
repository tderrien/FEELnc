package Orf;

$VERSION = v0.0.1;

# This package is meant to create an ORF object
# Get longest one

use strict;
use warnings;
use Bio::DB::Fasta;
use Bio::Seq;
use Data::Dumper;

use Utils;


$| = 1;


# Select best ORF depending on strict mode (with Met and stop) and relaxed (without) and cutoff size
sub chooseORF{
    my ($seq, $cutoff) = @_;
    $cutoff //=0.75; #  strict ORF will be considered if sizeORF_strict/sizeORF_relaxed > this cutoff

    die "Orf::chooseORF : Cannot read sequence '$seq'\n..." unless ($seq =~/[actg]/gi);

    my $ostrict  = longestORF2($seq,".",0,0);
    my $orelaxed = longestORF2($seq,".",1,1);
	
    my $orf_selected = $orelaxed;
    
    # if exists a strict orf
	if ($ostrict){
    	#  force ORF to be strict if size > cutoff% of ORF relaxed
	   if ( ($ostrict->{'orflength'}/$orelaxed->{'orflength'})  > $cutoff ) {
			$orf_selected = $ostrict;
		}
	}
    return $orf_selected;
}


# adapted from bioperl example : script/longorf.pl
# (c) Dan Kortschak
sub longestORF {
    my ($seq, $strand, $allow_no_start, $allow_no_stop, $reverse_strand) = @_;
    $strand         //= "."; # to put in ORF object
    $allow_no_start //= 0;   # 1 to allow longest ORF not starting with Met
    $allow_no_stop  //= 0;   # 1 to allow longest ORF not finishing with stop_codon
    $reverse_strand //= 0;   # 1 if reverse complement input sequence

    die "Orf::longestORF: Your input sequence is empty \n" if ($seq eq "");
    die "Orf::longestORF: Your input sequence is not DNA \n" if ($seq =~ /!(a|c|t|g)/gi);

    #  revcomp sequence
    if ($reverse_strand){
	$seq = Utils::getRevComp($seq);
	if ($strand eq "+"){$strand = "-"};
	if ($strand eq "-"){$strand = "+"};
    }

    # variables
    my $seqlength = length($seq);
    my @starts=();
    my @ends=();
    my $best=0;
    my $bestorf	="";
    my $seqprot ="";

    # no_start
    if ($allow_no_start){
	for (my $frame=0;$frame<3;$frame++) {
	    unless ($seq=~m/^.{$frame}(taa|tga|tag)/i) {
		push @starts,$frame+1;
	    }
	}
    }

    # Get starts
    while ($seq =~m/(atg)/gi) {
	push @starts,pos($seq)-2;
    }
    # Get stops
    while ($seq=~m/(taa|tga|tag)/gi) {
	push @ends,pos($seq)-2;
    }

    # no_stop in last 3 nt of the cDNA sequence
    if ($allow_no_stop){
	push @ends,($seqlength-2, $seqlength-1, $seqlength);
    }


    # Inspired by Transdecoder
    my %last_delete_pos = ( 0=>-1,
			    1=>-1,
			    2=>-1); #store position of last chosen stop codon in spec reading frame.

    my $bests=0;
    my $beste=0;

    # Get longest ORF
    for my $s (@starts) {
	my $start_pos_frame = $s % 3;

	for my $e (@ends) {

	    if ($e%3==$s%3 &&
		$e>$s &&
		($s > $last_delete_pos{$start_pos_frame})){ #only count each stop once.

		print "e:$e - start: $s > $last_delete_pos{$start_pos_frame} $start_pos_frame\n";
		$last_delete_pos{$start_pos_frame} = $e;

		if ($e-$s>$best) {
		    $best=$e-$s;
		    ($bests,$beste)=($s-1,$e+2); # **** Here we use 0-based coordinates and include the stop codon ***
		    $bestorf	= substr($seq, $bests, $beste-$bests);
		}
		last
	    } else {
		next
	}
	}
    }
    # Get translation
    my $seq_obj = Bio::Seq->new(-seq => $bestorf, -alphabet => 'dna' );
    $seqprot	= $seq_obj->translate->seq;

    # Check start and stop
    my $checkstart = check_start_codon($seqprot);
    my $checkstop  = check_stop_codon($seqprot);

    # Store all data in hash
    my %orfobj =(start		=> $bests,
		 end			=> $beste,
		 seqlength 	=> $seqlength,
		 cds_seq 	=> $bestorf,
		 prot_seq 	=> $seqprot,
		 strand 		=> $strand,
		 check_start => $checkstart,
		 check_stop  => $checkstop
	);

    return \%orfobj;

}


sub longestORF2{
    my ($seq, $strand, $allow_no_start, $allow_no_stop, $reverse_strand, $maxstopskip) = @_;
    $strand         //= "."; # to put in ORF object
    $allow_no_start //= 0;   # 1 to allow longest ORF not starting with Met
    $allow_no_stop  //= 0;   # 1 to allow longest ORF not finishing with stop_codon
    $reverse_strand //= 0;   # 1 if reverse complement input sequence
    $maxstopskip    //= 1;   # number of stop codon to skip : 1 means we stop ORF at the first stop codon

    die "Orf::longestORF: Your input sequence is empty \n" if ($seq eq "");
    die "Orf::longestORF: Your input sequence is not DNA \n" if ($seq =~ /!(a|c|t|g)/gi);

    #  revcomp sequence
    if ($reverse_strand){
	$seq = Utils::getRevComp($seq);
	if ($strand eq "+"){$strand = "-"};
	if ($strand eq "-"){$strand = "+"};
    }

    # variables
    my $seqlength = length($seq);
    my @starts=();
    my @ends=();
    my $best=0;
    my $bestorf	="";
    my $seqprot ="";

    # no_start
    if ($allow_no_start){
	for (my $frame=0;$frame<3;$frame++) {
	    unless ($seq=~m/^.{$frame}(taa|tga|tag)/i) {
		push @starts,$frame+1;
	    }
	}
    }

    # Get starts
    while ($seq =~ /atg/gi) {
	push @starts, $-[0]+1;
    }
    # Get stops
    while ($seq =~ /taa|tga|tag/gi) {
	push @ends, $-[0]+1 ;
    }

    # no_stop in last 3 nt of the cDNA sequence
    if ($allow_no_stop){
	push @ends,($seqlength-2, $seqlength-1, $seqlength);
    }

    # Inspired by Transdecoder
    my %last_delete_pos = ( 0=>-1,
			    1=>-1,
			    2=>-1); #store position of last chosen stop codon in spec reading frame.

    my $bests=0;
    my $beste=0;

    # Get longest ORF
    for my $s (@starts) {
	my $start_pos_frame = $s % 3;

	my $nbfalsestop = 0;
	for my $e (@ends) {

	    if ($e%3==$s%3 && $e>$s) {

		if ($s > $last_delete_pos{$start_pos_frame}){ #only count each stop once.

		    if ($e-$s>$best) {
			$best=$e-$s;
			($bests, $beste) = ($s-1, $e+2); # **** Here we use 0-based coordinates and include the stop codon ***
			$bestorf	= substr($seq, $bests, $beste-$bests);
		    }
		    # We only store the last stop codon if we reach the nb of stop codon to pass by
		    $nbfalsestop++;
		    if ($nbfalsestop == $maxstopskip){
			$last_delete_pos{$start_pos_frame} = $e;
		    }
		}
	    } else {
		next;
	    }
	}
    }

    # Check for the lenght of the orf seq
    my $rest = length($bestorf)%3;

    if($rest != 0)
    {
	if(check_start_codon_dna($bestorf))
	{
	    $bestorf = substr($bestorf, 0, -$rest);
	}
	elsif(check_stop_codon_dna($bestorf))
	{
	    $bestorf = substr($bestorf, $rest);
	}
    }

    # Get translation
    my $seq_obj = Bio::Seq->new(-seq => $bestorf, -alphabet => 'dna' );
    $seqprot    = $seq_obj->translate->seq;

    # Check start and stop
    my $checkstart = check_start_codon($seqprot);
    my $checkstop  = check_stop_codon($seqprot);

    # Store all data in hash
    my %orfobj =(start       => $bests,
		 end         => $beste,
		 seqlength   => $seqlength,
		 cds_seq     => $bestorf,
		 orflength   => length($bestorf),
		 prot_seq    => $seqprot,
		 strand      => $strand,
		 check_start => $checkstart,
		 check_stop  => $checkstop
	);

    if(length($bestorf) > 0)
    {
	return \%orfobj;
    }
    else
    {
	return undef;
    }

}

sub orfSeq2orfOb{
    my ($seq, $strand, $verbosity) = @_;
    $strand    //= ".";
    $verbosity //= 1;
    
    # Test parsing i.e empty reafarray
    die "Orf::orfSeq2orfOb => $seq is weird...\n" if (not defined $seq || $seq eq "");

    # Get translation
    my $seq_obj = Bio::Seq->new(-seq => $seq, -alphabet => 'dna' );
    my $seqprot	= $seq_obj->translate->seq;

    # Check start and stop
    my $checkstart = check_start_codon($seqprot);
    my $checkstop  = check_stop_codon($seqprot);

    # Put all data in hash
    my %orfobj =(start       => undef,
		 end         => undef,
		 seqlength   => undef, # we do not known the length of the mRNA
		 cds_seq     => $seq,
		 prot_seq    => $seqprot,
		 strand      => $strand,
		 check_start => $checkstart,
		 check_stop  => $checkstop
	);

    return \%orfobj;
}

sub checkDnaStartStop {
    my ($seq, $verbosity) = @_;
    $verbosity //= 1;
    
    print STDERR "Orf::checkDnaStartStop => Sequence is empty" if ($seq eq "" && $verbosity > 1);

    if (defined ($seq) && $seq=~m/^ATG/i && $seq =~m/TAA$|TAG$|TGA$/i){
	return 1;
    } else {
	return 0;
    }
}

# Function that takes a protein sequence and
# return 1 if sequence start with Met, 0 otherwise
sub check_start_codon{
    my ($seqprot) = @_;

    if ($seqprot ne "" && $seqprot=~m/^M/){
	return 1;
    } else {
	return 0;
    }
}

# Function that takes a dna sequence and
# return 1 if sequence start with ATG, 0 otherwise
sub check_start_codon_dna{
    my ($seqdna) = @_;

    if ($seqdna ne "" && $seqdna=~m/^atg/gi){
	return 1;
    } else {
	return 0;
    }
}


# Function that takes a protein sequence and
# return 1 if sequence stop with *, 0 otherwise
sub check_stop_codon{
    my ($seqprot) = @_;

    if ($seqprot ne "" && $seqprot =~m/\*$/){
	return 1;
    } else {
	return 0;
    }
}

# Function that takes a dna sequence and
# return 1 if sequence stop with taa|tga|tag, 0 otherwise
sub check_stop_codon_dna{
    my ($seqdna) = @_;

    if ($seqdna ne "" && $seqdna =~m/taa|tga|tag$/gi){
	return 1;
    } else {
	return 0;
    }
}

1;
