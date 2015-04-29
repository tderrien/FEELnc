package  ExtractCdnaOrf;

$VERSION = v0.0.1;

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


# Compare length of two or three orf objects and return the longest
sub compOrfLen
{
    my ($obj1, $ty1, $obj2, $ty2, $obj3, $ty3) = @_;
    $obj3 ||= undef;
    $ty3  ||= undef;

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
sub getTypeOrf2
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


# # Return $typeOrf if there is an ORF extracted, -1 if no ORF found (regarding the parameters)
# sub getTypeOrf2
# {
#     my ($name, $seq, $str, $refOrf, $type, $kmerMax) = @_;
#     my $orfob0;
#     my $orfob1;
#     my $orfob2;
#     my $orfob4;
#     my $flag0 = 0;
#     my $flag1 = 0;
#     my $flag2 = 0;

#     return(-1) if($seq eq "");
#     # -- if the sequence is empty, return -1

#     $orfob0 = Orf::longestORF2($seq,$str, 0, 0, undef, 1);
#     $orfob1 = Orf::longestORF2($seq,$str, 0, 1, undef, 1);
#     $orfob2 = Orf::longestORF2($seq,$str, 1, 0, undef, 1);
#     $orfob4 = Orf::longestORF2($seq,$str, 1, 1, undef, 1);


#     if($orfob0->{'check_start'} && $orfob0->{'check_stop'} && (length($orfob0->{'cds_seq'}) >= (2*$kmerMax)))
#     {
# 	$flag0 = 1;
#     }
#     if($orfob1->{'check_start'} && (length($orfob1->{'cds_seq'}) >= (2*$kmerMax)))
#     {
# 	$flag1 = 1;
#     }
#     if($orfob2->{'check_stop'} && (length($orfob2->{'cds_seq'}) >= (2*$kmerMax)))
#     {
# 	$flag2 = 1;
#     }


#     # Type 0
#     if($type==0 && $flag0==1)
#     {
# 	$refOrf->{$name} = $orfob0->{'cds_seq'};
# 	return(0);
#     }
#     elsif($type==0 && $flag0==0)
#     {
# 	return(-1);
#     }

#     # Type 1
#     if($type==1)
#     {
# 	if($flag0==1 && $flag1==1)
# 	{
# 	    if(length($orfob0->{'cds_seq'}) >= length($orfob1->{'cds_seq'}))
# 	    {
# 		$refOrf->{$name} = $orfob0->{'cds_seq'};
# 		return(0);
# 	    }
# 	    else
# 	    {
# 		$refOrf->{$name} = $orfob1->{'cds_seq'};
# 		return(1);
# 	    }
# 	}
# 	elsif($flag0==1 && $flag1==0)
# 	{
# 	    $refOrf->{$name} = $orfob0->{'cds_seq'};
# 	    return(0);
# 	}
# 	elsif($flag0==0 && $flag1==1)
# 	{
# 	    $refOrf->{$name} = $orfob1->{'cds_seq'};
# 	    return(1);
# 	}
# 	else
# 	{
# 	    return(-1);
# 	}
#     }

#     # Type 2
#     if($type==2)
#     {
# 	if($flag0==1 && $flag2==1)
# 	{
# 	    if(length($orfob0->{'cds_seq'}) >= length($orfob2->{'cds_seq'}))
# 	    {
# 		$refOrf->{$name} = $orfob0->{'cds_seq'};
# 		return(0);
# 	    }
# 	    else
# 	    {
# 		$refOrf->{$name} = $orfob2->{'cds_seq'};
# 		return(2);
# 	    }
# 	}
# 	elsif($flag0==1 && $flag2==0)
# 	{
# 	    $refOrf->{$name} = $orfob0->{'cds_seq'};
# 	    return(0);
# 	}
# 	elsif($flag0==0 && $flag2==1)
# 	{
# 	    $refOrf->{$name} = $orfob2->{'cds_seq'};
# 	    return(2);
# 	}
# 	else
# 	{
# 	    return(-1);
# 	}
#     }


#     # Type 3
#     if($type==3)
#     {
# 	if($flag0==1 && $flag1==1 && $flag2==1)
# 	{
# 	    if(length($orfob0->{'cds_seq'}) >= length($orfob1->{'cds_seq'}) && length($orfob0->{'cds_seq'}) >= length($orfob2->{'cds_seq'}))
# 	    {
# 		$refOrf->{$name} = $orfob0->{'cds_seq'};
# 		return(0);
# 	    }
# 	    elsif(length($orfob0->{'cds_seq'}) <= length($orfob1->{'cds_seq'}) && length($orfob1->{'cds_seq'}) >= length($orfob2->{'cds_seq'}))
# 	    {
# 		$refOrf->{$name} = $orfob1->{'cds_seq'};
# 		return(1);
# 	    }
# 	    else
# 	    {
# 		$refOrf->{$name} = $orfob2->{'cds_seq'};
# 		return(2);
# 	    }

# 	}
# 	elsif($flag0==1 && $flag1==1 && $flag2==0)
# 	{
# 	    if(length($orfob0->{'cds_seq'}) >= length($orfob2->{'cds_seq'}))
# 	    {
# 		$refOrf->{$name} = $orfob0->{'cds_seq'};
# 		return(0);
# 	    }
# 	    else
# 	    {
# 		$refOrf->{$name} = $orfob1->{'cds_seq'};
# 		return(1);
# 	    }
# 	}
# 	elsif($flag0==1 && $flag1==0 && $flag2==1)
# 	{
# 	    if(length($orfob0->{'cds_seq'}) >= length($orfob2->{'cds_seq'}))
# 	    {
# 		$refOrf->{$name} = $orfob0->{'cds_seq'};
# 		return(0);
# 	    }
# 	    else
# 	    {
# 		$refOrf->{$name} = $orfob2->{'cds_seq'};
# 		return(2);
# 	    }
# 	}
# 	elsif($flag0==0 && $flag1==1 && $flag2==1)
# 	{
# 	    if(length($orfob1->{'cds_seq'}) >= length($orfob2->{'cds_seq'}))
# 	    {
# 		$refOrf->{$name} = $orfob1->{'cds_seq'};
# 		return(1);
# 	    }
# 	    else
# 	    {
# 		$refOrf->{$name} = $orfob2->{'cds_seq'};
# 		return(2);
# 	    }
# 	}
# 	elsif($flag0==1 && $flag1==0 && $flag2==0)
# 	{
# 	    $refOrf->{$name} = $orfob0->{'cds_seq'};
# 	    return(0);
# 	}
# 	elsif($flag0==0 && $flag1==1 && $flag2==0)
# 	{
# 	    $refOrf->{$name} = $orfob1->{'cds_seq'};
# 	    return(1);
# 	}
# 	elsif($flag0==0 && $flag1==0 && $flag2==1)
# 	{
# 		$refOrf->{$name} = $orfob2->{'cds_seq'};
# 		return(2);
# 	}
# 	else
# 	{
# 	    return(-1);
# 	}
#     }

#     return(-1);
#     # if no ORF found, return -1
# }


# Return $typeOrf if there is an ORF extracted, -1 if no ORF found (regarding the parameters)
sub getTypeOrf
{
    my ($name, $seq, $str, $refOrf, $type, $kmerMax) = @_;
    my $orfob;
    my $orfob2;

    return(-1) if($seq eq "");
    # -- if the sequence is empty, return -1

    # Type 0
    $orfob = Orf::longestORF2($seq,$str, 0, 0, undef, 1);

    if($orfob->{'check_start'} && $orfob->{'check_stop'} && length($orfob->{'cds_seq'}) >= (2*$kmerMax))
	# -- if an ORF is found with a start and a stop codon
    {
	$refOrf->{$name} = $orfob->{'cds_seq'};
	return(0);
    }
    # Type 1
    if($type==1 || $type==3 || $type==4)
	# -- if type 1, 3 or 4, check for an ORF with a start codon
    {
	$orfob = Orf::longestORF2($seq,$str, 0, 1, undef, 1);
	if($type==1 && $orfob->{'check_start'} && length($orfob->{'cds_seq'}) >= (2*$kmerMax))
	    # -- if type 1 and a start codon is found, get this ORF
	{
	    $refOrf->{$name} = $orfob->{'cds_seq'};
	    return(1);
	}
    }
    # Type 2
    if($type==2 || $type==3 || $type==4)
	# -- if type 2, 3 or 4, check for an ORF with a stop codon
    {
	$orfob2 = Orf::longestORF2($seq,$str, 1, 0, undef, 1);
	if($type==2 && $orfob2->{'check_stop'} && length($orfob2->{'cds_seq'}) >= (2*$kmerMax))
	    # -- if type 2 and a stop codon is found, get this ORF
	{
	    $refOrf->{$name} = $orfob2->{'cds_seq'};
	    return(2);
	}
    }
    # Type 3
    if($type==3 || $type==4)
	# -- if type 3 or 4, take the longest ORF between type 1 and 2 (orfob and orfob2)
    {
	if(length($orfob->{'cds_seq'}) > length($orfob2->{'cds_seq'}) && length($orfob->{'cds_seq'}) >= (2*$kmerMax))
	    # if ORF with start codon >= ORF with stop codon, take ORF start, else take ORF stop
	{
	    $refOrf->{$name} = $orfob->{'cds_seq'};
	    return(3);
	}
	elsif(length($orfob2->{'cds_seq'}) >= length($orfob->{'cds_seq'}) && length($orfob2->{'cds_seq'}) >= (2*$kmerMax))
	{
	    $refOrf->{$name} = $orfob2->{'cds_seq'};
	    return(3);
	}
    }
    # Type 4
    if($type==4 && length($orfob->{'cds_seq'}) >= (2*$kmerMax))
	# -- if type 4, take the longest ORF whatever there is a start/stop codon
    {
	$orfob = Orf::longestORF2($seq,$str, 1, 1, undef, 1);
	$refOrf->{$name} = $orfob->{'cds_seq'};
	return(4);
    }

    return(-1);
    # if no ORF found, return -1
}


sub CreateORFcDNAFromGTF
{
    my ($gtfFile, $cdnaFile, $orfFile, $nbtx, $minnumtx, $genome, $lineType, $refBiotype, $orfType, $verbosity, $kmerMax) = @_;
    # Note if $nbtx is undefined, we extract all ORF and cDNA

    # die if genome not specified
    pod2usage("Error: Cannot read your genome file '$genome' (-g option)...\nFor help, run with --help option\n") if (! -r $genome && !-d $genome);

    # Parse the GTF file
    my $refmrna  = Parser::parseGTF($gtfFile, $lineType, undef , $refBiotype , $verbosity);
    my $sizeh = keys(%{$refmrna});

    # Die if not enough transcript for training
    die "Your input mRNA file '", basename($gtfFile),"' contains only *$sizeh* transcripts.\nNot enough to train the program with the '--nbtx|n $nbtx' option (default option == 3000)\n" if (defined $nbtx && $sizeh < $minnumtx);
    print STDERR "\tYour input training file '", basename($gtfFile),"' contains *$sizeh* transcripts\n" if ($verbosity > 0 );

    my $orfob;
    my %h_orf;              # for storing and printing ORF sequence
    my %h_cdna;             # for storing and printing cDNA sequence
    my $countseqok     = 0; # counter on good ORF (start and end found)
    my $filterforCDS   = 0; # get only line with CDS level
    my $orfFlag        = 0; # get the result of getTypeOrf

    for my $tr (sort keys(%{$refmrna}))  # for reproducibility
    {
	# shortcut for feature2seq sub
	my $chr    = $refmrna->{$tr}->{'chr'};
	my $strand = $refmrna->{$tr}->{'strand'};

	# get cDNA sequence for transcript tr
	$filterforCDS = 0; # do we filter seq for CDS
	my $cdnaseq   = ExtractFromFeature::feature2seq($refmrna->{$tr}->{'feature'}, $genome, $chr , $strand, $filterforCDS, $verbosity);
	die "ERROR: Tx '$tr' returns an empty sequence...\n" if (!defined $cdnaseq);

	### Get ORF
	my $containCDS = ExtractFromFeature::checkCDS($refmrna->{$tr}->{'feature'});
	if (! $containCDS )
	{
	    warn "\tYour input GTF file does not contain CDS information... the program will extract the longest one for each transcript...\n" if ($countseqok < 1 && $verbosity > 5);
	    $orfFlag = &getTypeOrf2($tr, $cdnaseq, $strand, \%h_orf, $orfType, $kmerMax);

	    # Print accordingly to getTypeOrf result
	    if ($orfFlag != -1)
	    {
		# Get the sequence only if an ORF is found
		$h_cdna{$tr} = $cdnaseq;
		print STDERR "\tExtracting ORFs&cDNAs ", $countseqok++,"/$nbtx...\r"                 if( defined $nbtx);
		print STDERR "\tExtracting ORFs&cDNAs ", $countseqok++,"/".keys(%{$refmrna})."...\r" if(!defined $nbtx);
	    }
	    else
	    {
		warn "Tx: $tr without CDS features: $containCDS is not complete...skipping for training\n" if ($verbosity > 10);
		next; # next if ORF is not OK
	    }
	}
	else
	{
	    warn "\tYour input GTF file does contain CDS information...\n" if ($countseqok < 1 && $verbosity > 5);
	    $filterforCDS = 1; # we activate filter to get only CDS and stop codon DNA sequence
	    my $orfseq    = ExtractFromFeature::feature2seq($refmrna->{$tr}->{'feature'}, $genome, $chr , $strand, $filterforCDS, $verbosity);
	    # we create an ORF hash
	    $orfob        = Orf::orfSeq2orfOb($orfseq, $strand, $verbosity);
	    $h_orf{$tr}   = $orfob->{'cds_seq'};
	    $h_cdna{$tr}  = $cdnaseq;
	    print STDERR "\tExtracting ORFs&cDNAs ", $countseqok++,"/$nbtx...\r"                 if( defined $nbtx);
	    print STDERR "\tExtracting ORFs&cDNAs ", $countseqok++,"/".keys(%{$refmrna})."...\r" if(!defined $nbtx);
	}

	if (defined $nbtx && $countseqok == $nbtx) # Check for transcript extraction limit
	{
	    print STDERR "\tMax ORF/cDNAs sequences '$nbtx' reached..ending!\n";
	    last;
	}
    }

    my $sizehorf = keys(%h_orf);
    die "The number of complete ORF found with computeORF mode is *$sizehorf* transcripts... That's not enough to train the program\n" if (defined $nbtx && $sizehorf < $minnumtx);

    # Write fasta file
    &writefastafile(\%h_orf,  $orfFile, $verbosity);
    &writefastafile(\%h_cdna, $cdnaFile, $verbosity);

    return \%h_cdna;
}

# sub OldCreateORFcDNAFromGTF{

#     my ($h, $cdnafile, $orffile, $nbtx, $genome, $orfType, $verbosity) = @_;

#     # Note if $orffile is not defined, we just extract cDNA


# 	#######################################
# 	# ADD cDNA only (if ORF is OK) : see next in above block
# 	# store cDNA seq
# 	if (!defined $orffile){
# 	    print STDERR "\tExtracting cDNAs ", $countseqok++,"/$nbtx...\r" if( defined $nbtx);
# 	    print STDERR "\tExtracting cDNAs ", $countseqok++,"...\r"       if(!defined $nbtx);
# 	}
# 	$h_cdna{$tr} = $cdnaseq;


# 	if (defined $nbtx && $countseqok == $nbtx){
# 	    print STDERR "\tMax ORF/cDNAs sequences '$nbtx' reached..ending!\n";
# 	    last;
# 	}


#     }
#     # if dedfined ORFfile, we write ORF and cDNA file
#     if (defined $orffile){
# 	# Final Check if the number of complete ORF is ok


# 	# we write only  cDNA file
#     } else {

# 	my $sizeh = keys(%h_cdna);
# 	die "The number of cDNA sequences is *$sizeh* transcripts... That's not enough to train the program\n" if ($sizeh < $minnumtx);
# 	&writefastafile(\%h_cdna, $cdnafile, $verbosity);
#     }

#     return \%h_cdna;

# }



sub CreateORFcDNAFromFASTA
{
    my ($fastaFile, $cdnaFile, $orfFile, $nbtx, $minnumtx, $orfType, $verbosity, $kmerMax) = @_;
    # If $nbtx is undef, extract all sequences

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
    die "Your input FASTA '$fastaFile' contains only *$nbseq* sequences.\nNot enough to train the program (default option --ntx|-n)\n" if ($nbseq < $minnumtx);

    # weird have to recreate a seqio object
    $seqin = Bio::SeqIO->new(-file => $fastaFile, -format => "fasta");

    # Go through each sequences
    while(my $seq = $seqin->next_seq()) {

	my $tr = $seq->id();

	$orfFlag = &getTypeOrf2($tr, $seq->seq(), $strand, \%h_orf, $orfType, $kmerMax);

	# Print according to getTypeOrf result
	if ($orfFlag != -1)
	{
	    # Add cDNA only if an ORF is found
	    $h_cdna{$tr} = $seq->seq();
	    print STDERR "\tExtracting ORFs&cDNAs ", $countseqok++,"/$nbtx...\r"  if( defined $nbtx);
	    print STDERR "\tExtracting ORFs&cDNAs ", $countseqok++,"/$nbseq...\r" if(!defined $nbtx);
	}
	else
	{
	    warn "Tx: $tr : ORF is not complete...skipping for training\n" if ($verbosity > 5);
	    next; # next if ORF is not OK
	}

	# Check if nbtx is reached
	if (defined $nbtx && $countseqok == $nbtx)
	{
	    print STDERR "\tMax cDNAs/ORF sequences '$nbtx' reached..ending!\n";
	    last;
	}
    }

    # Final Check if the number of complete ORF is ok
    my $sizehorf = keys(%h_orf);
    die "The number of complete ORF found with computeORF mode is *$sizehorf* ... That's not enough to train the program\n" if (defined $nbtx && $sizehorf < $minnumtx);

    &writefastafile(\%h_orf,  $orfFile, $verbosity);
    &writefastafile(\%h_cdna, $cdnaFile, $verbosity);
}

# sub OldCreateORFcDNAFromFASTA{

#     my  ($fastafile, $cdnafile, $orffile, $nbtx, $orfType, $verbosity)	=	@_;


#     print STDERR "Extract ORF/cDNA from fasta file '$fastafile'..\n";

#     my %h_orf;              # for storing and printing ORF sequence
#     my %h_cdna;             # for storing and printing cDNA sequence
#     my $orfFlag = 0;        # get getTypeOrf result

#     # counter for seq with ORF ok
#     my $countseqok = 0;
#     my $strand     = ".";

#     # Create SeqIO objects
#     my $seqin = Bio::SeqIO->new(-file => $fastafile, -format => "fasta");

#     # count the nb of sequences
#     my $nbseq = 0;
#     $nbseq++ while( my $seq = $seqin->next_seq());
#     die "Your input FASTA '$fastafile' contains only *$nbseq* sequences.\nNot enough to train the program (default option --ntx|-n)\n" if ($nbseq < $minnumtx);

#     # weird have to recreate a seqio object
#     $seqin = Bio::SeqIO->new(-file => $fastafile, -format => "fasta");

#     # Go through each sequences
#     while(my $seq = $seqin->next_seq()) {

# 	my $tr = $seq->id();


# 	# if not orf
# 	if (defined $orffile){ # get also ORF
# 	    $h_cdna{$tr} = $seq->seq();
# 	    $orfFlag     = &getTypeOrf($tr, $seq->seq(), $strand, \%h_orf, $orfType);

# 	    # Print according to getTypeOrf result
# 	    if ($orfFlag != -1){
# 		print STDERR "\tExtracting ORFs&cDNAs ", $countseqok++,"/$nbtx...\r" if( defined $nbtx);
# 		print STDERR "\tExtracting ORFs&cDNAs ", $countseqok++,"...\r"       if(!defined $nbtx);
# 	    } else {
# 		warn "Tx: $tr : ORF is not complete...skipping for training\n" if ($verbosity > 5);
# 		next; # next if ORF is not OK
# 	    }
# 	}
# 	else
# 	{
# 	    $h_cdna{$tr} = $seq->seq();
# 	    print STDERR "\tExtracting cDNAs from FASTA ", $countseqok++,"/$nbtx complete cDNA(s)...\r";
# 	}

# 	# Check if numtx is reached
# 	if (defined $nbtx && $countseqok == $nbtx){
# 	    print STDERR "\tMax cDNAs/ORF sequences '$nbtx' reached..ending!\n";
# 	    last;
# 	}
#     }
#     # if dedfined ORFfile, we write ORF and cDNA file
#     if (defined $orffile){

# 	# Final Check if the number of complete ORF is ok
# 	my $sizehorf = keys(%h_orf);
# 	die "The number of complete ORF found with computeORF mode is *$sizehorf* ... That's not enough to train the program\n" if (defined $nbtx && $sizehorf < $minnumtx);

# 	&writefastafile(\%h_orf,  $orffile, $verbosity);
# 	&writefastafile(\%h_cdna, $cdnafile, $verbosity);

# 	# we write only  cDNA file
#     } else {

# 	my $sizeh = keys(%h_cdna);
# 	die "The number of cDNA sequences is *$sizeh* transcripts... That's not enough to train the program\n" if ($sizeh < $minnumtx);
# 	&writefastafile(\%h_cdna, $cdnafile, $verbosity);
#     }
# }


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

    my ($h, $ref_cDNA_passed, $cdnafile, $orffile, $genome, $nbtx, $minnumtx, $sizecorrec, $orfType, $maxTries, $maxN, $verbosity, $kmerMax, $seed) = @_;

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

    #  hashref tx annotation sizes
    my $refannotsize = ExtractFromHash::getCumulSizeFromGtfHash ($h,$verbosity, 0);

    print STDERR "- Relocate Transcripts \n" if ($verbosity > 0);
    my $i = 0;
    my $h_transcript_size = keys(%{$h});

    my %h_cdna_rdm;  # to store correclty relocated sequences
    my %h_orf;       # to store the orf
    my %h_orf_tmp;   # to temporary store the orf
    my $orfFlag = 0; # to get the type of ORF
    my $cptok   = 0; # to get the number of extracted sequence

    if(defined $seed)
    {
    	srand($seed); # the seed is initiated to have reproducibility
    }

  TX:
    foreach my $tx (sort keys %{$refannotsize}){ # sort for reproducibility

	next if ( ! exists $ref_cDNA_passed->{$tx}); # only keep mRNA tx that are in the cDNA fasat file for sorting CPAT :  for reproducibility


	my $overlap    = 1; # Initialize variable for iterative search for selfoverlap
	my $includeN   = 1; # Initialize variable for iterative search for N
	my $countTries = -1; # Number of tries

	# data for new sequence
	my ($chrrdm, $beg, $end, $seq);
	$seq = "";       # new fasta sequence ==> initialize in case pb with bio::db index
	my $seqORF = ""; # to temporary stock ORF sequence
	my $id;          # transcript name

	if (defined $nbtx && $i == $nbtx){
	    print STDERR "- Max number of transcripts (--nbtx == $nbtx) reached... ending!\n";
	    last;
	}

	# while there is an overlap with known annotation, N in seq or no ORF
	while ($overlap || $includeN || $orfFlag==-1){

	    # maxTries
	    $countTries++;
	    if ( $countTries ==  $maxTries){
		print  STDERR "MaxTries reached ($maxTries) for $tx...skipping it\n";
		next TX;
	    }

	    # my $seed = $i+$countTries; # the seed is initiated accroding to the $i (tx) and the nb of try... if only $i, the same chr:pos will be returned...
	    # # Initialize srand foreach tx
	    # srand($seed);
	    # define a rand indice for all chr hash
	    my $randindex = int( rand(scalar keys %{$refgenomesize}) );
	    my @chrrdm    = sort keys(%{$refgenomesize}); # sort hash for reproducibility
	    $chrrdm       = $chrrdm[$randindex];

	    # define a random start/begin position on the random chr (and thus the end)
	    $beg = int(rand($refgenomesize->{$chrrdm}));
	    $end = $beg + int( $refannotsize->{$tx}->{size} * $sizecorrec);
	    # if the final length is < 200bp (smaller than lncRNA definition), then the size is set to 200
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
	    ($propN,$seq) = getPropN($chrrdm,$beg,$end, $db, 'N');
	    if ($propN == -1){
		warn "Try: $countTries -> Extract sequences for $tx ($chrrdm:$beg-$end) returns an undefined sequence... skipping it\n" if ($verbosity > 10);
	    } elsif ($propN > $maxN){
		warn "Try: $countTries -> Extract sequences for $tx ($chrrdm:$beg-$end) returns a $propN % with N!... skipping it\n" if ($verbosity > 10);
	    }else {
		$includeN = 0;
	    }

	    $id = $tx."_random_($chrrdm:$beg-$end)";
	    # Get ORF sequence
	    $orfFlag = &getTypeOrf2($id, $seq, "+", \%h_orf_tmp, $orfType, $kmerMax);
	    $seqORF = $h_orf_tmp{$id};
	}
	# Write New random sequence
	$h_cdna_rdm{$id} = $seq;
	$h_orf{$id}      = $seqORF;
	print STDERR "\tExtracting ORFs&cDNAs with random coordinate ", $cptok++,"/$nbtx...\r"                      if( defined  $nbtx);
	print STDERR "\tExtracting ORFs&cDNAs with random coordinate ", $cptok++,"/".keys(%{$refannotsize})."...\r" if( !defined $nbtx);

	# verbosity
	$i++;
	if ($verbosity > 0){
	    Utils::showProgress($nbtx, $i, "Print ".$tx.": ");
	}
    }

    my $sizeh = keys(%h_cdna_rdm);
    die "The number of RANDOMLY relocated cDNA sequences =  *$sizeh* transcripts... That's not enough to train the program\n" if ($sizeh < $minnumtx);
    &writefastafile(\%h_cdna_rdm, $cdnafile, $verbosity);
    &writefastafile(\%h_orf,      $orffile,  $verbosity);
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



1;
