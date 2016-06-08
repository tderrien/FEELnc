package RNAshuffle;

$VERSION = v0.0.1;

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use Utils;
use ExtractFromFeature; # getKeyFromFeature line 932
use Bio::DB::Fasta;
use Bio::SeqIO;

# Perl libs
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Data::Dumper;
use List::Util 'shuffle';

# lib directory : ~tderrien/bin/perl/lib/
use Parser;
use ExtractFromHash;
use ExtractFromFeature;
use Intersect;
use Utils;
use Orf;


=encoding UTF-8

=head1 RNAshuffle.pm

=cut


# Function to shuffle the mRNA sequences
#	$seqFile   = the file name containing the sequences to be shuffle
#	$permOut   = the output file for the suffled sequences
#	$nameTmp   = the prefix of temporary files
#	$verbosity = for the verbosity
#	$seed      = the seed to assure reproduciability
sub runFastaUshuffle
{
    my ($seqFile, $permOut, $nameTmp, $verbosity, $seed) = @_;

    # First need to make a one line FASTA file
    my $oneLine = $nameTmp.".oneLine.fa";
    my $flag = 0;
    open SEQFILE,   "$seqFile" or die "Error! Cannot access file '". $seqFile . "': ".$!;
    open ONELINE, "> $oneLine" or die "Error! Cannot access file '". $oneLine . "': ".$!;

    while(<SEQFILE>)
    {
	chop;
	# If first line
	if($flag==0)
	{
	    print ONELINE "$_\n";
	    $flag = 1
	}
	# Else if the line is the name of a sequence
	elsif($_ =~ /^>/)
	{
	    print ONELINE "\n$_\n";
	}
	# Else the line is a sequence line
	else
	{
	    print ONELINE "$_";
	}
    }
    close SEQFILE;
    close ONELINE;
    
    # Use best values get on human gencode
    my $kmerSizePerm = 7;
    my $nbrPerm      = 3;
    my $shufflePath  = Utils::pathProg("fasta_ushuffle");
    my $cmd          = "";
    
    $cmd = "$shufflePath -k $kmerSizePerm -s $seed -n $nbrPerm < $oneLine > $permOut";
    print STDERR "\t\tRun fasta_ushuffle:\n\t\t$cmd.\n" if($verbosity > 1);
    system($cmd);
    if ($? != 0)
    {
        die "\nFailed to run fasta_ushuffle:\n$cmd\n";
    }

    # Remove the one line FASTA
    unlink $oneLine;
}


# Function the check for each suffled sequence if they validate all criteria:
# i) different from the original sequence; ii) have on ORF
#	$codFile    = the coding sequence
#	$nonFile    = where the cDNA from non coding sequences need to be written
#	$nonOrfFile = where the ORF from non coding sequences need to be written
#	$nbtx       = number of transcripts to extract
#	$permOut    = the file with all permutation
#	$orfType    = the ORF type to extract
#	$verbosity  = for the verbosity
#	$kmerMax    = the maximum kmer size
sub getPermutSeq
{
    my ($codFile, $nonFile, $nonOrfFile, $nbtx, $minnumtx, $permOut, $orfType, $verbosity, $kmerMax) = @_;

    # First put in memory the mRNA sequences to check if the permuted sequences are really permuted
    my %codSeq;
    my $seq;
    my $multiFasta = new Bio::SeqIO(-file  => $codFile);

    while($seq = $multiFasta->next_seq())
    {
	$codSeq{$seq->seq()} = 1;
    }

    # Second, for each permuted sequences, check if the are different from mRNA sequences, if they get an ORF
    # until all sequences are passed or if the maximum number is reach
    $multiFasta    = new Bio::SeqIO(-file  => $permOut);
    my $countseqok = 0;
    my $orfFlag    = 0;
    my $nbseq      = 0;
    my $tr         = "";
    my %h_orf;          # for storing and printing ORF sequence
    my %h_cdna;         # for storing and printing cDNA sequence

    # Count the number of sequences if $nbtx is not defined
    if(!defined $nbtx)
    {
	$nbseq++ while( my $seq = $multiFasta->next_seq());
	die "Your permuted FASTA '$permOut' contains only *$nbseq* sequences.\nNot enough to train the program (default option --ntx|-n)\n" if (defined $minnumtx && $nbseq < $minnumtx);
	$multiFasta    = new Bio::SeqIO(-file  => $permOut);
    }

    # Get valid sequences
    while($seq = $multiFasta->next_seq())
    {
	$tr = $seq->id();

	if( exists $codSeq{$seq->seq()} )
	{
	    print STDERR "$tr has not been successfully shuffled... Skip it\n" if ($verbosity > 1);	    
	}
	else
	{

	    # Get the ORF for the sequence
	    $orfFlag = ExtractCdnaOrf::getTypeOrf($tr, $seq->seq(), ".", \%h_orf, $orfType, $kmerMax);

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
	    }

	    # Check if nbtx is reached
	    if (defined $nbtx && $countseqok == $nbtx)
	    {
		print STDERR "\n\tMax cDNAs/ORF sequences '$nbtx' reached.\n";
		last;
	    }
	}
    }

    # Final Check if the number of complete ORF is ok
    my $sizehorf = keys(%h_orf);
    die "The number of complete ORF found with computeORF mode is *$sizehorf*... That's not enough to train the program\n" if (defined $nbtx && $sizehorf < $minnumtx);

    if(defined $nbtx && $countseqok != $nbtx)
    {
	print STDERR "\n\t'$countseqok' shuffled sequences have been successfully validated.\n";
    }    
    if(!defined $nbtx)
    {
	print STDERR "\n\tExtracted '$countseqok' ORF/cDNAs sequences on '$nbseq'.\n";
    }

    # Write output FASTA files
    ExtractCdnaOrf::write2fastafile(\%h_cdna, $nonFile, \%h_orf,  $nonOrfFile, $verbosity);

    # Remove the $permOut file
    unlink $permOut;
}


sub makePerm
{
    my ($codFile, $nonFile, $nonOrfFile, $nbtx, $minnumtx, $orfType, $verbosity, $kmerMax, $nameTmp, $seed) = @_;

    # The file name with all permuted sequences
    my $permOut = $nameTmp.".ushuffle.3perm.fa";

    # Make the permutation
    &runFastaUshuffle($codFile, $permOut, $nameTmp, $verbosity, $seed);

    # Validate and write permutated sequences
    &getPermutSeq($codFile, $nonFile, $nonOrfFile, $nbtx, $minnumtx, $permOut, $orfType, $verbosity, $kmerMax);
}



