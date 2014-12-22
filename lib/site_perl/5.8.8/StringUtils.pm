package StringUtils;

$VERSION = v0.0.1;

use strict;
use warnings;
use Bio::DB::Fasta;
use Bio::Seq;
use File::Basename;

use Data::Dumper;


use Utils;


$| = 1;


# Get sequence using Bio::DB::Fasta
sub getSubSequenceFasta{
	my ($genome, $chr, $start, $end, $strand) = @_;
	
	# Create index
	my $db       = Bio::DB::Fasta->new($genome);
	
	my $seqstr   = $db->seq($chr, $start => $end);

	if ($strand eq "-"){
		$seqstr = getRevComp($seqstr);
	}
	
	return $seqstr;

}


# get reverseComplement of a sequence
sub getRevComp {

	my ($sequence, $verbosity)	= @_; 
	warn "Utils::getRevComp: Your input sequence is empty\n" if ($sequence eq "" && $verbosity > 5);

	$sequence 	=~tr/ACGTacgtNn/TGCAtgcaNn/;
	$sequence 	=  reverse($sequence);
	return ($sequence);	

}

# get a subsequence in a big multi-fasta file (faster than 
# uses Heng Li samtools faidx : see http://samtools.sourceforge.net/samtools.shtml
sub getSubSequenceSamtools{

	my ($fileName, $sequence_id, $start, $end) = @_;
	
	my $seq_return="";
	
	# Test if samtools installed and is executable by current user
	if (! grep { -x "$_/samtools"}split /:/,$ENV{PATH}){
	
		print STDERR "Error in getSubSequenceSamtools : You need install and put in your PATH the 'samtools' programs \n";
		print STDERR "See here: http://samtools.sourceforge.net/samtools.shtml\n";
		exit;
	}
	
	# open FASTA File
	open FILE, "$fileName" or die "Argh, Cannot open FASTA file '" . $fileName . "': " . $!;
	
	# test if index
	print STDERR "No index found for $fileName!\nBuilding it with samtools:\n" if (! -r $fileName.".fai");
	# Lauch samtools and remove header using grep -v
	# Note :
	# samtools faidx first tests whether there is an index in the same directory with .fai suffix
	# if not, it creates it with message "[fai_load] build FASTA index." and then run the request
	open(SAMTOOLS,"samtools faidx $fileName $sequence_id:$start-$end | grep -v '>' |");
	while(my $line = <SAMTOOLS>){
	     chomp ($line);
	     $seq_return .= $line;
	}
	# test if empty
	if ($seq_return eq "") {
		print STDERR "Error : Empty sub-sequence (call to getSubSequenceSamtools $fileName, $sequence_id, $start, $end)!\n";
		print STDERR "Check that start <= end coordinates and chromosome id match...\n";
	    exit;
	}
	return $seq_return;

}
        
# get a subsequence in a (Multi) Fasta file
sub getSubSequence {
  my ($fileName, $sequence, $start, $end) = @_;

  my $goodRegion = 0;
  my $newSequence = "";
  my $line;
  my $count = 0;

  open FILE, "$fileName" or die "Cannot open FASTA file '" . $fileName . "': " . $!;

  if ($end < $start) {
    my $tmp = $start;
    $start = $end;
    $end = $tmp;
  }

  $start--;

  $line = <FILE>;
  if ($line =~ /$sequence/) {
    my $pos1 = tell(FILE);
    $line = <FILE>;
    my $pos2 = tell(FILE);
    chomp $line;
    my $size = length $line;
    my $address = $pos1 + (($start - ($start % $size)) / $size) * ($pos2 - $pos1);

    $goodRegion = 1;
    $count      = Utils::max2(0, $start - ($start % $size));
    seek FILE, $address, 0 or die "Cannot seek FASTA file '" . $fileName . "' at adress " . $address . ": $!";
  }
  else {
    seek(FILE, 0, 0);
  }

  while ($line = <FILE>) {
    chomp $line;

    if ($goodRegion) {
      if ($line =~ /\>/) {
        close FILE;
        return $newSequence;
      }
      
      my $subStart = $start - $count;
      if ($subStart < 0) {
        $subStart = 0;
      }
      my $subEnd = $end - $count;
      my $subSize = $subEnd - $subStart;
      if ($subSize + $subStart > length($line)) {
        $subSize = length($line) - $subStart;
      }
      if ($subEnd < 0) {
        close FILE;
        return $newSequence;
      }
      if ($subStart <= length($line)) {
        $newSequence .= substr $line, $subStart, $subSize;
      }

      $count += length($line);
    }
    
    if ($line =~ /$sequence/) {
      $goodRegion = 1;
    }
  }
  close FILE;
 
  if ($newSequence eq "") {
    die "Error: Empty sub-sequence (call to getSubSequence $fileName, $sequence, $start, $end)!\n";
  }
  return $newSequence;
}





1;