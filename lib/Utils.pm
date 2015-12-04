package Utils;

$VERSION = v0.0.1;

use strict;
use warnings;
use Scalar::Util;
use File::Basename;
use POSIX;

$| = 1;

# get the minimum of two elements
sub min2 {
    my ($a, $b) = @_;

    if (! defined $a) {
	return $b;
    }
    if (! defined $b) {
	return $a;
    }
    return (($a < $b)? $a: $b);
}

# get the minimum of three elements
sub min3 {
    my ($a, $b, $c) = @_;

    return min2(min2($a, $b), $c);
}

# get the maximum of two elements
sub max2 {
    my ($a, $b) = @_;

    if (! defined $a) {
	return $b;
    }
    if (! defined $b) {
	return $a;
    }
    return (($a > $b)? $a: $b);
}

# get the maximum of three elements
sub max3 {
    my ($a, $b, $c) = @_;

    return max2(max2($a, $b), $c);
}

# get the average of elements
sub average {
    my (@elements) = @_;
    my $sum = 0;
    my $nb = 0;
    foreach my $element (@elements) {
	if (defined $element) {
	    $sum += $element;
	    $nb++;
	}
    }
    if ($nb == 0) {
	return 0;
    }

    return ($sum / $nb);
}

# write to a file
sub writeNewFile {
    my ($fileName, $content) = @_;

    if (-e $fileName) {
	unlink $fileName;
    }
    open (TMPFILE, "+> $fileName") or die $!;
    print TMPFILE $content;
    close TMPFILE;
}

# show a progress bar
sub showProgress {
    my ($aim, $current, $message) = @_;

    if (! defined $message) {
	$message = "";
    }
    my $percent = int ((($current+1) / ($aim+1)) * 100);
    print STDERR $message . " " x (30 - (length $message)) . "[" . "-" x $percent . " " x (100 - $percent) . "]\r";
}

# compute (G+C)%
sub getGC {
    my ($sequence) = @_;
    my $size = 0;
    my $gc = 0;

    for (my $i = 0; $i < length $sequence; $i++) {
	if ((substr $sequence, $i, 1) =~ /[GCgc]/) {
	    $gc++;
	}
	$size++;
    }
    return (($gc / $size) * 100);
}

# remove leading or trailing whitespace characters
sub trim {
    my ($string) = @_;

    $string =~ s/^\s+//;
    $string =~ s/\s+$//;

    return $string;
}

# Sub routines for computing overlap between 2 ranges
# return 1 if overlap, 0 otherwise
sub foverlap
{
    my ($beg1,$end1,$beg2,$end2, $s1, $s2, $mode) = @_;
    $mode //= 0; # mode 0 is unstranded (1 if required same strandness) -1 for different strand

    # test for stranded mode
    if ($mode != 0){
	die 'foverlap: Strand $s1 or $s2 not defined...' if (!defined $s1 || !defined $s2);
	warn "Warning: Mode stranded activated with strandA = '$s1' and strandB = '$s2' ($beg1-$end1 vs $beg2-$end2)\n" if ($s1 ne "+" && $s1 ne "-" && $s2 ne "+" && $s1 ne "-");
    }
    # if overlap
    return ( strandmode($s1,$s2,$mode) and ($end1>=$beg2)&&($beg1<=$end2));


}

# return 1 if the two strands are ok wrt to mode 
# mode  : mode 0 is unstranded // 1 if required same strandness // -1 for different strand
# if strand is not specified "." and mode is stranded or different strand, return 0 (stringent way)
sub strandmode{

    my ($s1,$s2,$mode) = @_;

    # if mode is stranded i.e (1 or -1)
    if ($mode){
	if ($s1 eq "." || $s2 eq "."){
	    return 0;
	}else{
	    if ($mode == 1) {		# if mode is stranded sense
		($s1 eq $s2) ? return 1:0;
	    } elsif ($mode == -1){ # if mode is stranded antisense
		($s1 ne $s2) ? return 1:0;		
	    }
	}
    } else { # if mode is unstranded, return 1
	return 1
    }
}


# Sub routines for computing overlap between 2 ranges
# return size of the overlap if overlap, 0 otherwise
# Note : to be tested extensively
sub foverlapsizerange
{
    my ($b1,$e1,$b2,$e2, $s1, $s2,$stranded) = @_;
    $stranded //= 0; # mode 0 is unstranded (1 if required same strandness)
    
    if ( foverlap ($b1,$e1,$b2,$e2, $s1, $s2, $stranded) ){
	
	if ($b1<=$b2 && $e1<=$e2) {return ($e1-$b2+1)}
	if ($b1>=$b2 && $e1>=$e2) {return ($e2-$b1+1)}
	if ($b2>=$b1 && $e1>=$e2) {return ($e2-$b2+1)}
	if ($b2<=$b1 && $e1<=$e2) {return ($e1-$b1+1)}						
    } else { 
	return 0;
    }
    
}


# Sub routines for computing overlap between 2 ranges
# return the minimal range of the intersection
# Note : to be tested extensively
sub foverlapmin
{
    my ($b1,$e1,$b2,$e2, $s1, $s2, $stranded) = @_;
    $stranded //= 0; 			# mode 0 is unstranded (1 if required same strandness)
    
    if ( foverlap ($b1,$e1,$b2,$e2, $s1, $s2, $stranded) ){
	
	if ($b1<=$b2 && $e1<=$e2) {return ($b2,$e1)}
	if ($b1>=$b2 && $e1>=$e2) {return ($b1,$e2)}
	if ($b1<=$b2 && $e1>=$e2) {return ($b2,$e2)}
	if ($b2<=$b1 && $e1<=$e2) {return ($b1,$e1)}						
    } else { 
	return (0,0);
    }
}



# Sub routines for computing distance between 2 features
# return size of the overlap if overlap, 0 otherwise
# Note : to be tested extensively
sub getdistance
{
    my ($b1,$e1,$b2,$e2) = @_;

    die "getdistance : input value $b1 is not integer\n" unless (  $b1 =~ /^-?\d+$/);
    die "getdistance : input value $e1 is not integer\n" unless (  $e1 =~ /^-?\d+$/);
    die "getdistance : input value $b2 is not integer\n" unless (  $b2 =~ /^-?\d+$/);
    die "getdistance : input value $e2 is not integer\n" unless (  $e2 =~ /^-?\d+$/);
    
    if (($e1>=$b2)&&($b1<=$e2)){return 0}
    else {
	if ($b1>$b2){
	    return ($b1-$e2);
	}else {
	    return ($b2-$e1);
	}
    }
}

# Get unique elements from an array
# http://perlmaven.com/unique-values-in-an-array-in-perl
sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}

# count the nb of lines in a file
# fastest way
# see : http://docstore.mik.ua/orelly/perl/cookbook/ch08_03.htm
sub countlinefile{
    my ($file) = @_;
    my $count = 0;
    open FILE, "$file" or die "Error! Cannot open File ". $file . ": ".$!;
    $count += tr/\n/\n/ while sysread(FILE, $_, 2 ** 16);
    close FILE;
    return $count;
}

#  guess format based on file suffix
sub guess_format {
    my ($filename) = @_;
    
    return 'fasta'   if ($filename =~ /\.fasta$/i || $filename =~ /\.fa$/i);
    return 'gff'     if ($filename =~ /\.gff3?$/i || $filename =~ /\.gff?$/i);
    return 'gtf'     if ($filename =~ /\.gtf$/i  || $filename =~ /\.gff2?$/i);
    return 'bed'     if ($filename =~ /\.bed$/i);

    return 'gtf'; #the default
}

sub renamefile {

    my ($filename, $suffix) = @_;
    my $basename = basename($filename);
    my $newname;
    
    ($newname = $basename) =~ s/\.[^.]+$/$suffix/;

    return $newname;
}

# Test if a program is in PATH and exec
# In : program name
# Out : full path to the program
# Example :     my $pathlogit       =   Utils::pathProg("make_logitModel.py");
sub pathProg{
    my ($prog) = @_;
    $prog    //= undef;
    my @path;

    if (defined $prog){
        @path = grep { -x "$_/".$prog}split /:/,$ENV{PATH};

	if (!@path){
	    die "pathProg: '$prog' not in your PATH...";
        }		

        if ( -x $path[-1]."/".$prog){
            return $path[-1]."/".$prog;
        } elsif (-x $prog){
            return $prog;
        } else {
	    die "pathProg: '$prog' not executable..";
        }
    }else{
	die "pathProg: undefined program '$prog'...\n";
    }
}



# getFastaNbr: return the number of sequence in a FASTA file
sub getFastaNbr
{
    my($input) = @_;

    # Open input file
    my $inputFasta = new Bio::SeqIO(-file => "$input", '-format' => 'fasta');
    my $seq;
    my $cpt = 0;

    while($seq = $inputFasta->next_seq())
    {
	$cpt++;
    }

    return($cpt);
}


# divFasta: cut a FASTA file in two files, with the first one equal to X% of the sequences of the input file (default 0.5%)
sub divFasta
{
    my($input, $name1, $name2, $perc, $verbosity) = @_;
    $verbosity //= 1;

    my $nbrSeq = &getFastaNbr($input);
    my $nbr1   = POSIX::floor($perc * $nbrSeq);
    my $nbr2   = $nbrSeq - $nbr1;

    # Open input file
    my $inputFasta = new Bio::SeqIO(-file => "$input", '-format' => 'fasta');
    my $seq;
    
    # Open output files
    my $out1 = new Bio::SeqIO(-file => "> $name1" , '-format' => 'fasta');
    my $out2 = new Bio::SeqIO(-file => "> $name2" , '-format' => 'fasta');

    print STDERR "\t\tDividing $input ($nbrSeq seq) file into two learning files ($nbr1 and $nbr2 seq)\n" if($verbosity > 1);
    # Write the first $nbr1 sequences into $out1
    for(my $i=1; $i<=$nbr1; $i++)
    {
	$seq = $inputFasta->next_seq();
	$out1->write_seq($seq);
    }
    
    # Write the last $nbr2 sequences into $out2
    for(my $i=1; $i<=$nbr2; $i++)
    {
	$seq = $inputFasta->next_seq();
	$out2->write_seq($seq);
    }
}


1;
