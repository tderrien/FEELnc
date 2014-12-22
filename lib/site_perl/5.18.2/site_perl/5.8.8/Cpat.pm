package Cpat;

$VERSION = v0.0.1;

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use List::Util 'shuffle';

use Utils;
use ExtractFromFeature; # getKeyFromFeature line 932
use Bio::DB::Fasta;






# Uses all commande of cpat and a small script r
# Input: 
#	- fasta file of ORF gene($trainORF)
#	- fasta file of mRNA gene($trainmRNA)
#	- fasta file of non coding gene ($trainlnc)
#	- name of the ouputfile ($outfile) 
#   - verbosity
# Ouput:
#	- result of cpat

sub runCPAT{
	my ($trainORF, $trainmRNA, $trainlnc, $outfile, $testfile, $verbosity) = @_;
    $trainORF   ||= undef;
    $trainmRNA  ||= undef;
    $trainlnc   ||= undef;
    $testfile   ||= undef; # if testfile is still undef, we only build the hexamer and logit model (not CPAT) in order to find the best cutoff
    $outfile    ||= "cpat_out";

    # test path program
    my $pathhexamer     =   Utils::pathProg("make_hexamer_tab.py");
    my $pathlogit       =   Utils::pathProg("make_logitModel.py");
    my $pathcpat        =   Utils::pathProg("cpat.py");
            
    # test infiles
    die "runCPAT: ORF training file not defined... exiting\n"    if (!defined $trainORF);
    die "runCPAT: mRNA training file not defined... exiting\n"   if (!defined $trainmRNA);
    die "runCPAT: lncRNA training file not defined... exiting\n" if (!defined $trainlnc);
    
    # emtpy
    die "runCPAT: ORF training file '$trainORF' is empty... exiting\n"    unless (-s $trainORF);
    die "runCPAT: mRNA training file '$trainmRNA' is empty... exiting\n"  unless (-s $trainmRNA);
    die "runCPAT: lncRNA training file '$trainlnc' is empty... exiting\n" unless (-s $trainlnc);
    
    
    ################
    # build the hexamer table
	print "CPAT:: Build hexamer table\n" if ($verbosity > 5);	
	my $cmd = "$pathhexamer -c $trainORF -n $trainlnc  > ".$outfile."_hexamer.table  2>/dev/null";
	print "$cmd\n" if ($verbosity > 10);
	system($cmd);
	
	# test that outfile are created
	die "Something went wrong with CPAT::$pathhexamer\n Check CPAT installation..."  unless (-s $outfile."_hexamer.table");

	# build the logit
	print "CPAT:: Build Logit model \n" if ($verbosity > 5);
	$cmd    = "$pathlogit -x ".$outfile."_hexamer.table -c $trainmRNA -n $trainlnc -o $outfile  2>/dev/null";
	print "$cmd\n" if ($verbosity > 10);
	system($cmd);
	
	# test that outfile are created
	die "Something went wrong with CPAT::$pathlogit\n Check CPAT installation..." unless ( -r $outfile.".logit.RData");

	
	# use CPAT if testfile if defined
	# els
	if (defined $testfile){
		if (-s $testfile){
			print "CPAT:: Launch CPAT\n" if ($verbosity > 5);
			$cmd    = "$pathcpat -d ".$outfile.".logit.RData -x ".$outfile."_hexamer.table -g $testfile -o ".$outfile.".cpat  2>/dev/null";
			print "$cmd\n" if ($verbosity > 10);
			system($cmd);
		} else {
			die "runCPAT: cpat test file '$testfile' is empty... exiting\n";
		}
	}
	# cleaning
	unlink $outfile.".logit.RData", $outfile."_hexamer.table", $outfile.".cpat.dat", $outfile.".make_logitModel.r", $outfile.".cpat.r";
		
}

# From a hashref of CPAT 
# write a randomize file for 10 fold cross validation

sub WriteRdmCPATFile{

	my ($h, $outfile, $verbosity) = @_;

	# Create a tmp file which shuffle lncRNA and mRNA used to compute 10 cross validation Sn Sp by R
	open(CPATVAL,"> ".$outfile) or die("Cannot open tmpfile for CPAT validation\n");
	
	# header
	print CPATVAL "\tmRNA_size\tORF_size\tFickett_score\tHexamer_score\tcoding_prob\tlabel\n";
	
	foreach my $id ( shuffle  keys (%$h) ){
		my @ar = ($id,
				 $h->{$id}->{'mRNA_size'}, 
				 $h->{$id}->{'ORF_size'}, 
				 $h->{$id}->{'Fickett_score'}, 
				 $h->{$id}->{'Hexamer_score'}, 
				 $h->{$id}->{'coding_prob'},
				 $h->{$id}->{'label'}
				 );
		print CPATVAL join ("\t", @ar);
		print CPATVAL "\n";
				 
	}	
	close CPATVAL;
	
	return 1;
}

sub permute {
# random permutation of array's elements
# only the "size" last elements are randomized

    my $array = shift;
    my $size = shift;
    my $i;
    my $j;

    for ($i = $size; $i > 1; $i--) {
        $j = 1 + int(rand(1) * $i);
        @$array[$i,$j] = @$array[$j,$i];
    }
}


1;