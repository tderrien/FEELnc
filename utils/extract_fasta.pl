#!/usr/bin/perl -w

#
# July 2014
#
# tderrien@univ-rennes1.fr
# extract a fasta sequences based on:
#	- a number range (1-3)
#	- a file of id
#	- ids seperated by commas 
# run 'extract_fasta.pl -h' for help
#
##########################################################################################


############################ library #####################################################
use warnings;
use strict;

use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

use Data::Dumper;
############################# Parameters #################################################


# mfasta file
my $infile	= "";

# Extract 10 first sequences by default
my $range	=	undef ;
my $filter 	=	undef;
my $orf 	= 0;

# help
my $help	=	0;
my $man		=	0;

################################ GetOptions###############################################
GetOptions (
	"infile|i=s" 	=>	\$infile,
	"filter|f=s"	=>	\$filter,
	"orf|o!"		=>	\$orf,	
	"help|?"		=>	\$help,
	"man|m"		=>	\$man,
	
) or pod2usage(2);	

# Print help if needed
pod2usage(1) if $help;
pod2usage(-exitval=>0, -verbose=>2) if $man;

# Mandatory options
###########################################################################################
# test infile
pod2usage ("-i|infile option: I need a non-empty input multifasta file '$infile'...\n") unless (-s $infile);
pod2usage ("-f|filter option: I need a filtering option...\n") unless (defined $filter);

# array of filter
my @filter;
my ($min, $max)= 0,0;

#############################################################################################	
# Check filter

if (defined $filter){

	# filter are in a file
	if (-s $filter){
		open (FILTER, $filter) or die "Cannot open your file of filter '$filter': $!\n";
		
		while (<FILTER>){
			chomp;
			push @filter, split /,/;
		}
		die "Your file of IDs '$filter' should contain Ids separated by commas ','\n" if (@filter<2);
		print STDERR "Extracting sequences by FILE OF IDs:\n", join(",",@filter),"\n";
	
	
	# else filter option is a string
	}elsif (index ($filter, ',') != -1){ # if IDs
		@filter = split /,/ , $filter;
		print STDERR "Extracting sequences by IDs:\n", join(",",@filter),"\n";
	
	
	} elsif (index($filter, '-') != -1) { # if range
		($min, $max) = split /-/, $filter;
		if ($max<$min){my $temp=$max; $max=$min;$min=$temp}
		print STDERR "Extracting sequences by NUMBER from: $min->$max\n";
	
	
	} else {
		print STDERR "Cannot understand filter options : '$filter'\n";
		pod2usage(-verbose => 0);
		exit 1;
	}
}
###############################
# Create SeqIO objects 
my $seqin  = Bio::SeqIO->new(-file => $infile,      -format => "fasta");



# put sequence in array to keep filter order
my @seqArr;
while (my $seq = $seqin->next_seq()) {
	# need to deal with spaces
	#$seq->desc( $seq->id . " ". $seq->desc);
	push(@seqArr, $seq);
}

############################################################

my $seqout = Bio::SeqIO->new( -format => "fasta");


# if filter is ids or file of ids
if (scalar @filter > 0){

	# Keep ordered according to @filter
	foreach my $fid (@filter){
		foreach my $seq (@seqArr) {

			next if ($fid ne $seq->id);

			if ($orf){

				if (checkDnaStartStop($seq->seq())){
					$seqout->write_seq($seq);
				}
			} else {
				$seqout->write_seq($seq);
			}
		}
	}
# if range
}else{
	if($min >0 && $max >0) {
	
		#print "$min $max ok\n";
		my $count=0;
	
		foreach my $seq (@seqArr) {

			$count++;
			if ($orf){

				if (checkDnaStartStop($seq->seq())){
					$seqout->write_seq($seq) if ($count>= $min && $count <=$max);
				}
			} else {
				$seqout->write_seq($seq) if ($count>= $min && $count <=$max);			
			}
		}
	}else {
		print STDERR "Error in extracting sequence...\n";
	}
}

sub checkDnaStartStop {
my ($seq) = @_;
	
	if (defined ($seq) && $seq=~m/^ATG/i && $seq =~m/TAA$|TAG$|TGA$/i){
		return 1;
	} else {
		return 0;
	}
}

__END__

=head1 NAME

extract_fasta.pl - Extract fasta sequences from a multi-fasta file


=head1 SYNOPSIS

perl extract_fasta.pl -i <MultiFasta_File>  -f <filters>

	* Required parameters:
		-i|--infile <MultiFasta_File>
		-f|--filter filter on sequences to extract		[default 'undef']

	* Optional Parameters:

		-o|--orf	get only sequence beginning with start codon and ending with stop.	[default '0']
		-h|--help	brief help message	[default '0']
		-m|--man 	full documentation	[default '0']
		

	* Filter Examples:

$ perl extract_fasta.pl -i infile.fa  B<-f ENSCAFT0001,ENSCAFT0002> 

: will extract sequences ENSCAFT0001 and ENSCAFT0002 from infile.fa

$ perl extract_fasta.pl -i infile.fa  B<-f 1-10>

: will extract 10 first sequences from infile.fa

$ perl extract_fasta.pl -i infile.fa  B<-f infile.txt>

: will extract id sequences in file infile.txt from infile.fa


=head1 Filter parameters

=over 12

=item B<--filter (or -f) = file> 

will extract id sequences in file infile.txt from MultiFasta_File

=item B<--filter (or -f) = range (1-5)> 

will extract 5 first sequences from MultiFasta_File

=item B<--filter (or -f) = ids (ENSCAFTA,ENSCAFTB)> 

will extract sequences ENSCAFT0001 and ENSCAFT0002 from infile.fa

=back

=head1 B<--orf (or -o)>

get only sequence beginning with start codon and ending with stop.	

=head1 B<--help (or -h)>

This Help message

=head1 B<--man (or -m)>

Manual page

