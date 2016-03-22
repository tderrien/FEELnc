#!/usr/bin/perl -w


use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Basename;

# database parameters and type are always required
# you  might want to give two others parameters, always going together, gff3 file and gtf file from cuffmerge
use Bio::SeqFeature::database_part;
use Bio::SeqFeature::LncRNAs_Factory;
use Bio::SeqFeature::InteractionCollection;

# Global Variables
my $lncrna_file	="";
my $mrna_file   ="";
my $window      = 10000;  # 10kb
my $maxwindow   = 100000; # 100kb
my $biotype     = 0;
my $verbosity   = 1;
my $help        = 0;
my $man         = 0;


# Parsing parameters
GetOptions(
    "window|w=i"    => \$window,
    "lncrna|i=s"    => \$lncrna_file,
    "mrna|a=s"      => \$mrna_file,
    "maxwindow|m=i" => \$maxwindow,
    "biotype|b"     => \$biotype,
    "verbosity|v=i" => \$verbosity,
    "help|?"        => \$help,
    "man"           => \$man,
    ) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

my $progname=basename($0);

# Test parameters
pod2usage("Error: Cannot read your input lncRNA GTF file '$lncrna_file'...\nFor help, see:\n$progname --help\n") unless( -r $lncrna_file);
pod2usage("Error: Cannot read your input annotation file '$mrna_file'...\nFor help, see:\n$progname --help\n") unless( -r $mrna_file);


unless (defined $maxwindow) {$maxwindow=$window};

my $collection =  Bio::SeqFeature::LncRNAs_Factory->DoItForMe($window,$lncrna_file,$mrna_file,$maxwindow,$biotype,$verbosity);

$collection->print_all_interactions($biotype);
exit(0);




__END__

=pod

=encoding UTF-8

=head1 NAME

FEELnc_classifier.pl - classifying new lncRNAs w.r.t to the annotation of mRNAs

=head1 VERSION

version 0.01

=head1 SYNOPSIS

FEELnc_classifier.pl -i lncRNA.gtf -a mRNA.gtf  > lncRNA_classes.txt

=head1 DESCRIPTION

FEELnc (Fast and Effective Extraction of Long non-coding RNAs) is dedicated to the annotation of lncRNAs
based on a set of transcripts as input (basically a cufflink transcripts.gtf file)
The last step if the pipeline (FEELnc_classifier) consists in classifying new lncRNAs w.r.t to the annotation of mRNAs in order to annotate :

* Intergenic lncRNAs i.e lincRNAs

	- divergent : when the lincRNA is transcribed in an opposite direction (head to head) w.r.t to the closest mRNA

	- convergent: when the lincRNA is transcribed in a convergent direction w.r.t to the closest mRNA

	- same_strand: when the lincRNA is transcribed in a same starnd w.r.t to the closest mRNA

* Genic lncRNAs:  lncRNAs overlapping mRNAs either

	- Exonic:
		antisense : at least one lncRNA exon overlaps in antisense an mRNA exon
		sense : there should not be since there are filtered in the first step

	- Intronic:
		antisense : lncRNA exon overlaps in antisense mRNA introns (but none exons)
		sense : lncRNA exon overlaps in sense mRNA introns (but none exons)

	- Containing:
		antisense : lncRNA intron overlaps antisense mRNA exons
		sense : lncRNA intron overlaps sense mRNA exons

=head1 OPTIONS

=head2 General

  -b, --biotype         Print the biotype of each transcripts in the output
  --help                Print this help
  --man                 Open man page
  -v,--verbosity	Level of verbosity


=head2 Mandatory arguments

  -i,--lncrna=file.gtf	Specify the lncRNA GTF file
  -a,--mrna=file.gtf	Specify the annotation GTF file (file of protein coding annotation)


=head2 Filtering arguments

  -w, --window=10000		Size of the window around the lncRNA to compute interactions/classification [default 10000]
  -m, --maxwindow=100000	Maximal size of the window during the expansion process [default 10000]

=head1 AUTHORS

=over 4

=item -
Valentin Wucher <vwucher@univ-rennes1.fr>

=item
Thomas Derrien <tderrien@univ-rennes1.fr>

=item
Fabrice Legeai <fabrice.legeai@inria.fr>

=back

=head1 COPYRIGHT AND LICENSE

To be done

=cut


