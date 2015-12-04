#!/usr/bin/perl -w

# Perl libs
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Data::Dumper;
use Parallel::ForkManager;


# Own lib : /home/genouest/umr6061/recomgen/tderrien/bin/ThomasPerl/lib/
use Parser;
use ExtractFromHash;
use ExtractFromFeature;
use Intersect;
use Utils;



my $progname=basename($0);

# Variables
my $infile    = '';
my $mRNAfile  = '';
my %biotype;
my $man       = 0;
my $help      = 0;
my $verbosity = 1;


my $minsize      = 200;
my $monoexonic   = 0;  # -1 keep monoexonicAS, 1 keep all monoexonic, 0 remove all monoexonic
                       # restricted by $linconly
my $linconly     = 0;  # bool : 1 to only extract intergenic tx
my $biexonicsize = 25; # minimum size of an exon in bp for transcript having 2 exons
my $minfrac_over = 0;  # min proportion of overlap to remove candidate lncRNAs
my $strandedmode = 1;  # default stranded 
my $proc         = 4;
my $outputlog;

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions(
    'i|infile=s'       => \$infile,
    'a|mRNAfile=s'     => \$mRNAfile,	
    's|size=i'         => \$minsize,
    'biex=i'           => \$biexonicsize,    
    'f|minfrac_over=f' => \$minfrac_over,
    'monoex=i'         => \$monoexonic,
    'l|linconly!'      => \$linconly,    
    'p|proc=i'         => \$proc,
    'b|biotype=s'      => \%biotype,	
    'o|outlog=s'       => \$outputlog,	
    'v|verbosity=i'    => \$verbosity,
    'help|?'           => \$help,
    'man'              => \$man
    ) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

# Test parameters
pod2usage("Error: Cannot read your input GTF file '$infile'...\nFor help, see:\n$progname --help\n") unless( -r $infile);
pod2usage("Error: Cannot read your input annotation file '$mRNAfile'...\nFor help, see:\n$progname --help\n") unless( -r $mRNAfile);
pod2usage ("- Error: \$minfrac_over option '$minfrac_over' should be a float between 0 and 1 [0-1] (e.g 0.5 if 50% overlap is required)\n") unless ($minfrac_over >= 0 and $minfrac_over <= 1);
pod2usage ("- Error: --monoex option '$monoexonic' should be -1 keep monoexonicAS, 1 keep all monoexonic, 0 remove all monoexonic \n") if ($monoexonic != 0 and $monoexonic != 1 and $monoexonic != -1);
pod2usage ("- Error: --monoex option '$monoexonic' (keep monoexonic antisense) cannot be activated in the same time as keep only lincRNA (-l|--linconly)...\n") if ($linconly && $monoexonic == -1);

#############################################################
my $commandline = qx/ps -o args $$/;

# Output log
my $basename = basename($infile);
if (!defined $outputlog){
    ($outputlog = $basename) =~ s/\.[^.]+$/.feelncfilter.log/;
}
open(LOG,">$outputlog") or die("Cannot open '$outputlog'");

print LOG $commandline;
print STDERR "Filtered transcripts will be available in file: '$outputlog'\n";


# Parsing candidate lncRNAs
my $splitchr = 0;
my $reflnc   = Parser::parseGTF($infile, 'exon', $splitchr, undef, $verbosity);


# counters
my ($ctminsize, $ctmonoexonic, $ctdubious) = (0,0,0);

# Filtering steps
foreach my $tx (keys %{$reflnc}){

    my $txfeatures  =   $reflnc->{$tx}->{'feature'};

    # size
    my $size = ExtractFromFeature::features2size($txfeatures, 0);
    if ($size <$minsize){
        print LOG "Filter Size ($minsize): $tx = $size nt...\n";
        $ctminsize++;
        delete $reflnc->{$tx};
        next;
    }
    
    # nb exon
    my $nbexon = ExtractFromFeature::features2nbExon($txfeatures, 0);
    if ($nbexon == 1){ 
    	
    	# if 0                    => we remove all monoexonic tx
    	# if -1 and strand is "." => we also remove the transcript (not strans information)
    	if ($monoexonic == 0 || ($monoexonic == -1 &&  ($reflnc->{$tx}->{'strand'} eq "."))) {
	    print LOG "Filter monoexonic (option $monoexonic): $tx =  $nbexon exon (with strand ",$reflnc->{$tx}->{'strand'},")...\n";
    	    $ctmonoexonic++;
	    delete $reflnc->{$tx};    
	    next;
	}
    }
    
    # dubious biexonic
    my $dubious = ExtractFromFeature::features2dubExon($txfeatures, $biexonicsize, 0);
    if ($dubious && $nbexon == 2){
	print LOG "Filter biexonic ($biexonicsize): $tx =  $nbexon exon...\n";
        $ctdubious++;
        delete $reflnc->{$tx};
        next;
    }
}

print STDERR "> Filter size ($minsize): $ctminsize\n" if ($verbosity > 0);
print STDERR "> Filter monoexonic ($monoexonic): $ctmonoexonic\n" if ($verbosity > 0);
print STDERR "> Filter biexonicsize ($biexonicsize): $ctdubious\n" if ($verbosity > 0);
print STDERR ">> Transcripts left after fitler(s): ",scalar keys (%{$reflnc}),"\n" if ($verbosity > 0);


#########################################################################################
# Split lnc by chr to speed up the overlap loop
my $reflncchr = Parser::splitbyChr($reflnc, $verbosity);

# mRNA ref storing datastructure
my $refmRNAchr = Parser::parseGTF($mRNAfile , 'exon', 1 , \%biotype, $verbosity);



# Launch parralel
my $pm = Parallel::ForkManager->new($proc);


my %refhchild; # will store overlaping lncRNA <=> mRNA
# Sub that will be executed at the end of each child process
$pm -> run_on_finish (
    sub {
	my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $result_ref) = @_;
	my $chr = $ident;

	warn("Child $chr killed by signal $exit_signal"), return if $exit_signal;
	warn("Child $chr exited with error $exit_code"),  return if $exit_code;
	warn("Child $chr encountered an unknown error"),  return if !$result_ref;

	$refhchild{$chr} = $result_ref;
    }
    );


##########
# Overlap
print STDERR "Computing overlap (exon level) with reference annotation...\n" if ($verbosity > 1);
# Process 2 hash chr
foreach my $chrlnc (keys %{$reflncchr} ) {

    my %lnctorm; # hash tmp storing tx IDs to remove, the correct hash is : $refhchild
    
    if (! exists $refmRNAchr->{$chrlnc} ) {
	
	# if user wants monoexonic antisense and the lncRNA chr does not belong to the annotation, we remove all monoexonic the lncRNA associated-chr
	if ($monoexonic == -1){
	    my %h_mono = ExtractFromHash::getMonoExonicFromGtfHash( $reflncchr->{$chrlnc});
	    my @idstorm =  keys %h_mono;
	    print LOG "Filter overlap ($minfrac_over-$monoexonic-$linconly): $_   not antisense...\n" for (@idstorm);
	    delete @{$reflnc}{@idstorm};
    	}
    } elsif (exists $refmRNAchr->{$chrlnc}) { # else we check for overlap
	
	print STDERR "$chrlnc\n";
	# start fork
	my $pid = $pm->start($chrlnc) and next;
     	
     	# get ref on hash per chromosome
        %lnctorm = Intersect::getOverlapping($reflncchr->{$chrlnc}, $refmRNAchr->{$chrlnc}, $strandedmode, $minfrac_over, $monoexonic, $linconly, $verbosity);

    }
    $pm->finish(0, \%lnctorm);
    
}

# wait all sub process
$pm->wait_all_children;


# remove matching transcripts from hash of process chr
foreach my $thread (keys %refhchild){
    my @idstorm =  keys %{$refhchild{$thread}};
    print LOG "Filter overlap ($minfrac_over-$monoexonic-$linconly): $_ overlap ${$refhchild{$thread}}{$_}...\n" for (@idstorm);
    delete @{$reflnc}{@idstorm} ;  
}


print STDERR ">> Transcripts left after overlap: ",scalar keys (%{$reflnc}),"\n" if ($verbosity > 1);
print STDERR "Printing candidates lncRNAs...\n" if ($verbosity > 1);
ExtractFromHash::printGTF($reflnc, 'all',  $verbosity)



__END__

=encoding UTF-8

=pod

=head1 NAME

FEELnc_filter.pl - Extract, filter candidate long non-coding RNAs

version 0.01

=head1 SYNOPSIS

FEELnc_filter.pl -i candidate.gtf -a mRNA.gtf  > candidate_lncRNA.gtf

=head1 DESCRIPTION

FEELnc (Fast and Effective Extraction of Long non-coding RNAs) is dedicated to the annotation of lncRNAs 
based on a set of transcripts as input (basically a cufflink transcripts.gtf file)
The first step if the pipeline (FEELnc_filter) is to filter unwanted/spurious transcripts and/or transcripts
overlapping in sense exons of the reference annotation.

=head1 OPTIONS

=head2 General

    --help			Print this help
    --man			Open man page
    --verbosity			Level of verbosity 0, 1 and 2 [default 1]
  

=head2 Mandatory arguments

    -i,--infile=file.gtf	Specify the GTF file to be filtered (such as a cufflinks transcripts/merged .GTF file)
    -a,--mRNAfile=file.gtf	Specify the annotation GTF file to be filtered on based on sense exon overlap (file of protein coding annotation)
    
=head2 Filtering arguments

    -s,--size=200		Keep transcript with a minimal size (default 200)
    -b,--biotype		Only consider transcript(s) from the reference annotation having this(these) biotype(s) (e.g : -b transcript_biotype=protein_coding,pseudogene) [default undef i.e all transcripts]
    -l,--linconly		Keep only long intergenic/interveaning ncRNAs [default FALSE]
    --monoex=-1|0|1		Keep monoexonic transcript(s): mode to be selected from : -1 keep monoexonic antisense (for RNASeq stranded protocol), 1 keep all monoexonic, 0 remove all monoexonic	[default 0]
    --biex=25			Discard biexonic transcripts having one exon size lower to this value (default 25)
    
=head2 Overlapping specification 


    -f,--minfrac_over=0		Minimal fraction out of the candidate lncRNA size to be considered for overlap [default 0 i.e 1nt]
    -p,--proc=4			Number of thread for computing overlap [default 4]


=head2 Log output

    -o,--outlog=file.log		Specify the log file of output which [default infile.log]



=head1 AUTHORS

=over 4

=item - Valentin WUCHER <vwucher@univ-rennes1.fr>

=item - Thomas DERRIEN <tderrien@univ-rennes1.fr>

=item - Fabrice Legeai <fabrice.legeai@inria.fr>

=back

=head1 COPYRIGHT AND LICENSE

    This software is Copyright (c) 2014 by IGDR - CNRS

=cut
