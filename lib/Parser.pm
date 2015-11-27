package Parser;

$VERSION = v0.0.1;

use strict;
use warnings;
use File::Basename;
use Data::Dumper;

use ExtractFromHash;
use Utils;

=head1 NAME

package B<Parser> - Parse input files for FEELnc datastructures

=head1 SYNOPSIS

	use Parser;
	my $refhgtf  = Parser::parseGTF('infile.gtf');
	my $refhgene = Parser::parseGTFgnlight('infile.gtf');

=head1 DESCRIPTION

This module provides different function to parse input files (especially gtf)
and return specific data structure: 3 data structures (dsc) are possible:
	* tx-based   : transcript-based with exon as arrayfeat
	* chrom-based: chr-based (see split variable)
	* gene-based : gene-based with transcript and exon levels (possibly with only gene ranges (see parseGTFgnlight))
                        - full
                        - light
=cut

=head2 parseGTF
	Title   : parseGTF
	Function: parse a .GTF file
	Example : $refmrna = Parser::parseGTF($infilegtf);
	Returns : A Feelnc tx-based DSC
	Args    :
                   - infile    : file - a gtf file (mandatory)
                   - levels    : string - specify the the level(s) of annotation to be parsed separated by "," (e.g exon,CDS)
                   - split     : boolean - 0 or 1 if we want to return a chrom-based dsc
                   - refhfilter: hashref - with key value to be extracted e.g 'transcript_id=ENSCAFT00000017943,ENSCAFT00000017911'
                   - verbosity : numeric - level of verbosity

=cut
sub parseGTF{
    my ($infile, $levels, $split, $refhfilter, $verbosity) = @_;
    $infile     //= '';
    $levels     //= 'exon';
    $split      //= 0; 	   # split data structure by chr
    $refhfilter //= undef; # filter some interesting tag or fields (beyond gene_id and transcript_id) for instance 'transcript_biotype' : use transcript_id,gene_id for just 12th fields
    # in a program with hash (get)option such as $> -f transcript_id=ENSCAFT00000017943,ENSCAFT00000017911
    #  this hashref will be : 'transcript_id' => 'ENSCAFT00000017943,ENSCAFT00000017911'
    $verbosity	//= 1;
    
    # infile
    my $basename  = basename($infile);
    my @all_level = split (",", $levels);

    # Open file
    open GTFFILE, "$infile" or die "Error! Cannot open GTF File ". $infile . ": ".$!;

    print STDERR "Parsing file '$basename'...\n";

    # count nb of lines
    my $lc=Utils::countlinefile($infile);

    # Store data
    my %h_transcript;

    # iterator line
    my $i = 0;

    # label LINE for next if filterFields activated
  LINE:
    while (<GTFFILE>){

	# verbose
	$i++;
	Utils::showProgress($lc, $i, "Parse input file: ") if ($verbosity > 0);

	chop;
	next if /^#/;
        next if /^$/;
        next if /^track/;

	# split line by tab
	my ($chr, $source, $level, $beg_feat, $end_feat, $score, $strand, $frame, $attributes) = split(/\t/);

	# check number of columns
	die "Error:\n $_ does not look like GTF... not 9 columns\n Check that fields are tab-separated.\n" if ( !defined $attributes );

	next if ($chr =~ /PATCH$/ ||  $chr =~ /TEST$/);
        # gene level line does not contain transcript_id so we removed them at the moment
	next if ($level eq "gene");

	# Patch for wrong biotype in field 2 (i.e Ensembl gtf file)
        my $correctbiotype = correctSourceBiotype($source, $i, $verbosity);
	$attributes        = $attributes.$correctbiotype if (defined $correctbiotype);


	# parse attributes and extract interesting information (filtertag)
	my $refattrib = parseAttributes($attributes, $refhfilter, $verbosity);
	next unless (defined $refattrib);

	# check for mandatory tags : gene_id and transcript_id
	die "Error:\n[$_] does not contain 'transcript_id' attribute...\n" unless(defined $refattrib->{'transcript_id'});
	die "Error:\n[$_] does not contain 'gene_id' attribute...\n"       unless(defined $refattrib->{'gene_id'});

	my $transcript = $refattrib->{'transcript_id'};
	my $gene_id    = $refattrib->{'gene_id'};

	# delete from ref since we already check they exist and we dont need them anymore
	delete $refattrib->{'transcript_id'};
	delete $refattrib->{'gene_id'};

	for my $lvl (@all_level) {

	    if ($lvl eq $level){

		$h_transcript{$transcript}->{"chr"}     = $chr;
		$h_transcript{$transcript}->{"source"}  = $source;
		$h_transcript{$transcript}->{"startt"}  = Utils::min2($h_transcript{$transcript}->{"startt"}, $beg_feat);
		$h_transcript{$transcript}->{"endt"}    = Utils::max2($h_transcript{$transcript}->{"endt"}, $end_feat);
		$h_transcript{$transcript}->{"score"}   = $score;
		$h_transcript{$transcript}->{"strand"}  = $strand;
		$h_transcript{$transcript}->{"gene_id"} = $gene_id;


		# Feature infos are stored in N hashes from a array reference
		# $h{$tr}->{"feature"} is the reference on a array
		# @{$h{$tr}->{"feature"}} to dereference it
		# To obtain a particular element in the array : ${ $h{$tr}->{"feature"} }[0]
		my %feature = ( "start" => $beg_feat, "end" => $end_feat, "feat_level" => $level, 'strand' => $strand, 'frame' => $frame);

		# add refeature
		my %merge = (%{$refattrib}, %feature);

		push (@{$h_transcript{$transcript}->{"feature"}}, \%merge);

		# Filter *official fields* with filterFields (-f option) on $refhfilter hashref
		if ( defined  $refhfilter && !filterFields($refhfilter, $h_transcript{$transcript}) ){
		    delete $h_transcript{$transcript};
		}
	    }
	}
    }
    print STDERR "\n" if ($verbosity > 0); # because of the showProgress

    # close handle
    close GTFFILE;

    # Test parsing i.e empty hash
    if (scalar keys(%h_transcript) == 0){
	print STDERR "Parser::parseGTF => Data Structure returns an empty hash\nPossible reasons:\n";
	print STDERR "\t*Feature level '$levels' is not present in 3rd field of '$infile'\n";
	print STDERR "\t*chromosome/seqname (chr) or patch chromosome...\n";
	print STDERR "\t*Filtering tag/Attributes (--filter|-f) option returns no results\nTry --help for help\n";
	exit;
    }
    # Sort exons
    my $refh = ExtractFromHash::sortExons(\%h_transcript, $verbosity);

    # Do we split by chr?
    if ($split){
	my $refhchr = splitbyChr($refh, $verbosity);

	return $refhchr;

    } else {
	return $refh;
    }
}

=head2 parseAttributes

    Title   : parseAttributes
  Function: parse attributes of a gtf file (i.e after 9th column)
    Returns : A hashref with key/val corresponding to the attributes
    Example : $href	= Parser::parseAttributes('gene_id "RLOC_00032935"; transcript_id "CFRNASEQ_IGNC_Spliced_00193147"; transcript_biotype "lncRNA";');
Args		:
    - string	: string - corresponding to attrributes (mandatory)
    - refhfilter: hashref - with key value to be extracted e.g 'transcript_id=ENSCAFT00000017943,ENSCAFT00000017911'
    - verbosity	: numeric - level of verbosity
  Note: transcript_biotype and gene_biotype (not transcript/gene_type) are the official nomenclature for 'type' of transcript (lincRNA, pseudogene...)

=cut
# parse gtf attributes after transcript_id
sub parseAttributes{

    my ($string, $refhfilter, $verbosity) = @_;
    $verbosity  //= 1;
    $refhfilter //= undef; # filter some interesting (beyond gene_id and transcript_id) for instance 'transcript_biotype' : use transcript_id,gene_id for just 12th fields
                           # in a program with hash (get)option such as $> -f transcript_id=ENSCAFT00000017943,ENSCAFT00000017911
                           # this hashref will be : 'transcript_id' => 'ENSCAFT00000017943,ENSCAFT00000017911'

    my %attribs; # hash that will be returned


    die "parseAttributes: attributes field '$string' not defined\n" unless defined $string;

    my @extrafields = split(";", $string);


    # store ids and additional information in second hash
    foreach my $attr ( @extrafields ) {

	next unless $attr =~ /^\s*(.+)\s(.+)$/;
	my $c_type  = $1;
	my $c_value = $2;
	$c_value=~ s/\"//g;
	if (exists($attribs{$c_type})){
	    print STDERR "WARNINGS: key '$c_type' already exists with ",$attribs{$c_type}," value...\n" if ($verbosity > 1);
	}

	# **** FEELnc nomenclature is to use : 'biotype' instead of 'type' *******
	$c_type =~ s/transcript_type/transcript_biotype/g;
	$c_type =~ s/gene_type/gene_biotype/g;


	# Filtering
	if (defined $refhfilter && exists ($refhfilter->{$c_type}) ){ # if exists tag key in line attrib

	    my @valfitler = split /,/, $refhfilter->{$c_type};
	    my %hash = map { $_ => 1 } @valfitler;

	    return undef if ( !exists $hash{$c_value} ); # we'll remove this line

	}

	# populate attribs
	$attribs{$c_type} = $c_value;
    }
    return \%attribs;
}




=head2 filterFields
    Title   : filterFields
    Function: parse attributes of a gtf file (i.e after 9th column)
    Returns : A boolean : 0 if we filter out the transcript and 1 otherwise
    Example : $bool = Parser::filterFields(\strand => '-', $h_transcript{$transcript}
Args		:
    - refhfilter: hashref - with key value to be extracted e.g 'transcript_id=ENSCAFT00000017943,ENSCAFT00000017911'
    - refhanonym: hasref - a tx-based or gene-based refh
    Note:  see also parseAttributes for attributes
    Note:  to check if to be deprecated
=cut
sub filterFields{

    my ($refhfilter, $refhanonym) = @_;
    $refhfilter //= undef; # filter some interesting fields for instance 'strand'
                           # in a program with hash (get)option such as $> -f strand=-,+
                           # this hashref will be : 'chr' => 'chrX,chr1'

    die "parseGtf::filterFields => fields '$refhfilter' not defined\n" unless(ref $refhfilter);
    die "parseGtf::filterFields => fields '$refhanonym' not defined\n" unless(ref $refhanonym);


    foreach my $k_refhanonym (keys %{$refhanonym}){
	if (exists ($refhfilter->{$k_refhanonym})){

	    my @valfitler = split /,/, $refhfilter->{$k_refhanonym};
	    my %hash = map { $_ => 1 } @valfitler;

	    return undef unless ( exists($hash{$refhanonym->{$k_refhanonym}})  );
	}
    }

    return 1;
}




=head2 correctSourceBiotype
    Title   : correctSourceBiotype
    Function: Correct the transcript_biotype if it is in 2nd fields of the gtf file (such as in Ensembl)
    Returns : a string to be added to the attribtues or undef otherwise : my $string " transcript_biotype \"$source\";";
    Example : my $correctbiotype = correctSourceBiotype($source, $i, $verbosity);
    Args:
          - source   : string - corresponding to 2nd fields of the gtf
          - i        : numeric - $i is just a counter for printing (efault is to print 5 occurences)
          - verbosity: numeric - level of verbosity
    Note: checks whether the 2nd columns contains known biotype extracted from Ensembl Homo_sapiens.GRCh37.70.gtf file (see list in the code)
=cut
sub correctSourceBiotype{

    my ($source, $i, $verbosity) = @_; # $i is just a counter for printing

    my $string = undef;

    # known tx_biotype were extracted from human Ens70 .gtf with: cat ~/dogomaha/DATA/hg19/annotation/Homo_sapiens.GRCh37.70.gtf | awk '{print $2}' | sort | uniq | transp.awk
    my @known_biot = qw/3prime_overlapping_ncrna ambiguous_orf antisense IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene lincRNA miRNA misc_RNA Mt_rRNA Mt_tRNA non_coding nonsense_mediated_decay non_stop_decay polymorphic_pseudogene processed_pseudogene processed_transcript protein_coding pseudogene retained_intron rRNA sense_intronic sense_overlapping snoRNA snRNA TEC transcribed_processed_pseudogene transcribed_unprocessed_pseudogene TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene TR_V_pseudogene unitary_pseudogene unprocessed_pseudogene/;
    my %hash = map { $_ => 1 } @known_biot;

    # 	if (grep {$_ =~ m/^\Q$source\E$/i} @known_biot){
    if (exists ($hash{$source})){
	$string = " transcript_biotype \"$source\";";
	print STDERR "WARNING : Field 2 of the input file matches a 'transcript_biotype' = '$source' (should be the source of the file)\n" if ($verbosity > 1 && $i < 5);
    }
    return $string;
}


=head2 splitbyChr
    Title   : splitbyChr
    Function: split a tx/gene-based feelnc dsc into a chrom-based (seqname) feelnc dsc to speed up overlap computation
    Returns : a hashref chrom-based feelnc dsc
    Example : my $refhchr = splitbyChr($refh, $verbosity);
    Args    :
              - hreftx   : href - tx/gene-based feelnc dsc
              - verbosity: numeric - level of verbosity
=cut
sub splitbyChr{

    my ($hreftx, $verbosity) = @_;
    $verbosity //= 1;

    die "splitbyChr: hash not defined" unless(ref($hreftx));

    print STDERR "Splitting transcripts by chromosomes\n" if ($verbosity > 1);

    my %h_chr;
    foreach my $tr (keys %{$hreftx}){

	my $chr = $hreftx->{$tr}->{"chr"};
	$h_chr{$chr}->{$tr} = $hreftx->{$tr};

    }

    # return ref hash split
    return \%h_chr;
}


=head2 parseGTFgene
    Title   : parseGTFgene
    Function: parse a .GTF file into gene-based feelnc dsc
    Example : $refmrna = Parser::parseGTF($infilegtf);
    Returns : A Feelnc tx-based DSC
    Args    :
    - infile    : file - a gtf file (mandatory)
    - levels    : string - specify the the level(s) of annotation to be parsed separated by "," (e.g exon,CDS)
    - split     : boolean - 0 or 1 if we want to return a chrom-based dsc
    - refhfilter: hashref - with key value to be extracted e.g 'transcript_id=ENSCAFT00000017943,ENSCAFT00000017911'
    - verbosity : numeric - level of verbosity

=cut
sub parseGTFgene{

    my ($infile, $levels, $split, $refhfilter, $verbosity) = @_;
    $infile     //='';
    $levels     //='exon';
    $split      //= 0; 	   # split data structure by chr
    $refhfilter //= undef; #  filter some interesting (beyond gene_id and transcript_id) for instance 'transcript_biotype' : use transcript_id,gene_id for just 12th fields
    $verbosity  //= 1;

    my $basename  = basename($infile);
    my @all_level = split (",", $levels);

    # Open file
    open GTFFILE, "$infile" or die "Error! Cannot open GTF File ". $infile . ": ".$!;

    print STDERR "Parsing parseGTFgene file '$basename'...\n";

    # count nb of lines
    my $lc=Utils::countlinefile($infile);

    # Store data
    my %h_gene;

    # counter line
    my $i = 0;

    # label LINE for next if filterFields activated
  LINE:
    while (<GTFFILE>){

	# verbose
	$i++;
	Utils::showProgress($lc, $i, "Parse input file: ") if ($verbosity > 0);

	#
	chomp;
	next if /^#/;
        next if /^$/;
	next if /^track/;

	# split by tab
	my ($chr, $source, $level, $beg_feat, $end_feat, $score, $strand, $frame, $attributes) = split(/\t/);


	# check number of columns
	die "Error:\n[$_] does not look like GTF... not 9 columns\n Check that fields are tab-separated.\n" if ( !defined $attributes );

	# Patch for wrong biotype in field 2 (i.e Ensembl gtf file)
        my $correctbiotype = correctSourceBiotype($source, $i, $verbosity);
	$attributes        = $attributes.$correctbiotype if (defined $correctbiotype);


	# parse attributes and extract interesting information (filtertag)
	my $refattrib = parseAttributes($attributes, $refhfilter, $verbosity);
	next unless (defined $refattrib);


	# check for mandatory tags : gene_id and transcript_id
	die "[$_] does not contain 'transcript_id' attribute...\n" if (!defined $refattrib->{'transcript_id'} && $level ne "gene");
	die "[$_] does not contain 'gene_id' attribute...\n" unless(defined $refattrib->{'gene_id'});

	my $transcript = $refattrib->{'transcript_id'};
	my $gene_id    = $refattrib->{'gene_id'};

	# delete from ref since we already check they exist and we dont need them anymore
	delete $refattrib->{'transcript_id'};
	delete $refattrib->{'gene_id'};

	for my $lvl (@all_level) {

	    if ($lvl eq $level){

		# gene level
		$h_gene{$gene_id}->{"chr"}    = $chr;
		$h_gene{$gene_id}->{"source"} = $source;
		$h_gene{$gene_id}->{"startg"} = Utils::min2($h_gene{$gene_id}->{"startg"}, $beg_feat);
		$h_gene{$gene_id}->{"endg"}   = Utils::max2($h_gene{$gene_id}->{"endg"}, $end_feat);
		$h_gene{$gene_id}->{"score"}  = $score;
		$h_gene{$gene_id}->{"strand"} = $strand;

		# Filter *official fields* with filterFields (-f option) on $refhfilter hashref
		if ( defined  $refhfilter && !filterFields($refhfilter, $h_gene{$gene_id}) ){
		    delete $h_gene{$gene_id};
		    next LINE;
		}

		# Patch for wrong biotype in field 2 (i.e Ensembl gtf file)
		my $hrefBiotype = correctSourceBiotype($source, $i, $verbosity);


		# tx level
		$h_gene{$gene_id}->{"transcript_id"}->{$transcript}->{"chr"}    = $chr;
		$h_gene{$gene_id}->{"transcript_id"}->{$transcript}->{"source"} = $source;
		$h_gene{$gene_id}->{"transcript_id"}->{$transcript}->{"startt"} = Utils::min2($h_gene{$gene_id}->{"transcript_id"}->{$transcript}->{"startt"}, $beg_feat);
		$h_gene{$gene_id}->{"transcript_id"}->{$transcript}->{"endt"}   = Utils::max2($h_gene{$gene_id}->{"transcript_id"}->{$transcript}->{"endt"}, $end_feat);
		$h_gene{$gene_id}->{"transcript_id"}->{$transcript}->{"score"}  = $score;
		$h_gene{$gene_id}->{"transcript_id"}->{$transcript}->{"strand"} = $strand;


		# Feature infos are stored in N hashes from a array reference
		# $h{$tr}->{"feature"} is the reference on a array
		# @{$h{$tr}->{"feature"}} to dereference it
		# To obtain a particular element in the array : ${ $h{$tr}->{"feature"} }[0]
		my %feature = ( "start" => $beg_feat, "end" => $end_feat, "feat_level" => $level, 'strand' => $strand, 'frame' => $frame);

		# add refeature
		my %merge = (%{$refattrib}, %feature);

                # Filter *official fields* with filterFields (-f option) on $refhfilter hashref
		if ( defined  $refhfilter && !filterFields($refhfilter, $h_gene{$gene_id}) ){
		    delete $h_gene{$gene_id};
		}
		push (@{$h_gene{$gene_id}->{"transcript_id"}->{$transcript}->{"feature"}}, \%merge);
	    }
	}
    }
    print STDERR "\n" if ($verbosity > 0); # because of the showProgress


    # close handle
    close GTFFILE;

    # Test parsing i.e empty hash
    if (scalar keys(%h_gene) == 0){
	print STDERR	"Parser::parseGTF => Data Structure returns an empty hash\nPossible reasons:\n";
	print STDERR	"\t*Feature level '$levels' is not present in 3rd field of '$infile'\n";
	print STDERR	"\t*chromosome/seqname (chr) or patch chromosome...\n";
	print STDERR	"\t*Filtering tag/Attributes (--filter|-f) option returns no results\nTry --help for help\n";
	exit;
    }
    return \%h_gene;
}


=head2 parsedoubleGTF
    Title   : parsedoubleGTF
    Function: parse a  double .GTF file i.e when 2 gtf line are overlapping for instance
    Example : $refmrna	= Parser:::parsedoubleGTF($infileB, 'exon',  undef, undef, $verbosity);
    Returns : A Feelnc tx-based DSC
    Args:		:
    - infile    : file - a gtf file (mandatory)
    - levels    : string - specify the the level(s) of annotation to be parsed separated by "," (e.g exon,CDS)
    - split     : boolean - 0 or 1 if we want to return a chrom-based dsc
    - refhfilter: hashref - with key value to be extracted e.g 'transcript_id=ENSCAFT00000017943,ENSCAFT00000017911'
    - verbosity	: numeric - level of verbosity
    Note:  only used in addAttribfromIntersect.pl (to be deprecated?)

=cut
sub parsedoubleGTF{

    my ($infile, $levels, $split, $refhfilter, $verbosity) = @_;
    $infile     //= '';
    $levels     //= 'exon';
    $split      //= 0; 	# split data structure by chr
    $refhfilter	//= undef; # filter some interesting (beyond gene_id and transcript_id) for instance 'transcript_biotype' : use transcript_id,gene_id for just 12th fields
                           # in a program with hash (get)option such as $> -f transcript_id=ENSCAFT00000017943,ENSCAFT00000017911
                           # this hashref will be : 'transcript_id' => 'ENSCAFT00000017943,ENSCAFT00000017911'
    $verbosity  //= 1;
    my $baseName  = basename($infile);
    my @all_level = split (",", $levels);

    # Open file
    open GTFFILE, "$infile" or die "Error! Cannot open GTF File ". $infile . ": ".$!;
    my @lines = <GTFFILE>;

    # Store data
    my %h_transcript;

    # counter line
    my $i = 0;

    # label LINE for next if filterFields activated
  LINE:
    foreach my $line  (@lines){

	# verbose
	$i++;
	Utils::showProgress(scalar @lines, $i, "Parse input file: ") if ($verbosity > 0);

	chomp $line;
	next if (($line =~ /^#/) || ($line =~ /^\s*$/) || $line =~ /^track/); # remove weird line

	# split by tab
	my ($chr, $source, $level, $beg_feat, $end_feat, $score, $strand, $frame,  $attributes, $chrB, $sourceB, $levelB, $beg_featB, $end_featB, $scoreB, $strandB, $frameB,  $attributesB) = split(/\t+/, $line);

	# check number of columns
	die "Error:\n[$line] does not look like *double* GTF... not 18 columns\n" if ( !defined $attributes );

	# parse attributes and extract interesting information (filtertag)
	my ($lineremove, $refattrib)   = parseAttributes($attributes, $refhfilter, $verbosity);
	my ($lineremoveB, $refattribB) = parseAttributes($attributesB, $refhfilter, $verbosity);

	next if ($lineremove or $lineremoveB);

	# check for mandatory tags : gene_id and transcript_id
	die "Error:\n[$line] does not contain 'transcript_id' attribute...\n" unless(defined $refattrib->{'transcript_id'});
	die "Error:\n[$line] does not contain 'gene_id' attribute...\n" unless(defined $refattrib->{'gene_id'});
	die "Error:\n[$line] does not contain 2nd 'transcript_id' attribute...\n" unless(defined $refattribB->{'transcript_id'});
	die "Error:\n[$line] does not contain 2nd 'gene_id' attribute...\n" unless(defined $refattribB->{'gene_id'});

	my $transcript = $refattrib->{'transcript_id'};
	my $gene_id    = $refattrib->{'gene_id'};

	my $transcriptB = $refattribB->{'transcript_id'};
	my $gene_idB    = $refattribB->{'gene_id'};

	# delete from ref since we already check they exist and we dont need them anymore
	delete $refattrib->{'transcript_id'};
	delete $refattrib->{'gene_id'};

	delete $refattribB->{'transcript_id'};
	delete $refattribB->{'gene_id'};

	for my $lvl (@all_level) {

	    if ($lvl eq $level){

		my %h_txB;
		$h_txB{$transcript}                      = $transcriptB;
		$h_transcript{$transcript}{$transcriptB} = $gene_idB;
	    }
	}
    }
    print STDERR "\n" if ($verbosity > 0); # because of the showProgress


    # close handle
    close GTFFILE;

    # Test parsing i.e empty hash
    die "Parser::parsedoucleGTF => Data Structure returns an empty hash...Exiting\n" if (scalar keys(%h_transcript) == 0);

    return \%h_transcript;
}


=head2 parseBed
    Title   :	parseBed
    Function:	parse a .BED file
    Example :	$refmrna	= Parser::parseBed($infilegtf);
    Returns : 	A Feelnc tx-based DSC
    Args:
    - infile            : file - a gtf file (mandatory)
    - transcript_biotype: string - since it is not present in .BED it could be added (default 'NA')
    - verbosity         : numeric - level of verbosity
    Note: to be TESTED and for conformity with parseGTF

=cut
sub parseBed{

    my ($infile, $transcript_biotype, $verbosity) = @_;
    $infile             //= "";
    $verbosity          //= 1;
    $transcript_biotype //= "NA";

    my $baseName = basename($infile);

    # Open file
    open BEDFILE, "$infile" or die "Error! Cannot open BED File ". $infile . ": ".$!;
    my @lines = <BEDFILE>;

    # Store data
    my %h_transcript;

    # counter line
    my $i = 0;

    # test if transcript_id ($4) is seen multiple time
    my %seen;

    # bool bed12
    my $bed12 = 0;
    foreach my $line  (@lines){
	chomp $line;

	##########
	# if bed12
	if ($line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+([+-])\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)/) {

	    $bed12 = 1;

	    my $chr        = $1;
	    my $level      = "transcript";
	    my $beg_feat   = $2;
	    my $end_feat   = $3;
	    my $transcript = $4;
	    my $score	   = $5;
	    my $strand     = $6;
	    my $cds_start  = $7;
	    my $cds_end    = $8;
	    my $color      = $9;
	    my $frame      = ".";
	    my $nbExons	   = $10;
	    my $exonSizes  = $11;
	    my $exonStarts = $12;
	    my $extrafield = $13;

	    my @exonSizes  = split(/,/, $exonSizes);
	    my @exonStarts = split(/,/, $exonStarts);

	    # Dedup transcript id
	    $seen{$transcript}++;
	    if ($seen{$transcript}>1) {
		$transcript = $transcript."_dup".$seen{$transcript};
	    }

	    # Transcript levels
	    $h_transcript{$transcript}->{"chr"}	               = $chr;
	    $h_transcript{$transcript}->{"transcript_biotype"} = $transcript_biotype;
	    $h_transcript{$transcript}->{"startt"}             = Utils::min2($h_transcript{$transcript}->{"startt"}, $beg_feat);
	    $h_transcript{$transcript}->{"endt"}               = Utils::max2($h_transcript{$transcript}->{"endt"}, $end_feat);
	    $h_transcript{$transcript}->{"score"}	       = $score;
	    $h_transcript{$transcript}->{"strand"}             = $strand;
	    $h_transcript{$transcript}->{"frame"}	       = $frame;
	    $h_transcript{$transcript}->{"gene_id"}	       = $transcript;


	    # Feature infos are stored in N hashes from a array reference
	    # $h{$tr}->{"feature"} is the reference on a array
	    # @{$h_transcript{$transcript}->{"feature"}} to dereference it
	    # To obtain a particular element in the array : ${ $h{$tr}->{"feature"} }[0]
	    for (my $j = 0; $j < $nbExons; $j++) {

		my %feature = ( "start"     => $h_transcript{$transcript}->{"startt"} + $exonStarts[$j] -1,
				"end"       => $h_transcript{$transcript}->{"startt"} + $exonStarts[$j] + $exonSizes[$j],
				"feat_level" => 'exon',
				"extrafield" => $extrafield);
		push (@{$h_transcript{$transcript}->{"feature"}}, \%feature);
	    }

	    ##########
	    # if bed6
	} elsif ($line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+([+-])/) {


	    my $chr         = $1;
	    my $level	    = "transcript";
	    my $beg_feat    = $2;
	    my $end_feat    = $3;
	    my $transcript  = $4;
	    my $score	    = $5;
	    my $strand      = $6;


	    # Dedup transcript id
	    $seen{$transcript}++;
	    if ($seen{$transcript}>1) {
		$transcript = $transcript."_dup".$seen{$transcript};
	    }

	    # Transcript levels
	    $h_transcript{$transcript}->{"chr"}			= $chr;
	    $h_transcript{$transcript}->{"transcript_biotype"}	= $transcript_biotype;
	    $h_transcript{$transcript}->{"startt"}		= Utils::min2($h_transcript{$transcript}->{"startt"}, $beg_feat);
	    $h_transcript{$transcript}->{"endt"}		= Utils::max2($h_transcript{$transcript}->{"endt"}, $end_feat);
	    $h_transcript{$transcript}->{"score"}		= $score;
	    $h_transcript{$transcript}->{"strand"}              = $strand;
	    $h_transcript{$transcript}->{"gene_id"}		= $transcript;


	} elsif (($line =~ /^#/) || ($line =~ /^\s*$/) || $line =~ /^track/)  {
	    # skip
	} else {
	    die "\nError! Cannot parse this line in '$infile':\n", $line,"\n! Aborting.\nCheck the format of your input file\n";
	}
	# 		verbose
	if ($verbosity > 0){
	    $i++;
	    Utils::showProgress(scalar @lines, $i, "Parse input file: ");
	}
    }
    print STDERR "\n" if ($verbosity > 0); # because of the showProgress

    # close handle
    close BEDFILE;

    # Test parsing i.e empty hash
    die " Error! Invalid parsing step (empty hash)...Check that your file is in .bed format...\n" unless (keys(%h_transcript)>0);

    ################################################################################
    # Sort exons if bed12
    my $refh = \%h_transcript;
    if ($bed12){
	$refh = ExtractFromHash::sortExons(\%h_transcript, $verbosity)
    }
    # return hash
    return $refh;

}


=head2 parseGTFgnlight
	Title		:	parseGTFgnlight
	Function	:	parse a .GTF file into gene-based feelnc dsc with only gene ranges (wihout transcript and exon levels)
	Example		:	$refgene	= Parser::parseGTF($infilegtf);
	Returns		: 	A Feelnc gene-based DSC
	Args		:
					- infile	: file - a gtf file (mandatory)
					- split		: boolean - 0 or 1 if we want to return a chrom-based dsc
					- verbosity	: numeric - level of verbosity

=cut
sub parseGTFgnlight {

    my ($infile, $split,  $verbosity) = @_;
    $infile      //='';
    $split       //= 0;	# split data structure by chr
    $verbosity   //= 1;
    my $basename   = basename($infile);

    # Open file
    open GTFFILE, "$infile" or die "Error! Cannot open GTF File ". $infile . ": ".$!;
    print STDERR "Parsing parseGTFlight file '$basename'...\n";

    # count nb of lines
    my $lc=Utils::countlinefile($infile);

    # Store data
    my %h_gene;

    # counter line
    my $i = 0;


    while (<GTFFILE>){

	# verbose
	$i++;
	Utils::showProgress($lc, $i, "Parse input file: ") if ($verbosity > 0);

	#
	chop;
	next if /^#/;
        next if /^$/;
        next if /^track/;

	# split by tab
	my ($chr, $source, $level, $beg_feat, $end_feat, $score, $strand, $frame, $attributes) = split(/\t/);


	# check number of columns
	if ( !defined $attributes ) {
	    die "Error:\n[$_] does not look like GTF... not 9 columns\n Check that fields are tab-separated.\n";
	}

	# parse attributes and extract interesting information (filtertag)
	my $refattrib = parseAttributes($attributes, undef, $verbosity);
	next unless (defined $refattrib);


	# check for mandatory tags : gene_id and transcript_id
	die "[$_] does not contain 'gene_id' attribute...\n" unless (defined $refattrib->{'gene_id'});
	my $transcript = $refattrib->{'transcript_id'};
	my $gene_id    = $refattrib->{'gene_id'};

	# delete from ref since we already check they exist and we dont need them anymore
	delete $refattrib->{'transcript_id'};
	delete $refattrib->{'gene_id'};

	# gene level
	$h_gene{$gene_id}->{"chr"}    = $chr;
	$h_gene{$gene_id}->{"source"} = $source;
	$h_gene{$gene_id}->{"startg"} = Utils::min2($h_gene{$gene_id}->{"startg"}, $beg_feat);
	$h_gene{$gene_id}->{"endg"}   = Utils::max2($h_gene{$gene_id}->{"endg"}, $end_feat);
	$h_gene{$gene_id}->{"score"}  = $score;
	$h_gene{$gene_id}->{"strand"} = $strand;


    }
    print STDERR "\n" if ($verbosity > 0); # because of the showProgress

    # close handle
    close GTFFILE;

    # Test parsing i.e empty hash
    if (scalar keys(%h_gene) == 0){
	die "Parser::parseGTFlight => Data Structure returns an empty hash!\n";
    }

    # Do we split by chr?
    if ($split){
	my $refhchr = splitbyChr(\%h_gene, $verbosity);

	return $refhchr;

    } else {
	return \%h_gene;
    }

}


=head2 GTF2GTFgnlight
	Title	: GTF2GTFgnlight
	Function: convert a full gene-based feelnc dsc into light gene-based feelnc dsc with only gene ranges (wihout transcript and exon levels)
	Example	: GTF2GTFgnlight ($h, $split, $verbosity);
	Returns	: A Feelnc light gene-based DSC
	Args	:
		  - href     : ref - hashref of full gene-based feelnc dsc
		  - split    : boolean - 0 or 1 if we want to return a chrom-based dsc
		  - verbosity: numeric - level of verbosity
	Note: see FEELnc_codpot.pl

=cut
sub GTF2GTFgnlight{
    my ($href, $split, $verbosity) = @_;
    $split     //= 0; # split data structure by chr
    $verbosity //= 0;

    die "GTF2GTFgnlight: refhash not defined " unless (ref($href));

    print STDERR "Parser::GTF2GTFgnlight with split $split\n" if ($verbosity > 1);

    my %h_gene;
    foreach my $tx (keys %{$href}){

	my $gene_id = $href->{$tx}->{"gene_id"};
	# gene level
	$h_gene{$gene_id}->{"chr"}    = $href->{$tx}->{'chr'};
	$h_gene{$gene_id}->{"source"} = $href->{$tx}->{'source'};
	$h_gene{$gene_id}->{"startg"} = Utils::min2($h_gene{$gene_id}->{"startg"}, $href->{$tx}->{'startt'});
	$h_gene{$gene_id}->{"endg"}   = Utils::max2($h_gene{$gene_id}->{"endg"}  , $href->{$tx}->{'endt'});
	$h_gene{$gene_id}->{"score"}  = $href->{$tx}->{'score'};
	$h_gene{$gene_id}->{"strand"} = $href->{$tx}->{'strand'};

    }
    # Test conversion
    if (scalar keys(%h_gene) == 0){
	die "Parser::GTF2GTFgnlight => Data Structure returns an empty hash...\n";
    }


    # Do we split by chr?
    if ($split){
	my $refhchr = splitbyChr(\%h_gene, $verbosity);

	return $refhchr;

    } else {
	return \%h_gene;
    }

}


1;
