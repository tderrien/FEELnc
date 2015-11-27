package ExtractFromFeature;

# AIM :
# Extract data from Feature data_structure
# 'ENSCAFT00000018069' => {
#                                     'gene_id' => 'ENSCAFG00000011362',
#                                     'chr' => 'X',
#                                     'frame' => '.',
#                                     'score' => '.',
#                                     'endt' => 6433061,
#                                     'source' => 'protein_coding',
#                                     'feature' => [
#                                                    {
#                                                      'transcript_name' => 'TBL1Y-201',
#                                                      'end' => '6300870',
#                                                      'exon_number' => '1',
#                                                      'gene_biotype' => 'protein_coding',
#                                                      'exon_id' => 'ENSCAFE00000291557',
#                                                      'feat_level' => 'exon',
#                                                      'gene_name' => 'TBL1Y',
#                                                      'start' => 6300797
#                                                    },
#                                                    {
#                                                      'transcript_name' => 'TBL1Y-201',
#                                                      'end' => '6369718',
#                                                      'exon_number' => '2',
#                                                      'gene_biotype' => 'protein_coding',
#                                                      'exon_id' => 'ENSCAFE00000306730',
#                                                      'feat_level' => 'exon',
#                                                      'gene_name' => 'TBL1Y',
#                                                      'start' => 6369630
#                                                    },



$VERSION = v0.00001;

use warnings;
use strict;
# use List::MoreUtils qw( minmax );

use Data::Dumper;
use Utils;
use StringUtils;
use Bio::DB::Fasta;


$| = 1;

# Check whether all features from a tx contain at least one CDS feat_level
sub checkCDS {
    my ($refarray) = @_;

    die "ExtractFromFeature::checkCDS => not features for refarray...\n" if (scalar @{$refarray}==0);

    my $hasCDS = 0;

    foreach my $exon (@{$refarray}) {
	if ($exon->{'feat_level'} =~ /CDS/i){
	    return 1;
	}
    }
    return $hasCDS;
}

# Sort features according to strand
sub sortFeatures{
    my ($refarray, $strand) = @_;

    die "ExtractFromFeature::feature2CDS => transcript  not defined...\n" if (not defined $refarray);

    my @features;
    # reverse order if negative strand
    if ($strand eq "-" || $strand eq "-1"){

	@features = sort { $b->{"start"} <=> $a->{"start"}} @{$refarray};

    } elsif ($strand eq "+" || $strand eq "1" || $strand eq "."){

	@features = sort { $a->{"start"} <=> $b->{"start"}} @{$refarray};

    }else {
	die "ExtractFromFeature::sortFeatures => strand '$strand' unknown...";
    }

    return  @features;

}

# get the cumulative size of features from a array reference
sub features2size{
    my $refonarray   = shift;
    my $verbosity    = shift;
    $verbosity     //= 1;

    # print if verbosity
    print STDERR "Getting cumulative size of exon levels...\n" if ($verbosity > 1);

    # cumulsize
    my $cumulsize = 0;

    # Parse gtfHash to be printed
    foreach my $feat1 (@{$refonarray}) {

	$cumulsize += ($feat1->{"end"} - $feat1->{"start"} +1) if ($feat1->{'feat_level'} eq "exon");

    }
    return $cumulsize;
}

sub features2nbExon{
    my $refonarray   = shift;
    my $verbosity    = shift;
    $verbosity     //= 1;

    # print if verbosity
    print STDERR "Getting number of exon(s)...\n" if ($verbosity > 1);

    # cumulsize
    my $nbex	=	0;

    # Parse gtfHash to be printed
    foreach my $feat1 (@{$refonarray}) {

	$nbex++ if ($feat1->{'feat_level'} eq "exon");

    }

    return $nbex;

}

# return true if dubious biexonic tx (dubious means a tx with a small exon
sub features2dubExon{
    my ($refonarray, $minsize, $verbosity) = @_;
    $minsize   //= 25;
    $verbosity //= 1;

    # print if verbosity
    print STDERR "Getting dubious exonic transcript...\n" if ($verbosity > 1);

    # cumulsize
    my $dubious = 0;

    # Parse gtfHash to be printed
    foreach my $feat1 (@{$refonarray}) {

        my $size = ($feat1->{"end"} - $feat1->{"start"} +1) if ($feat1->{'feat_level'} eq "exon");
        if ($size < $minsize){
            return 1;
            last;
        }

    }

    return $dubious;

}

# Add CDS via ORF object to a reference on array of exons
# Frame is not implemented yet
sub Tx2CDS{
    my ($tx, $refORF, $verbosity) = @_;
    $verbosity //= 1;


    # Test parsing i.e empty hash
    die "ExtractFromFeature::feature2CDS => transcript  not defined...\n" if (not defined $tx);
    die "ExtractFromFeature::feature2CDS => ref ORF not defined...\n" if (not defined $refORF);

    # ORF inf
    my $orf_length = $refORF->{'orflength'} ;  # for ORF/sequence, we are at 0-based coordinates

    # Tx info
    my $strand = $tx->{"strand"};

    my $frame = ".";

    my @newfeatures = @{$tx->{"feature"}};

    # exon number
    my $nbexons = scalar (@{$tx->{"feature"}});

    my $exonfromCDS=0;
    # Iterate over ****SORTED**** wrt strand exons
    foreach my $exon (ExtractFromFeature::sortFeatures($tx->{"feature"}, $strand)) {

	# Current Exon length
	my $exon_length = $exon->{"end"} - $exon->{"start"} +1;

	# CDS start and end in genomic context
	my ($cdss, $cdse);

	# *** Most important trick here ***
	if ($strand eq "+" || $strand eq "." ){
	    $cdss = $exon->{"start"} + $refORF->{'start'};
	    $cdse = $cdss            + $orf_length;
	} else {
	    $cdse = $exon->{"end"}   - $refORF->{'start'};
	    $cdss = $cdse            - $orf_length;
	}

	# Compute the interval matched between exon and projected ORF
	my ($startmatch, $endmatch) = Utils::foverlapmin($exon->{"start"}, $exon->{"end"}, $cdss, $cdse, '+', '+');
	# Size of the matched region
	my $matchlength = $endmatch - $startmatch +1;

	# store features leve i.e CDS 5'UTR and 3'UTR at the moment
	my (%FpUTR, %TpUTR, %CDS, %UTR, %start_codon, %stop_codon);


        # DEBUG
        # main variables are :
        #   * $startmatc
        #   * $refORF->{'start'}
        #   * $matchlength
        #   * $exon_length
        #	* $orf_length : what is the rest of ORF length after deducing
	#  		print STDERR "$nbexons : CDSEx $exonfromCDS: refORF START: ",$refORF->{'start'},"  ORFLEngth: $orf_length - ($startmatch, $endmatch) $matchlength : (",$exon->{'start'},", ",$exon->{'end'},") $exon_length :  $cdss, $cdse\n";
	my $startexon = $exon->{'start'};
	my $endexon   = $exon->{'end'};

	# if match between ORFgenomic and exon => we are in CDS
	if ($startmatch){

            $exonfromCDS++;

            # PLUS STRAND
            if ($strand  eq "+" || $strand eq "."){

		my $FpUTR_beg = $startexon;
		my $FpUTR_end = $startmatch -1;
    	    	my $CDS_beg   = $startmatch;
		my $CDS_end   = $endmatch;
		my $TpUTR_beg = $endmatch;
		my $TpUTR_end = $endexon;

		# if first exon
		if ($exonfromCDS == 1 ) {

                    %FpUTR       = ( "start" => $FpUTR_beg   , "end" => $FpUTR_end       , "feat_level" => 'UTR', 'frame' => $frame);
                    %start_codon = ( "start" => $startmatch  , "end" => $startmatch +2   ,  "feat_level" => 'start_codon', 'frame' => $frame);


		    if ( $refORF->{'check_start'} ) {
			if ( $FpUTR_beg < $FpUTR_end ){
                            push (@newfeatures, \%FpUTR, \%start_codon);

			} else {
                            push (@newfeatures,  \%start_codon);

			}
		    }
		}

		# last exon belonging to CDS
		if ( $orf_length <= $matchlength ){

		    # weirdly, somteimes threre is is a shift between orf length and matchlength on thel ast exons
		    my $diff_RestORFlength = $matchlength - $orf_length;
		    $CDS_end               = $CDS_end - $diff_RestORFlength;


		    #    		        $CDS_beg
		    %TpUTR      = ( "start" => $TpUTR_beg, "end" => $TpUTR_end, "feat_level" => 'UTR', 'frame' => $frame);
		    %stop_codon = ( "start"  => $CDS_end-2, "end" => $CDS_end,  "feat_level" => 'stop_codon', 'frame' => $frame);

		    if ( $refORF->{'check_stop'} ) {

			push (@newfeatures, \%stop_codon);

			# adapt CDS coordinates with stop for + strand
                        $CDS_end = $CDS_end - 3 ;
		    }

            	    # add UTR
                    if ( $TpUTR_beg <  $TpUTR_end ){
                        push (@newfeatures, \%TpUTR)
                    }
		}

		# CDS
                %CDS = ( "start" => $CDS_beg, "end" => $CDS_end,  "feat_level" => 'CDS', 'frame' => $frame);

                push (@newfeatures,  \%CDS);

		##############
		# MINUS STRAND
		##############
            } elsif ($strand  eq "-"){

		my $FpUTR_beg = $endmatch;
		my $FpUTR_end = $startmatch;
    	    	my $CDS_beg   = $startmatch;
		my $CDS_end   = $endmatch;
		my $TpUTR_beg = $startexon;
		my $TpUTR_end = $startmatch;

		# if first exon
		if ($exonfromCDS == 1 ) {

                    %FpUTR       = ( "start" => $FpUTR_beg, "end" => $FpUTR_end, "feat_level" => 'UTR', 'frame' => $frame);
                    %start_codon = ( "start" => $endmatch -2, "end" => $endmatch,  "feat_level" => 'start_codon', 'frame' => $frame);


		    if ( $refORF->{'check_start'} ) {
			if ( $FpUTR_beg < $FpUTR_end ){
                            push (@newfeatures, \%FpUTR, \%start_codon);

			} else {
                            push (@newfeatures,  \%start_codon);

			}
		    }
		}

		# last exon belonging to CDS
		if ( $orf_length <= $matchlength ){

		    # weirdly, somteimes threre is is a shift between orf length and matchlength on thel ast exons
		    my $diff_RestORFlength  = $matchlength - $orf_length;
		    $CDS_beg                = $CDS_beg + $diff_RestORFlength;


		    #    		        $CDS_beg
		    %TpUTR      = ( "start" => $TpUTR_beg, "end" => $TpUTR_end, "feat_level" => 'UTR', 'frame' => $frame);
		    %stop_codon = ( "start"  => $CDS_beg, "end" => $CDS_beg + 2,  "feat_level" => 'stop_codon', 'frame' => $frame);

		    if ( $refORF->{'check_stop'} ) {

			push (@newfeatures, \%stop_codon);

			# adapt CDS coordinates with stop for + strand
			$CDS_beg = $CDS_beg +3;
		    }

            	    # add UTR
                    if ( $TpUTR_beg <  $TpUTR_end ){
                        push (@newfeatures, \%TpUTR)
                    }
		}

		# CDS
                %CDS = ( "start" => $CDS_beg, "end" => $CDS_end,  "feat_level" => 'CDS', 'frame' => $frame);

                push (@newfeatures,  \%CDS);

            }

	    #            print Dumper \@newfeatures;

	    ##################################################
	    # Reinitialiez ORF start and length	for next exon
	    $refORF->{'start'} = 0;            # we are in a match with ORF,
	    $orf_length        = $orf_length - $matchlength;


	    # UTR only
	} else {
	    %FpUTR = ( "start" => $exon->{'start'}, "end" => $exon->{'end'}  ,  "feat_level" => 'UTR', 'frame' => $frame);
	    push (@newfeatures, \%FpUTR);

	    # Reinitialiez UTR size or ORF start
	    $refORF->{'start'}	= $refORF->{'start'} - $exon_length;

	}
    }

    # assign new feature instead of features
    @{$tx->{'feature'}} = @newfeatures;


    return $tx;
}





# Get fasta sequence from a set of features coordinates
sub feature2seq{
    my ($refarray,  $genome, $chr, $strand, $filterCDS, $verbosity)	= @_;

    $chr       //= undef; # since,we are at the feature level, we don't have the chromosome information (need to be defined)
    $filterCDS //= 0; # get only CDS and stop codon sequence
    $verbosity //= 1;

    # initialize sequence string
    my $seqstring = "";

    # Test parsing i.e empty reafarray
    die "ExtractFromFeature::feature2seq => ref not defined...\n" if (not defined $refarray);

    # check $chr is defined
    die "ExtractFromFeature::feature2seq => chromosome '$chr' not defined...\n" if (not defined $chr);


    # create genome DB if does not exist
    my $db = Bio::DB::Fasta->new($genome);


    foreach my $exon (@{$refarray}) {

	## VW modification
	# next if don't want CDS ***
	next if ((!$filterCDS) && ($exon->{'feat_level'} =~ /(CDS|stop_codon|start_codon)/i));

	# we want only CDS
	# if we want CDS but the feature level we look at is not a CDS or the codon stop
	next if (($filterCDS)  && ($exon->{'feat_level'} !~ /(CDS|stop_codon)/i));

	# get coordinates
	my $s       = $exon->{"start"};
	my $e       = $exon->{"end"};
	my $featseq =  ""; # initialize in case pb with bio::db index
	$featseq    =  $db->seq($chr, $s => $e);
	if (!defined $featseq || $featseq eq ""){
	    warn "ExtractFromFeature::feature2seq: your seq $chr:$s-$e returns an empty string!...Check 'chr' prefix between your annotation and genome files or remove your genome index file ('".$genome.".index')...\n";
	    return undef;
	}
	$seqstring .= $featseq;
    }

    #RevComp if strand -
    if ($strand eq '-'|| $strand eq '-1') {
        $seqstring = StringUtils::getRevComp($seqstring, $verbosity);
    }
    return $seqstring;
}


# Get min and max positions from a set of exons coordinates
sub getMinMax{
    my ($refa,  $verbosity) = @_;
    $verbosity //= 1;

    # Test parsing i.e empty hash
    die " Error! ExtractFromFeature: getMinMax return empty hash...\n" if (scalar(@{$refa}) == 0);

    my @array =();

    # put all dat in @array
    foreach my $exon (@{$refa}) {
	push (@array, $exon->{"start"}, $exon->{"end"});
    }
    my @arraysorted = sort { $a <=> $b } @array;

    my $mincoord = $arraysorted[0];
    my $maxcoord = $arraysorted[-1];

    return ($mincoord, $maxcoord);
}

# Get introns coordinates from a array of hashes of exons ($h{$tr}->{'features'}
# return a array of hash(es) of introns
sub getIntrons{
    my ($refarray,  $verbosity)	= @_;
    $verbosity //= 1;

    # Test parsing i.e empty hash
    die "ExtractFromFeature::getIntrons => ref not defined...\n" if (not defined $refarray);

    my @introns =();

    # if transcript mono exonic
    if (scalar (@{$refarray}) == 1){
    	return \@introns;
    } else {
	for (my $i=1; $i<scalar(@{$refarray}); $i++) {
	    my $intron_start	=	${$refarray}[$i-1]->{"end"} + 1;
	    my $intron_end		=	${$refarray}[$i]->{"start"} - 1;
	    my $intron_size		= 	$intron_end - $intron_start + 1;
	    my %feature_intron = ( "start" => $intron_start, "end" => $intron_end, "feat_level" => 'intron', "intron_size" => $intron_size);
	    push (@introns, \%feature_intron);
	}
	return \@introns;
    }
}

# get gene or transcript transcript_biotype from feature
sub getKeyFromFeature{
    my ($refh, $tx, $keyfeat, $verbosity) = @_;
    my %h        = %{$refh};
    $tx        //= ''; # $tx id
    $keyfeat   //= ''; # key feature could be (depending on parsing) 'transcript_name' 'end' 'exon_number' 'gene_biotype' 'exon_id' 'feat_level' 'gene_name' 'start'...
    $verbosity //= 1;

    # Test parsing i.e empty hash
    die " Error! ExtractFromFeature: getKeyFromFeature return empty hash...\n" unless (scalar(keys(%h))>0);

    # die if transcript is not in hash
    die "Error! ExtractFromFeature: $tx not in hash..." unless exists ($h{$tx});

    # return keyfeat or 0
    if ( exists(${$h{$tx}->{"feature"}}[0]->{$keyfeat}) ){
	return ${$h{$tx}->{"feature"}}[0]->{$keyfeat};
    } else {
	return "NA";
    }
}




sub addKeyValAttrib{
    my ($tx, $refhmatchTxGn, $keylabel, $vallabel, $verbosity)	= @_;
    $tx            //= undef; # $tx id
    $refhmatchTxGn //= undef; # ref on a hash that contain matching transcript as key and their geneid as values (from Parser::parsedoubleGtf) { 'ENSCAFT00000041202' => 'ENSCAFG00000026919' },
    $keylabel      //= 'keylabel';
    $vallabel      //= "vallabel";
    $verbosity     //= 1;

    my $txids = 'NA';
    my $gnids = 'NA';

    # 	print "($tx, $refhmatchTxGn, $keylabel, $vallabel, $verbosity)\n";
    # Test parsing i.e empty hash
    die "ExtractFromFeature::addKeyValAttrib => transcript  not defined...\n" if (not defined $tx);

    if (defined $refhmatchTxGn && scalar (keys %{$refhmatchTxGn})){
	my @txids = keys   (%{$refhmatchTxGn});
	my @gnids = values (%{$refhmatchTxGn});
	# Only get Unique gene
	@gnids = Utils::uniq(@gnids);

	$txids = join ",", @txids;
	$gnids = join ",", @gnids;
    }

    my @newfeatures;
    foreach my $feat (@{$tx->{"feature"}}){

	# create new hash to be added
	my %addinfos = ( $keylabel => $txids, $vallabel => $gnids);
	my %merge    = (%addinfos, %{$feat});

	# add to feature
	push (@newfeatures, \%merge);
    }
    # replace old feat by new feat
    @{$tx->{'feature'}} = @newfeatures;

    return $tx;
}



# Intersect features (i.e exons or introns) from a 2 references of an array of hashes
sub intersectFeatures{

    my ($refA, $refB, $stranded, $verbosity)	= @_;
    $refA      //= undef; # $ref on exon array A
    $refB      //= undef; # $ref on exon array B
    $stranded  //= 0;     # mode 0 is unstranded, 1 if required same strande and -1 different strand
    $verbosity //= 1;

    die "ExtractFromFeature::intersectFeatures => refA not defined...\n" unless (defined $refA);
    die "ExtractFromFeature::intersectFeatures => refB not defined...\n" unless (defined $refB);

    my $cumul_overlap_size =0;

    # Loop on exon of A
    foreach my $featureA (@{$refA}) {

	# and B
	foreach my $featureB (@{$refB}) {

	    #  END exons of B if if current exons of B start is greater than end exon A
	    last if ( $featureB->{"start"} > $featureA->{"end"});

	    # compute the size of the overlap with sub from Utils;
	    my $sizeover 	= Utils::foverlapsizerange($featureA->{"start"}, $featureA->{"end"}, $featureB->{"start"}, $featureB->{"end"}, $featureA->{"strand"}, $featureB->{"strand"}, $stranded);

	    # for mode wholeseq
	    $cumul_overlap_size += $sizeover;
	}
    }
    return $cumul_overlap_size;
}

# Intersect features (i.e exons or introns) from a 2 references of an array of hashes
sub intersectFeaturesBool{

    my ($refA, $refB, $stranded, $verbosity)	= @_;

    $refA      //= undef; # $ref on exon array A
    $refB      //= undef; # $ref on exon array B
    $stranded  //= 0;     # mode 0 is unstranded, 1 if required same strande and -1 different strand
    $verbosity //= 1;

    die "ExtractFromFeature::intersectFeatures => refA not defined...\n" unless (defined $refA);
    die "ExtractFromFeature::intersectFeatures => refB not defined...\n" unless (defined $refB);

    # Loop on exon of A
    foreach my $featureA (@{$refA}) {

	# and B
	foreach my $featureB (@{$refB}) {

	    #  END exons of B if if current exons of B start is greater than end exon A
	    last if ( $featureB->{"start"} > $featureA->{"end"});

	    # compute the size of the overlap with sub from Utils;
	    return 1 if ( Utils::foverlap($featureA->{"start"}, $featureA->{"end"}, $featureB->{"start"}, $featureB->{"end"}, $featureA->{"strand"}, $featureB->{"strand"}, $stranded) );

	}
    }
    return undef;
}


1;
