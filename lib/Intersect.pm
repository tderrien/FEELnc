package Intersect;

# 
$VERSION = v0.0.1;

use warnings;
use strict;
use Utils;
use ExtractFromHash;
use ExtractFromFeature;
use Data::Dumper;
use Parallel::ForkManager;

$| = 1;

# Input:
# - chr 	= chr name
# - reftx	= a reference on a FEELnc tx Datastructure
# Output
# - a reference on a filtered FEELnc tx Datastructure by $chr
sub filterChrFromHash{
    my ($chr, $refhtx) = @_;
    
    die "filterChrFromHash: $chr not defined " unless (defined $chr);
    
    # Get the transcript id that match the condition
    my @txok 	= grep { $refhtx->{$_}->{'chr'} eq $chr } keys %{$refhtx};
    
    # Extract/ filter out the refhtx wrt to the array of ok tx
    my %hfilter = map { $_ => $refhtx->{$_} } @txok;

    return \%hfilter;
}

# Input:
# - reftx	= a reference on a FEELnc tx Datastructure hash
# - nbpart 	= the number of sub datastructrue to be returned (i.e corresponding to the number of cpus of rinstance)
# Output
# - an array of 'nbpart' arrays of tx ids
sub splithashref{

    my ($refh, $nbpart, $verbosity) = @_;
    
    my @k = keys %{$refh};
    my $slice ;
    
    if ($nbpart > scalar @k){
	$slice = 1;
    } else {
	$slice = int (scalar @k /$nbpart);
    }
    my @aoa ; # array of arrays of keys
    push @aoa, [ splice @k, 0, $slice ] while @k;
    
    return \@aoa;

}



# from 2 hash ref of split annotation => return an array which contains the list of A chrom that are not in B
sub compare_chr_distrib{
    my ($refha, $refhb) = @_;
    
    my %cmp = map { $_ => 1 } keys %{$refha};
    
    for my $key (keys %{$refhb}) {
        delete $cmp{$key} if exists $cmp{$key};
    }
    # 	my @arkey_justinA = keys %cmp;

    return \%cmp;
}




# Input:
#	- $refh1		: ref on a FEELnc tx Datastructure hash
# 	- $refh2		: ref on a FEELnc tx Datastructure hash
# 	- $stranded		: mode 0 is unstranded // 1 if required same strandness // -1 for oppposite strand
# 	- $fraction		: min proportion of overlap from tx1 to be considered as overlapping
# 	- $monoexonic	: whether to keep mono exonic as: -1 keep monoexonicAS, 1 keep all monoexonic (if 0 their should be no monoexonic at this stage)
# 	- $verbosity	: level of verbosity
# Output:
#	- %matchingtx 	: hash of matching txs such as (tx1 => tx2)
sub getOverlapping{

    my ($refh1, $refh2, $stranded, $fraction, $monoexonic, $linconly, $verbosity)	= @_;

    $stranded 		//= 0; 	# default unstranded (see also monoexonic test)
    $fraction 		//= 0;	# if the fraction of the lncRNA candidate overlap with >0 nt (ie 1bp) a mRNA, it will be removed
    $monoexonic		//= 0;	
    $linconly		//= 0;
    $verbosity 		//= 1;
    
    # hash containing the matching lncRNA <=> mRNA
    my %matchingtx = ();
    my %matchingtxAS;
    
    # nb tx1 and
    my $nb_tx1 = keys( %{$refh1} );
    my $i = 0;
    
    foreach my $tr1 ( ExtractFromHash::sortTxsStartt($refh1) ) {
	
	$i++;
	# verbose on file1
	if ($verbosity > 0){ Utils::showProgress( $nb_tx1, $i, "Intersect fileA: ");}

	# tx1 data
	my $start1  = $refh1->{$tr1}->{"startt"};
	my $end1    = $refh1->{$tr1}->{"endt"};
	my $strand1 = $refh1->{$tr1}->{"strand"};
	
	# we compute nb 
	my $tr1nbexon 	= ExtractFromFeature::features2nbExon($refh1->{$tr1}->{"feature"});

	next if ($tr1nbexon == 1 && $monoexonic == 1); # keep all monoexonic if == 1
	my %monoAS;
	
	# 2nd tx
	foreach my $tr2 (  ExtractFromHash::sortTxsStartt($refh2) ){			

	    # tx1 data
	    my $start2	 	= $refh2->{$tr2}->{"startt"};
	    my $end2	 	= $refh2->{$tr2}->{"endt"};
	    my $strand2 	= $refh2->{$tr2}->{"strand"};

	    # trick to speed  loop
	    if ($monoexonic != -1 && $tr1nbexon != 1){
		next if ($end2   < $start1 );            
		last if ($start2 > $end1   );
	    }		
	    # If linconly we only test overlap at the tx level (not exon)
	    if ($linconly){
		if ( Utils::foverlap( $start1, $end1, $start2, $end2, $strand1, $strand2, $stranded) ){
		    $matchingtx{$tr1} = $tr2;
		    last;				
		}
	    }
	    ##### Overlap window  ########

	    
	    # get cumulsize of lncRNA cDNA (exon level)
	    my $tr1_cdna_size = ExtractFromHash::cumulSize($refh1->{$tr1}->{"feature"});
	    
	    # Get the number  of bp overlap
	    my $cumul_overlap_size = 0;
	    my $fractionoverexon1  = 0;

	    
	    # compute all sense overlapping (whether monoexonic or not)
	    $cumul_overlap_size = ExtractFromFeature::intersectFeatures($refh1->{$tr1}->{'feature'}, $refh2->{$tr2}->{'feature'}, '1', $verbosity);
	    $fractionoverexon1  = $cumul_overlap_size/$tr1_cdna_size;


	    my $cumul_overlap_sizeAS = 0;
	    my $fractionoverexon1AS  = 0;			
	    
	    # compute antisense overlapping for monoexonic
	    if ($monoexonic == -1 && $tr1nbexon == 1){
		$cumul_overlap_sizeAS = ExtractFromFeature::intersectFeatures($refh1->{$tr1}->{'feature'}, $refh2->{$tr2}->{'feature'}, '-1', $verbosity);
		$fractionoverexon1AS  = $cumul_overlap_sizeAS/$tr1_cdna_size;			
	    }
	    
	    
	    # intersectFeatures returns > 0 if there is an ANTISENSE OVERLAP (otherwise it is  0 for sense or if strand is ".")
	    # if no AS match or uncertainty with either lncRNA or mRNA strand we remove the lncRNA

	    # If Sense Overlap
	    if ($cumul_overlap_size > 0 && $fractionoverexon1 > $fraction){
		$matchingtx{$tr1} = $tr2  ;	
		last;
	    } else {
		# If NO AS overlap
		if ($cumul_overlap_sizeAS == 0 ){
		    
		    # if monoexonic, we remove it
		    if ($tr1nbexon == 1 && $monoexonic == -1){ 
			$matchingtx{$tr1} = $tr1."_woAS"  ;
			last if ($start2 > $end1);
		    }
		    # if AS overlap with a mRNA, we exclude it from the remove list :)
		} else {
		    
		    delete $matchingtx{$tr1};
		}	
	    }
	}

    }
    print STDERR "\n" if ($verbosity > 0); # because of the showProgress
    
    return %matchingtx;
}

sub getIntersect{

    my ($refh1, $refh2, $stranded, $fraction, $verbosity) = @_;

    $stranded 		//= 0; 	# default unstranded (see also monoexonic test)
    $fraction 		//= 0;	# if the fraction of the lncRNA candidate overlap with >0 nt (ie 1bp) a mRNA, it will be removed
    $verbosity 		//= 1;
    
    # nb tx1 and
    my $nb_tx1 = keys( %{$refh1} );
    my $i = 0;

    # array that will contain the hashes of matching transcripts
    my @matchingpairs = ();	

    foreach my $tr1 ( ExtractFromHash::sortTxsStartt($refh1) ) {
	
	$i++;
	# verbose on file1
	if ($verbosity > 0){ Utils::showProgress( $nb_tx1, $i, "Intersect fileA: ");}

	# tx1 data
	my $start1  = $refh1->{$tr1}->{"startt"};
	my $end1    = $refh1->{$tr1}->{"endt"};
	my $strand1 = $refh1->{$tr1}->{"strand"};
	

	# 2nd tx
	foreach my $tr2 (  ExtractFromHash::sortTxsStartt($refh2) ){			

	    # tx1 data
	    my $start2  = $refh2->{$tr2}->{"startt"};
	    my $end2    = $refh2->{$tr2}->{"endt"};
	    my $strand2 = $refh2->{$tr2}->{"strand"};

	    # trick to speed  loop
	    last if ($start2 > $end1);
	    next if ($end2   < $start1);            

	    
	    ##### Overlap window  ########
	    
	    # get cumulsize of lncRNA cDNA (exon level)
	    my $tr1_cdna_size = ExtractFromHash::cumulSize($refh1->{$tr1}->{"feature"});
	    
	    # Get the number  of bp overlap
	    my $cumul_overlap_size = ExtractFromFeature::intersectFeatures($refh1->{$tr1}->{'feature'}, $refh2->{$tr2}->{'feature'}, $stranded, $verbosity);
	    my $fractionoverexon1  = $cumul_overlap_size/$tr1_cdna_size;
	    
	    # intersectFeatures returns > 0 if there is an OVERLAP with respect to strand
	    # 	* if s = 0 [default] and A and B on SAME or DIFFERENT strand : return size of overlap or 0 if not overlap
	    # 	* if s = 1           and A and B on SAME              strand : return size of overlap or 0 if not overlap or of overlap but different strand
	    if ($cumul_overlap_size > 0 && $fractionoverexon1 > $fraction){
		my %matchingtx;
		$matchingtx{$tr1}	=	$tr2;
		push @matchingpairs, \%matchingtx;

		# 				last if $uniquexa;
	    }
	}

    }
    print STDERR "\n" if ($verbosity > 0); # because of the showProgress
    
    return \@matchingpairs;
}

# if we consider all features from transcripts
sub Intersect2Hsplit{


    my ($refh1, $refh2, $strand, $fraction, $wb, $mode, $reciprocal, $verbosity)	= @_;
    my %h1      = %{$refh1}; # parsing a hash in sub need dereference in shift
    my %h2      = %{$refh2}; # parsing a hash in sub need dereference in shift
    $strand   //= 0; #require  same strandness for overlap
    $fraction //= 0; # min fraction overlap 
    $wb       //= 0; #Write the original entry in B for each overlap : WARN if mode= transcript : and 2 exons of txB overlap 1ex of txA => print only the first exon of B (see printDoubleGtfHash_Tx)
    $mode     //= "exon"; # could be exon: report all exons of A that are overlap by exon of B whenever some exons of A are not overlaped
                          # or txexon    : report all exons of A that are overlap by exon of B if all exons of A are overlapped  
                          # or txseq     : report Transcript of A that are overlapped by Tx B with fraction of $fraction  (independently of a match number of exons)

    $reciprocal //= 0; # Require that the fraction overlap be reciprocal for A and B. # restricted by $mode
    $verbosity  //= 1;
    
    my $nb_tx1 = keys( %h1 );
    my $i = 0;

    my $nb_tx2 = keys( %h2 );
    my $j = 0;	
    
    
    # Sort both hash by chr and transcript start ie startt
    foreach my $tr1 (
	sort {
	    $h1{$a}->{'chr'} cmp $h1{$b}->{'chr'} or
		$h1{$a}->{'startt'} <=> $h1{$b}->{'startt'} }
	keys %h1
	) {
       
        # verbose on file1
        $i++;
	if ($verbosity > 0){ Utils::showProgress( $nb_tx1, $i, "Analyse fileA: ");}


	# Sort both hash by chr and transcript start ie startt
        foreach my $tr2 (
	    sort { 
		$h2{$a}->{'chr'} cmp $h2{$b}->{'chr'} or
                    $h2{$a}->{'startt'} <=> $h2{$b}->{'startt'} }
	    keys %h2
	    ) {

	    if ($tr1 eq $tr2 ){print STDERR "SAME\t"}
	    # since hashes are ordered by chr => delete tx2 if its chr is lower than chr1
            if ($h2{$tr2}->{"chr"}  lt   $h1{$tr1}->{"chr"}){
                delete $h2{$tr2};
                next;
            }
            
	    # since hashes are ordered by chr => last if chr2 > chr1
            last if ($h2{$tr2}->{"chr"}  gt   $h1{$tr1}->{"chr"} );

	    # end 2nd loop hash if tx2 > tx1
            if ($h2{$tr2}->{"chr"}  eq   $h1{$tr1}->{"chr"}){
		
		last if ($h2{$tr2}->{"startt"} > $h1{$tr1}->{"endt"});
		next if ($h2{$tr2}->{"endt"} < $h1{$tr1}->{"startt"});
	    }            
	    
	    # skip if strand required is defined and different strand
            next if ( $strand && ($h1{$tr1}->{"strand"}  ne   $h2{$tr2}->{"strand"} ));
	    
	    # Go at the feature levels
	    #######################
	    
	    my $sizeover =0;  # intialize overlap
	    
	    my $feat2_over_feat1 = 0; # nb of feature/exon that will be overlapped for tx1
	    my $feat1_over_feat2 = 0; # nb of feature/exon that will be overlapped for tx2
    	    
	    my $cumul_overlap_size = 0;
	    
	    $cumul_overlap_size	= ExtractFromFeature::intersectExons($h1{$tr1}->{'feature'}, $h2{$tr2}->{'feature'});

            
    	    my $tr1_size	= ExtractFromHash::cumulSize($h1{$tr1}->{"feature"});
    	    my $tr2_size	= ExtractFromHash::cumulSize($h2{$tr2}->{"feature"});
	    
	    my $fractionover1	= $cumul_overlap_size/$tr1_size;
	    my $fractionover2	= $cumul_overlap_size/$tr2_size;

	    if ($mode eq "txseq"){			
		
		print $tr1,"\t",$h1{$tr1}->{"gene_id"},"\t", $h1{$tr1}->{"chr"},":",$h1{$tr1}->{"startt"},"-",$h1{$tr1}->{"endt"};
		print "\t",$tr2,"\t",$h2{$tr2}->{"gene_id"},"\t", $h2{$tr2}->{"chr"},":",$h2{$tr2}->{"startt"},"-",$h2{$tr2}->{"endt"},"\t";
		printf ("\t%.2f\t%.2f\n",$fractionover1, $fractionover2);
	    }
        }
    }
    print STDERR "\n" if ($verbosity > 0); # because of the showProgress
}


# Fork cleaning
sub Intersect2HsplitFork_clean{

    my ($refh1, $refh2,  $keepduptx, $keepdupint, $fh, $verbosity) = @_;
    $fh	       //= undef;
    $verbosity //= 1;

    # variable that should be options
    my $stranded 		= 1; #require  same strandness for overlap	
    my $fraction 		= 1;
    
    # nb tx1 and
    my $nb_tx1 = keys( %{$refh1} );
    my $i = 0;

    my $nb_tx2 = keys( %{$refh2} );
    my $j = 0;   
    
    my %seentx=();
    
    # Sort both hash by chr and transcript start ie startt
    foreach my $tr1 (sort { $refh1->{$a}->{'startt'} <=> $refh1->{$b}->{'startt'} } keys %{$refh1}) {
	

	############
	# skip tx if removed from 1st loop in 2nd hash;
	###########
	if ($seentx{$tr1}){
	    delete $refh1->{$tr1};
	    next;
	}
	
	# infos
	$i++;
	
	# verbose on file1
	if ($verbosity > 0){ Utils::showProgress( $nb_tx1, $i, "Analyse fileA: ");}


	# Sort both hash by chr and transcript start ie startt
	foreach my $tr2 (sort { $refh2->{$a}->{'startt'} <=> $refh2->{$b}->{'startt'} } keys %{$refh2}){			
	    
	    # trick to speed (?) loop
	    last if ($refh2->{$tr2}->{"startt"} > $refh1->{$tr1}->{"endt"});
	    next if ($refh2->{$tr2}->{"endt"} < $refh1->{$tr1}->{"startt"});            

	    ##### Overlap window  ########
	    
	    # ############## Exon overlap
	    if (!$keepduptx){
		my $tr1_cdna_size = ExtractFromHash::cumulSize($refh1->{$tr1}->{"feature"});
		my $tr2_exon_size = ExtractFromHash::cumulSize($refh2->{$tr2}->{"feature"});
		
		# Exon overlap
		my $cumul_overlap_size = ExtractFromFeature::intersectFeatures($refh1->{$tr1}->{'feature'}, $refh2->{$tr2}->{'feature'}, $stranded, $verbosity);
		
		
		# proportion
		my $fractionoverexon1	= $cumul_overlap_size/$tr1_cdna_size;
		my $fractionoverexon2	= $cumul_overlap_size/$tr2_exon_size;
		
		# remove same tx exons in 2
		if ( ($tr1 ne $tr2) && ($fractionoverexon1 == 1) && ($fractionoverexon2 ==1) ){

		    if ($fh){
			print $fh "SameExons: $tr1 - $tr2... removing $tr2\n";
		    }else {
			print STDERR "SameExons: $tr1 - $tr2... removing $tr2\n";
		    }
		    $seentx{$tr2}++;
		}
	    }
	    
	    # ############## Intron overlap			
	    if (!$keepdupint && scalar(@{$refh1->{$tr1}->{'feature'}})>1 && scalar(@{$refh2->{$tr2}->{'feature'}})>1){
		my $fractionoverintron1 = 0;
		my $fractionoverintron2 = 0;

		my $refintrons1	= ExtractFromFeature::getIntrons($refh1->{$tr1}->{'feature'});
		my $refintrons2	= ExtractFromFeature::getIntrons($refh2->{$tr2}->{'feature'});	
		
		my $tr1_intron_size = ExtractFromHash::cumulSize($refintrons1);
		my $tr2_intron_size = ExtractFromHash::cumulSize($refintrons2);
		
		# intron overlap
		my $cumul_intron_size = ExtractFromFeature::intersectFeatures($refintrons1, $refintrons2, $stranded, $verbosity);
		
		#  intron overlap
		$fractionoverintron1 = $cumul_intron_size/$tr1_intron_size;
		$fractionoverintron2 = $cumul_intron_size/$tr2_intron_size;
		
		# remove same tx introns in 2
		if ($tr1 ne $tr2 && $fractionoverintron1 == 1 && $fractionoverintron2 ==1){
		    
		    if ($fh){
			print $fh "SameIntrons: $tr1 - $tr2... removing $tr2\n";
		    }else {
			print STDERR "SameIntrons: $tr1 - $tr2... removing $tr2\n" ;
		    }					
		    $seentx{$tr2}++;
		}
	    }	
        }
    }
    print STDERR "\n" if ($verbosity > 0); # because of the showProgress
    
    return %{$refh1};
} 


# if we compare all features from transcripts
sub Intersect2HsplitFork_compare{

    my ($refh1, $refh2, $stranded, $uniquexa, $uniquexb, $common, $fraction, $verbosity)	= @_;

    $fraction 		//= 1;
    $verbosity 		//= 1;

    # nb tx1 and
    my $nb_tx1 = keys( %{$refh1} );
    my $i = 0;

    my $nb_tx2 = keys( %{$refh2} );
    my $j = 0;	    
    
    # array that will contain the hashes of matching transcripts
    my @matchingpairs = ();
    
    my $count=0;
    
    # Sort both hash by chr and transcript start ie startt
    foreach my $tr1 (sort { $refh1->{$a}->{'startt'} <=> $refh1->{$b}->{'startt'} } keys %{$refh1}) {

	# infos
	$i++;
	
	# verbose on file1
	if ($verbosity > 0){ Utils::showProgress( $nb_tx1, $i, "Analyse fileA: ");}


	# Sort both hash by chr and transcript start ie startt
	foreach my $tr2 (sort { $refh2->{$a}->{'startt'} <=> $refh2->{$b}->{'startt'} } keys %{$refh2}){			
	    # trick to speed (?) loop
	    last if ($refh2->{$tr2}->{"startt"} > $refh1->{$tr1}->{"endt"});
	    next if ($refh2->{$tr2}->{"endt"} < $refh1->{$tr1}->{"startt"});            

	    ##### Overlap window  ########
	    
	    # ############## Exon overlap
	    my $tr1_cdna_size		= ExtractFromHash::cumulSize($refh1->{$tr1}->{"feature"});
	    my $tr2_exon_size		= ExtractFromHash::cumulSize($refh2->{$tr2}->{"feature"});
	    
	    # Exon overlap
	    my $cumul_overlap_size	= ExtractFromFeature::intersectFeatures($refh1->{$tr1}->{'feature'}, $refh2->{$tr2}->{'feature'}, $stranded, $verbosity);
	    
	    # proportion
	    my $fractionoverexon1	= $cumul_overlap_size/$tr1_cdna_size;
	    my $fractionoverexon2	= $cumul_overlap_size/$tr2_exon_size;

	    print STDERR "($tr1 ne $tr2) && $fractionoverexon1 >= $fraction  && $fractionoverexon2 >= $fraction\n" if($verbosity > 1);
	    
	    #  same tx exons in 2
	    if ( $fractionoverexon1 >= $fraction  && $fractionoverexon2 >= $fraction ){
		my %matchingtx = ();
		$matchingtx{$tr1} = $tr2;
		push @matchingpairs, \%matchingtx;
		
		print STDERR "SameExonsStructure: $tr1 - $tr2...\n" if ($verbosity > 1);
		
	    }
        }

    }
    print STDERR "\n" if ($verbosity > 0); # because of the showProgress
    
    ################	
    my %hfinal;
    # return 
    if ($common){
	return @matchingpairs;
    } else {
	if ($uniquexa){
	    
	    # get all keys i.e matching transcripts from A from array of hashes and put in $matchfromA;
	    my %matchfromA;
	    for my $hashref (@matchingpairs) {
		for my $key (keys %$hashref) {
		    $matchfromA{$key}++;
		}
	    }
	    # Exclude matching A txs
	    for my $tx1 (keys %{$refh1}){
		$hfinal{$tx1} = $refh1->{$tx1}  if (!exists($matchfromA{$tx1})),
	    }
	}
	if ($uniquexb){
	    # get all values i.e matching transcripts from B from array of hashes and put in $matchfromB;
	    my %matchfromB;
	    for my $hashref (@matchingpairs) {
		for my $key (keys %$hashref) {
		    $matchfromB{$hashref->{$key}}++;
		}
	    }
	    # Exclude matching B txs
	    for my $tx2 (keys %{$refh2}){
		$hfinal{$tx2} = $refh2->{$tx2} if (!exists($matchfromB{$tx2}));
	    }
	}
    }
    return %hfinal;
}


1;
