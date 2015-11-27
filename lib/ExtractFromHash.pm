package ExtractFromHash;

$VERSION = v0.0.1;

use warnings;
use strict;
use List::Util 'shuffle';
use Utils;
use StringUtils;
use Data::Dumper;
use Bio::DB::Fasta;
use Data::Dumper;
use Bio::SeqIO;


$| = 1;

sub hash2fasta{

    # parsing a hash in sub need dereference in shift
    my ($h, $genome, $outfasta,  $verbosity) = @_;

    die "Error:hash2fasta: Cannot read the genome file/dir...\n" if (! -r $genome && !-d $genome);
    $outfasta  //= 'hash2fasta';
    $verbosity //= 1;

    # create genome DB
    my $db = Bio::DB::Fasta->new($genome);

    my $seqOUT;
    # Output
    if (defined $outfasta){
	$seqOUT	= Bio::SeqIO ->new(-format => 'fasta', -file => ">$outfasta", -alphabet =>'dna');
    }

    for my $tr (keys %{$h}){
	
	#Initalize sequence:
	my $seqstring	= "";
	my $id_sequence	= "";
	my $cpt         = 0;
	
	foreach my $exon (@{$h->{$tr}->{"feature"}}) {
	    
	    $cpt++;
	    my $chr     = $h->{$tr}->{"chr"};
	    my $s       = $exon->{"start"};
	    my $e       = $exon->{"end"};
	    $seqstring .=  $db->seq($chr, $s => $e); # else no slop
	}
	
	#RevComp if strand -
	if ( $h->{$tr}->{"strand"} eq '-' || $h->{$tr}->{"strand"} eq '-1') {
	    $seqstring = StringUtils::getRevComp($seqstring);
	}

	# Summarize data e.g >TCONS_00005869 XLOC_001028_-_1:2753268-2784339_Cufflinks
	# and fasta sequence
	# header
	# 		my $tx_biot	=	ExtractFromFeature::getKeyFromFeature($h, $tr, 'transcript_biotype', $verbose);
	my $id      = "$tr";
	my $new_seq = Bio::Seq->new(-id => $id, -seq => $seqstring);

	$seqOUT->write_seq($new_seq);
    }
}

# Print double gtf transcript having the same number of exons/features
sub printDoubleGTF{

    my ($refh1, $refh2, $arefmatch, $mode,  $verbosity)	= @_;
    $mode      //= "all"; # number of fields to be printed
    $verbosity //= 1;

    # 	print Dumper $arefmatch;
    # 	$VAR1 = [
    #           {
    #             'ENSCAFT00000000121' => 'ENSCAFT00000000121',
    #             'ENSCAFT00000000971' => 'ENSCAFT00000000971',
    # 		}
    # 		]
    # print if verbosity	
    print STDERR "\nPrinting transcripts GTF...\n" if ($verbosity > 1);

    foreach my $href (@$arefmatch){
	
	for my $tx1 (keys %$href){
	    
	    my $tx2 = $href->{$tx1};
	    # if not the same number of exons in 2 transcripts...
	    # WARNING : updated
	    # see 		    next if (!Utils::foverlap($feat1->{"start"}, $feat1->{"end"}, $feat2->{"start"}, $feat2->{"end"})); and last
	    warn "warnings: $tx1 has ",scalar (@{$refh1->{$tx1}->{"feature"}}),"feature(s) whereas $tx2 has",scalar (@{$refh2->{$tx2}->{"feature"}}),"feature(s)\n" unless (scalar (@{$refh1->{$tx1}->{"feature"}}) == scalar (@{$refh2->{$tx2}->{"feature"}}));
	    
	    
	    # Printing both hash with indices $j 
	    ##########
	    # 		foreach my $feat1 (@{$refh1->{$tx1}->{"feature"}}) {
	    for (my $j=0; $j < scalar @{$refh1->{$tx1}->{"feature"}}; $j++){
		
		my $feat1 = ${ $refh1->{$tx1}->{"feature"} }[$j];
		
		
		# 1st part
		##########
		print join("\t",$refh1->{$tx1}->{"chr"}, $refh1->{$tx1}->{"source"}, $feat1->{"feat_level"}, $feat1->{"start"}, $feat1->{"end"}, $refh1->{$tx1}->{"score"}, $refh1->{$tx1}->{"strand"}, $feat1->{"frame"});
		print "\tgene_id \"".$refh1->{$tx1}->{"gene_id"}."\"; transcript_id \"$tx1\";";			 			
		
		my @ar= qw/feat_level start end strand frame/;
		
		if ($mode eq "all"){
		    my %tmph = %{$feat1};
		    delete @tmph{@ar};
		    for (keys %tmph){
			print " $_ \"$tmph{$_}\";" if (defined $tmph{$_}); # id defined in order to test if parsinf $extrafield is OK (i.e in case people select FPKM and it does not exists in the file)
		    }
		}				
		
		print "\t";
		# 2nd part
		##########
		
		# Check if different number of exons: At the moment, we only report all exon of A
		# if B tx has less exons, we print simple GTF of A for the extra-line(s) of A
		if ( $j >= scalar @{$refh2->{$tx2}->{"feature"}} ){
		    print "\n";
		    next;
		}
		
		my $feat2 = ${ $refh2->{$tx2}->{"feature"} }[$j];

		print join("\t",$refh2->{$tx2}->{"chr"}, $refh2->{$tx2}->{"source"}, $feat2->{"feat_level"}, $feat2->{"start"}, $feat2->{"end"}, $refh2->{$tx2}->{"score"}, $refh2->{$tx2}->{"strand"}, $feat2->{"frame"});
		print "\tgene_id \"".$refh2->{$tx2}->{"gene_id"}."\"; transcript_id \"$tx2\";";			 			
		
		if ($mode eq "all"){
		    my %tmph = %{$feat2};
		    delete @tmph{@ar};
		    for (keys %tmph){
			print " $_ \"$tmph{$_}\";" if (defined $tmph{$_}); # id defined in order to test if parsinf $extrafield is OK (i.e in case people select FPKM and it does not exists in the file)
		    }
		}	
		print "\n";
	    }			
	} # END 1st HASH
    }
}

# Sort Tx according to start coordinates
# return an array of sorted transcripts
sub sortTxsStartt{

    my ($hreftx, $verbosity) = @_;
    $verbosity //= 1;
    
    die "sortTxs: hash not defined " unless (ref($hreftx));
    
    return sort { $hreftx->{$a}->{'startt'} <=> $hreftx->{$b}->{'startt'} } keys %{$hreftx};
    


}

# Sort Gn according to start coordinates
# return an array of sorted genes
sub sortGnsStartg{

    my ($hrefgn, $verbosity) = @_;
    $verbosity //= 1;
    
    die "sortGnsStartg: hash not defined " unless (ref($hrefgn));
    
    return sort { $hrefgn->{$a}->{'startg'} <=> $hrefgn->{$b}->{'startg'} } keys %{$hrefgn};

}

# Sort Tx exons according to start coordinates
sub sortExons{

    my ($hreftx, $verbosity) = @_;
    $verbosity //= 1;
    
    die "sortExons: hash not defined " unless (ref($hreftx));

    # Sort exons
    for my $t (values %{$hreftx}) {
    	@{$t->{"feature"}} = sort { $a->{"start"} <=> $b->{"start"}} @{$t->{"feature"}};
    }
    return $hreftx;

}

# print hash in GTF format
sub printGTF{
    # parsing a hash in sub need dereference in shift
    my ($refh, $mode, $verbosity) = @_;
    $mode      //= "all"; # number of fields to be printed
    $verbosity //= 0;

    # print if verbosity	
    print STDERR "\nPrinting transcripts GTF...\n" if ($verbosity > 1);
    
    # 	print Dumper $refh;
    # Parse gtfHash to be printed
    for my $tr (keys %{$refh}){
	foreach my $feat1 (@{$refh->{$tr}->{"feature"}}) {
	    print join("\t",$refh->{$tr}->{"chr"}, $refh->{$tr}->{"source"}, $feat1->{"feat_level"}, $feat1->{"start"}, $feat1->{"end"}, $refh->{$tr}->{"score"}, $refh->{$tr}->{"strand"}, $feat1->{"frame"});
	    print "\tgene_id \"".$refh->{$tr}->{"gene_id"}."\"; transcript_id \"$tr\";";
	    
	    if ($mode eq "all"){
		my %tmph = %{$feat1};
		# delete unnecesserary keys from hash %tmph
		delete @tmph{qw/feat_level start end strand frame/};
		for (sort keys %tmph){
		    print " $_ \"$tmph{$_}\";" if (defined $tmph{$_}); # id defined in order to test if parsinf $extrafield is OK (i.e in case people select FPKM and it does not exists in the file)
		}
	    }	
	    print "\n";
	}		
    }
}

sub printGTFrdmTx{
    # parsing a hash in sub need dereference in shift
    my ($refh, $mode, $verbosity) = @_;
    $mode      //= "all"; # number of fields to be printed
    $verbosity //= 0;

    # print if verbosity	
    print STDERR "\nPrinting transcripts GTF...\n" if ($verbosity > 1);
    
    my %seen_gn;
    
    # print Dumper $refh;
    # Parse gtfHash to be printed
    
    for my $tr (keys %{$refh}){

        # increment hash gene counter
        $seen_gn{$refh->{$tr}->{"gene_id"}}++;
        
        # next if already seen gene
        next if ($seen_gn{$refh->{$tr}->{"gene_id"}} > 1 );
	
	
	# print normally as printGTF
	foreach my $feat1 (@{$refh->{$tr}->{"feature"}}) {
	    print join("\t",$refh->{$tr}->{"chr"}, $refh->{$tr}->{"source"}, $feat1->{"feat_level"}, $feat1->{"start"}, $feat1->{"end"}, $refh->{$tr}->{"score"}, $refh->{$tr}->{"strand"}, $feat1->{"frame"});
	    print "\tgene_id \"".$refh->{$tr}->{"gene_id"}."\"; transcript_id \"$tr\";";
	    
	    if ($mode eq "all"){
		my %tmph = %{$feat1};
		# delete unnecesserary keys from hash %tmph
		delete @tmph{qw/feat_level start end strand frame/};
		for (sort keys %tmph){
		    print " $_ \"$tmph{$_}\";" if (defined $tmph{$_}); # id defined in order to test if parsinf $extrafield is OK (i.e in case people select FPKM and it does not exists in the file)
		}
	    }	
	    print "\n";
	}		
    }
}



# print hash in GTF format with transcript line
sub printGTFtx{
    # parsing a hash in sub need dereference in shift
    my ($refh, $mode, $verbosity) = @_;
    $mode      //= "all"; 		# number of fields to be printed
    $verbosity //= 1;

    # print if verbosity	
    print STDERR "\nPrinting transcripts GTF...\n" if ($verbosity > 1);
    
    # Parse gtfHash to be printed
    for my $tr (keys %{$refh}){
	
	my $count=0;
	foreach my $feat1 (@{$refh->{$tr}->{"feature"}}) {

	    $count++;
	    if ($count == 1){
		print join("\t",$refh->{$tr}->{"chr"}, $refh->{$tr}->{"source"}, "transcript", $refh->{$tr}->{"startt"}, $refh->{$tr}->{"endt"}, $refh->{$tr}->{"score"}, $refh->{$tr}->{"strand"}, ".");
		print "\tgene_id \"".$refh->{$tr}->{"gene_id"}."\"; transcript_id \"$tr\";";
		
		if ($mode eq "all"){
		    my %tmph = %{$feat1};
		    # delete unnecesserary keys from hash %tmph
		    delete @tmph{qw/feat_level start end strand frame/};
		    for (sort keys %tmph){
			print " $_ \"$tmph{$_}\";" if (defined $tmph{$_}); # id defined in order to test if parsinf $extrafield is OK (i.e in case people select FPKM and it does not exists in the file)
		    }
		}	
		print "\n";
	    }
	    
	    print join("\t",$refh->{$tr}->{"chr"}, $refh->{$tr}->{"source"}, $feat1->{"feat_level"}, $feat1->{"start"}, $feat1->{"end"}, $refh->{$tr}->{"score"}, $refh->{$tr}->{"strand"}, $feat1->{"frame"});
	    print "\tgene_id \"".$refh->{$tr}->{"gene_id"}."\"; transcript_id \"$tr\";";
	    
	    if ($mode eq "all"){
		my %tmph = %{$feat1};
		# delete unnecesserary keys from hash %tmph
		delete @tmph{qw/feat_level start end strand frame/};
		for (sort keys %tmph){
		    print " $_ \"$tmph{$_}\";" if (defined $tmph{$_}); # id defined in order to test if parsinf $extrafield is OK (i.e in case people select FPKM and it does not exists in the file)
		}
	    }	
	    print "\n";
	}		
    }
}


# MergeSmallIntrons
sub MergeSmallIntrons{

    my ($refh,  $minintron, $fh, $verbosity) = @_;
    $minintron //= 5;
    $fh        //= undef;
    $verbosity //= 0;


    # Test parsing i.e empty hash
    die "ExtractFromHash::MergeIntrons => ref not defined...\n" if (not defined $refh);
    my %h = %{$refh};
    
    # Filter on exon only
    print STDERR "- Get exon levels first...\n" if ($verbosity > 1);
    %h = filterGtfHash(\%h, 'exon', $verbosity);
    
    # Test parsing i.e empty hash
    die "Error! getIntronFromGtfHash return empty hash ...\nCannot compute introns Check that your infile contains 'exon' levels\n" unless (scalar(keys(%h))>0);

    # Get introns from has
    for my $tr (keys %h){
	# next if monoexonic
	next if (scalar (@{$h{$tr}->{"feature"}}) ==1);

	my %mergeexon;
	my @exons= @{$h{$tr}->{"feature"}};

	
	for (my $i=1; $i< @exons ; $i++) {
	    
	    last unless  ($i < @exons); # WARNING : since we reedit the @ we need to check it within the loop	: see http://stackoverflow.com/questions/6537838/loop-through-two-arrays-deleting-overlaps-in-perl	
	    last if (scalar(@exons) ==1); # If merging produce a  monoexonic
	    
	    # comute intron charac
	    my $intron_start = $exons[$i-1]->{"end"} + 1;
	    my $intron_end   = $exons[$i]->{"start"} - 1;
	    my $intron_size  = $intron_end - $intron_start + 1;
	    if ($intron_size < $minintron){
		if ($fh){
		    print $fh "Merging intron for $tr ",$h{$tr}->{"chr"},":",$intron_start-1,"-",$intron_end+1," (",$h{$tr}->{"strand"},")...\n";
		}else {
		    print STDERR "Merging intron for $tr ",$h{$tr}->{"chr"},":",$intron_start-1,"-",$intron_end+1," (",$h{$tr}->{"strand"},")...\n";
		}					
		# we reassign the start of exon ith-1 to exon ith
		$exons[$i]->{"start"} = $exons[$i-1]->{"start"};
		# we delete exon i
		splice @exons, $i-1,1;
		# relaunch with new exon boundaries
		redo;
	    }
	}
	
	# Reassign new array intron in place of exon
	@{$h{$tr}->{"feature"}} = @exons;
	
    }
    return %h;		
}



# Filter gff hash according to a string of level: it could be : exon, mRNA, transcript, intron, gene, CDS...
# 	#######################################################################################################
sub filterGtfHash{

    # parsing a hash in sub need dereference in shift
    my %h           = %{shift()};
    my $level       = shift;	
    my $verbosity   = shift;
    $verbosity    //= 1;
    
    # split levels accroding to ,
    my @all_level = split (",", $level);
    
    # 	my @authorized_levels = ('exon' ,'mRNA' , 'transcript', 'intron', 'gene', 'CDS', 'start_codon', 'stop_codon', 'UTR');
    die "Error! Level filter '$level' invalid...\n" unless ($level eq 'exon' || $level eq 'mRNA' || $level eq  'transcript'|| $level eq  'intron'|| $level eq  'gene'|| $level eq  'CDS'|| $level eq  'start_codon'|| $level eq  'stop_codon'|| $level eq  'UTR');	
    
    # print if verbosity	
    print STDERR "- Filter GTF file for '$level' level...\n" if ($verbosity > 1);
    
    # Where to store only filtered feature
    my @filtered_feature;
    
    for my $tr (keys %h){
	# $feat1 is a hash 
	# We push in a filtered array 'filtered_feature' only for good '$level'
	foreach my $feat1 (@{$h{$tr}->{"feature"}}) {
	    if ($feat1->{'feat_level'} eq $level) {
		push (@filtered_feature, $feat1);
	    }
	}
	# we then re assign the filtered array to the array ref $h{$tr}->{"feature"}
	@{$h{$tr}->{"feature"}} = @filtered_feature;
	
	# we reinitialize the filtred array for next transcript
	@filtered_feature = ();
    }
    # Warn if empty
    if (keys (%h) == 0){
	die "Error! Filter '$level' returns an empty elements...\n";
    }else {
	return %h;
    }
}


# Return a hash of introns
sub getIntronFromGtfHash{

    # parsing a hash in sub need dereference in shift
    my %h           = %{shift()};
    my $verbosity   = shift;
    $verbosity    //= 1;

    # print if verbosity	
    print STDERR "- Extract Introns levels...\n" if ($verbosity > 1);	
    
    # hash for monoexonique transcript i.e no intron possible
    my %h_countIntron;
    
    # Filter on exon only
    print STDERR "- Get exon levels first...\n" if ($verbosity > 1);
    %h = filterGtfHash(\%h, 'exon', $verbosity);
    
    # Test parsing i.e empty hash
    die "Error! getIntronFromGtfHash return empty hash ...\nCannot compute introns Check that your infile contains 'exon' levels\n" unless (scalar(keys(%h))>0);
    
    # hash that will contain intron data
    my %feature_intron ;
    # array that will contain hash of intron
    my @intron ;

    # Get introns from has
    for my $tr (keys %h){
	# if monoexonic
	if (scalar (@{$h{$tr}->{"feature"}}) ==1){
	    $h_countIntron{$tr}++;
	    delete $h{$tr};
	} else {
	    # If pluriexoniq
	    # we loop on every hash  in the ref array
	    # and get the exon end value + 1 = start of intron
	    # and get the exon start value - 1 = end of intron
	    for (my $i=1; $i<scalar(@{$h{$tr}->{"feature"}}); $i++) {
		# 			foreach my $feat1 (@{$h{$tr}->{"feature"}}) {
		my $intron_start	=	${$h{$tr}->{"feature"}}[$i-1]->{"end"} + 1;
		my $intron_end		=	${$h{$tr}->{"feature"}}[$i]->{"start"} - 1;
		my $intron_size		= 	$intron_end - $intron_start + 1;
		my %feature_intron = ( "start" => $intron_start, "end" => $intron_end, "feat_level" => 'intron', "intron_size" => $intron_size);
		push (@intron, \%feature_intron);
	    }
	    # Reassign new array intron in place of exon
	    @{$h{$tr}->{"feature"}} = @intron;
	    
	    # reinitialize intron array
	    @intron = ();			
	}
    }
    warn "# ".keys(%h_countIntron)." monoexonic transcript(s) (i.e without intron)\n" if ($verbosity>0);
    return %h;		
}




sub getTranscriptFromGtfHash{

    # parsing a hash in sub need dereference in shift
    my %h           = %{shift()};
    my $verbosity   = shift;
    $verbosity    //= 1;

    # print if verbosity	
    print STDERR "- Extract Transcript levels...\n" if ($verbosity > 1);
    

    # Filter on exon only
    print STDERR "- Get exon levels first...\n" if ($verbosity > 1);
    %h = filterGtfHash(\%h, 'exon', $verbosity);
    
    # Test parsing i.e empty hash
    die "Error! getPromoterFromGtfHash return empty hash ...\nCannot compute  transcript levels...\nCheck that your infile contains 'exon' levels\n" unless (scalar(keys(%h))>0);
    
    my @transcript=();
    
    # Get transcript from hash
    for my $tr (keys %h){
	my %feature_transcript = ( "start" => $h{$tr}->{"startt"}, "end" => $h{$tr}->{"endt"}, "feat_level" => 'transcript', "extrafield" => undef);
	push (@transcript, \%feature_transcript);

	# Reassign new array transcript (with 1 element) in place of exon feature
	@{$h{$tr}->{"feature"}} = @transcript;
	
	# reinitialize intron array
	@transcript = ();			
    }
    return %h;		
}

# Get the gene level from a hash of transcripts
sub getGeneFromGtfHash{

    # parsing a hash in sub need dereference in shift
    my %h           = %{shift()};
    my $verbosity   = shift;
    $verbosity    //= 1;

    # print if verbosity	
    print STDERR "- Extract Gene levels...\n" if ($verbosity > 1);

    # Filter on exon only
    print STDERR "- Get exon levels first...\n" if ($verbosity > 1);
    %h = filterGtfHash(\%h, 'exon', $verbosity);
    
    # Test parsing i.e empty hash
    die "Error! getGeneFromGtfHash return empty hash ...\nCannot compute gene levels...\nCheck that your infile contains 'exon' levels\n" unless (scalar(keys(%h))>0);
    
    #new hash
    my %h_gene ;
    
    for my $tr (keys %h){

	my $locus  = $h{$tr}->{"gene_id"};
	my $beg_ex = $h{$tr}->{"startt"};
	my $end_ex = $h{$tr}->{"endt"};
	
	# reassignment to a new hash of gene
        $h_gene{$locus}->{"chr"}     = $h{$tr}->{"chr"};
        $h_gene{$locus}->{"source"}  = $h{$tr}->{"source"};
        $h_gene{$locus}->{"strand"}  = $h{$tr}->{"strand"};
        $h_gene{$locus}->{"score"}   = ".";        
        $h_gene{$locus}->{"frame"}   = ".";
        $h_gene{$locus}->{"gene_id"} = $h{$tr}->{"gene_id"};        
        $h_gene{$locus}->{"startt"}  = Utils::min2($h_gene{$locus}->{"startt"}, $beg_ex);
        $h_gene{$locus}->{"endt"}    = Utils::max2($h_gene{$locus}->{"endt"}, $end_ex);
        
        
        # to be compatible with method "printGtfHash", we have to put in feature
	my %feature = ( "start" => $h_gene{$locus}->{"startt"}, "end" => $h_gene{$locus}->{"endt"}, "feat_level" => 'gene', "extrafield" => undef);
	
	# Warning since, there is only one refarray for a gene, we dont want to push 
	@{$h_gene{$locus}->{"feature"}} = \%feature;
    }
    return %h_gene;		
}


# Get promoter region with a size of $size
sub getPromoterFromGtfHash{

    # parsing a hash in sub need dereference in shift
    my %h           = %{shift()};
    my $size        = shift;
    my $verbosity   = shift;
    $size         //= 500;
    $verbosity    //= 1;
    
    # print if verbosity	
    print STDERR "- Extract '".$size."'bp Promoter levels...\n" if ($verbosity > 1);
    

    # Filter on exon only
    print STDERR "- Get exon levels first...\n" if ($verbosity > 1);
    %h = filterGtfHash(\%h, 'exon', $verbosity);
    
    # Test parsing i.e empty hash
    die "Error! getPromoterFromGtfHash return empty hash ...\nCannot compute  transcript levels...\nCheck that your infile contains 'exon' levels\n" unless (scalar(keys(%h))>0);
    
    my @promoter=();
    
    # Get promoter from hash
    for my $tr (keys %h){
	# intialize feature promoter
	my %feature_promoter;

	# Promoter region is dependant on strand
	if ($h{$tr}->{"strand"} eq "+" || $h{$tr}->{"strand"} eq "."){
	    %feature_promoter = ( "start" => $h{$tr}->{"startt"}-$size, "end" => $h{$tr}->{"startt"}, "feat_level" => 'promoter', "extrafield" => 'size "'.$size.'";');
	    push (@promoter, \%feature_promoter);
	} elsif ($h{$tr}->{"strand"} eq "-"){
	    %feature_promoter = ( "start" => $h{$tr}->{"endt"}, "end" => $h{$tr}->{"endt"} + $size, "feat_level" => 'promoter', "extrafield" => 'size "'.$size.'";');
	    push (@promoter, \%feature_promoter);			
	}
	# Reassign new array promoter in place of exon
	@{$h{$tr}->{"feature"}} = @promoter;
	
	# reinitialize intron array
	@promoter = ();			
    }
    return %h;
}





sub printGtfHash_Line{
    # parsing a hash in sub need dereference in shift
    my ($refh, $tx, $ex, $tab, $verbosity) = @_;
    my %h        = %{$refh}; # parsing a hash in sub need dereference in shift
    $tx        //= '';
    $ex        //= '';
    $tab       //= 0;
    $verbosity //= 1;

    # if transcript is not in hash
    die "$tx not in hash..." unless exists ($h{$tx});
    
    # print if verbosity	
    print STDERR "Printing printGtfHash...\n" if ($verbosity > 1);
    
    # Parse gtfHash to be printed
    foreach my $feat1 (@{$h{$tx}->{"feature"}}) {
	if ($feat1 eq $ex){
	    print join("\t",$h{$tx}->{"chr"}, $h{$tx}->{"source"}, $feat1->{"feat_level"}, $feat1->{"start"}, $feat1->{"end"}, $h{$tx}->{"score"}, $h{$tx}->{"strand"}, $h{$tx}->{"frame"});
	    print "\tgene_id \"".$h{$tx}->{"gene_id"}."\"; transcript_id \"$tx\";";
	    if (defined($feat1->{"extrafield"})){
		print " ",$feat1->{"extrafield"};
	    }
	    if (defined ($h{$tx}->{"size"})){
		print " feature_size \"",$h{$tx}->{"size"},"\";";
	    } else { 
		if ($tab){print "\t";}else{print "\n";}
	    }
	}
    }	
}

# get the cumulative size of features from a array reference
sub cumulSize{
    my $refonarray   = shift;	
    my $verbosity    = shift;
    $verbosity     //= 1;

    # print if verbosity	
    print STDERR "Getting cumulative size...\n" if ($verbosity > 1);
    
    # cumulsize
    my $cumulsize = 0;
    
    # Parse gtfHash to be printed
    foreach my $feat1 (@{$refonarray}) {

	$cumulsize += ($feat1->{"end"} - $feat1->{"start"} +1) if ($feat1->{'feat_level'} eq "exon");

    }
    return $cumulsize;
}

# Get cumulative size of feature from a gtf 
# Option : $longest = BOOL to return only longest per locus (default 0)
sub getCumulSizeFromGtfHash{

    # parsing a hash in sub need dereference in shift
    my ($refh,  $verbosity, $longest) = @_;
    $verbosity //= 1;
    $longest   //= 0;

    my %store_size;
    
    # print if verbosity	
    print STDERR "- Extract size from GtfHash (longest option set to $longest) ...\n" if ($verbosity > 1);	
    
    # Mode : only compute size of feature 
    if ($longest == 0 ) {
	
	# Get size of the array ref and assign to transcript
	for my $tr (keys %{$refh}){
	    
	    #compute size
	    my $size               = cumulSize($refh->{$tr}->{"feature"});
	    $refh->{$tr}->{"size"} = $size;
	}
	return $refh;

    } else { # MODE: get longest transcript per locus

	# print if verbosity	
	print STDERR "- Extract Longest feature per gene_id from GtfHash ...\n" if ($verbosity > 1);
	
        my $refgene;
        my %longest_tx;
	my $i=0;
	# Warning : sort keys to always get the same tx for reproducibility due to hash random storing and 2 isoforms of a same size
        for my $tr (sort keys %{$refh}){

	    if ($verbosity > 0) {
		Utils::showProgress(scalar(keys(%{$refh})), $i++, "Get longest Tx per gene: ");
	    }
	    
	    # Get tx infos
            my $size    = ExtractFromHash::cumulSize($refh->{$tr}->{"feature"});
            my $gene_id = $refh->{$tr}->{"gene_id"};

            # initialize gene hash size
            $refgene->{$gene_id}->{"size"} = 0 if (!exists ($refgene->{$gene_id}->{"size"} ));

            if ($size > $refgene->{$gene_id}->{"size"} ){
		$refgene->{$gene_id}->{"size"} = $size;
		$longest_tx{$gene_id}          = $tr;
            }
        }
	print STDERR "\n" if ($verbosity > 0); # because of the showProgress

	# map values 
        my %new_hash = map { $_ => $refh->{$_} } values(%longest_tx);

        return \%new_hash;
	
    }
}

# Get Transcripts with a size ge $min and le $max
sub getMinMaxTxSizeFromGtfHash{

    my ($refh, $min, $max, $verbosity) = @_;
    my %h        = %{$refh};
    $min       //= 200;    # min: default value is 200bp
    $max       //= 100000; # max: default value is 100,000bp
    $verbosity //= 1;
    

    # print if verbosity    
    print STDERR "- Extract transcripts with  '$min' < tx size > '$max'...\n" if ($verbosity > 1);

    # Test parsing i.e empty hash
    die "Error! Parsing getMinMaxTxSizeFromGtfHash return empty hash...\nCheck that your infile contains 'exon' levels\n" unless (scalar(keys(%h))>0);

    for my $tr (keys %h){

	#compute size
	my $size          = cumulSize($h{$tr}->{"feature"});
	$h{$tr}->{"size"} = $size;
	
	if ($size < $min){
	    delete $h{$tr};
	} elsif ($size > $max) {
	    delete $h{$tr};
	}
    }
    # Test parsing i.e empty hash
    die "Error! getMinMaxExonFromGtfHash return empty hash...\nCheck that min ='$min' and max='$max' size limits are valid\n" unless (scalar(keys(%h))>0);

    # Return %h
    return %h;
}



sub minStartFromRefArray{
    my $refonarray   = shift;	
    my $verbosity    = shift;
    $verbosity     //= 1;

    # print if verbosity	
    print STDERR "Getting min value from refArray...\n" if ($verbosity > 1);
    
    # min intialize to very big value
    my $min = 10000000000;
    
    # Parse gtfHash to be printed
    foreach my $feat1 (@{$refonarray}) {
	if ($feat1->{"start"} < $min){
	    $min = $feat1->{"start"};
	}
    }
    return $min;
}

sub maxEndFromRefArray{
    my $refonarray   = shift;	
    my $verbosity    = shift;
    $verbosity	   //= 1;

    # print if verbosity	
    print STDERR "Getting max value from refArray...\n" if ($verbosity > 1);
    
    # max intialize to very small value
    my $max = 0;
    
    # Parse gtfHash to be printed
    foreach my $feat1 (@{$refonarray}) {
	if ($feat1->{"end"} > $max){
	    $max = $feat1->{"end"};
	}
    }
    return $max;
}




# Get Transcripts with a number of exon ge $min and le $max
sub getMonoExonicFromGtfHash{

    # parsing a hash in sub need dereference in shift
    my %h = %{shift()};

    # Test parsing i.e empty hash
    die "Error! Parsing getMinMaxExonFromGtfHash return empty hash...\nCheck that your infile contains 'exon' levels\n" unless (scalar(keys(%h))>0);

    for my $tr (keys %h){
	if (scalar (@{$h{$tr}->{"feature"}}) != 1){
	    delete $h{$tr};
    	}
    }

    # Return %h
    return %h;
}



# Get Transcripts with a number of exon ge $min and le $max
sub getMinMaxExonFromGtfHash{

    # parsing a hash in sub need dereference in shift
    my %h         = %{shift()};
    my $min       = shift;
    my $max       = shift;
    my $verbosity = shift;
    $min          //= 1; # default monoexonic
    $max          //= 1; # default monoexonic
    $verbosity      = 1;

    # print if verbosity    
    print STDERR "- Extract transcripts with  '$min' < nb exon(s) > '$max'...\n" if ($verbosity > 1);

    # Test parsing i.e empty hash
    die "Error! Parsing getMinMaxExonFromGtfHash return empty hash...\nCheck that your infile contains 'exon' levels\n" unless (scalar(keys(%h))>0);

    for my $tr (keys %h){
	if (scalar (@{$h{$tr}->{"feature"}}) < $min){
	    delete $h{$tr};
	} elsif (scalar (@{$h{$tr}->{"feature"}}) > $max) {
	    delete $h{$tr};
	}
    }

    # Return %h
    return %h;
}

# Get N random Transcripts 
sub getNTxFromGtfHash{

    # parsing a hash in sub need dereference in shift
    my %h         = %{shift()};
    my $nbtx      = shift;
    my $verbosity = shift;
    $nbtx         //= 1; # default 1 tx
    $verbosity    //= 1;

    # print if verbosity    
    print STDERR "- Extract '$nbtx' random transcripts...\n" if ($verbosity > 1);

    # counter for nb tx
    my $count=0;

    # newh
    my %newh;

    for my $tr (shuffle keys %h){
	$count++;
	$newh{$tr} = $h{$tr};
	last if ($count ==  $nbtx);
    }

    # return good one(s)
    return %newh;
}


1;
