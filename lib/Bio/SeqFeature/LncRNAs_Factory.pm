package Bio::SeqFeature::LncRNAs_Factory;
use strict;
use warnings;
=head2 package : LncRNAs_Factory


Bio::SeqFeature::InteractionCollection

Use : Interaction , InteractionCollection

=head1 ABOUT
Authors: Audrey DAVID - M2 Informatique opt Bioinfo Nantes
		 Fabrice LEGEAI - IRISA/INRIA Rennes / INRA Rheu
		 Thomas	 DERRIEN - CNRS Rennes
If you have questions contact authors.

=head1 DESCRIPTION

Implement several function to have an classification collection of interaction come from 2 file (lncRNA and mRNA)
Implementation of classification rules of Derrien et al 2012.

=cut
use Bio::SeqFeature::database_part;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Interaction;
use Bio::SeqFeature::InteractionCollection;
use Bio::SeqFeature::Genic;
use Bio::SeqFeature::InterGenic;
use Bio::SeqFeature::Empty;

use Utils; # guess_format

use Carp;



=cut
 NOTICE: TO extract a part of the DB wich are specific tag source or\/and specific position
 		 please use the Bio::DB::SeqFeature::Store methods (get_seq_stream for example)
=cut

=head1 type_of_interaction
	return 1 (genic) or 0 (intergenic)
	object : SeqFeature element
	subject : SeqFeature element
=head1 DESCRIPTION
	what is the type of the interaction shared by an object SeqFeature element and a subject SeqFeature element?
	Does the interaction is genic or intergenic?
=cut

sub type_of_interaction{
	shift; #pkg
	my $object = shift;
	my $subject = shift;

	#control if both are SeqFeature object
	if ( _arethey_BioSeqFeatureI($object, $subject) == 0 ) {
		croak (" Type of interaction :  SeqFeatures objects are required to make this action . ", $object, "  sub ", $subject," \n");
	}

	# overlapping or not?
	if ($subject->overlaps($object) == 1){
		return(1); # genic
	}else {
		return(0); #intergenic
	}
}

# private method : _no_interaction
# function : design to answer to : what happens if an lncRNA do not have any interaction? -> report it

sub _no_interaction{
	shift;
	my $object = shift;
    my $interaction;
 	my $subject = Bio::SeqFeature::Generic->new(-seq_id =>  '-555', -start => 0, -end => 1 );

		# create an interaction intergenic interaction
		$interaction = Bio::SeqFeature::Empty->new('-object' => $object , '-subject' => $subject);
		#$interaction->printer();
	#sleep(1);

	return $interaction;
}


=cut
=head1 create_interaction
	return : Bio::DB::Interaction
=head1 PARAMETERS
	object : SeqFeature element
	subject : SeqFeature element
=head1 DESCRIPTION
	create an interaction
=cut

sub create_interaction{
	shift;
	my $object = shift;
	my $subject = shift;
	my $interaction;

	my $type = type_of_interaction("pkg",$object,$subject);

	if ( $type == 1 ){


		# create an genic interaction
		$interaction = Bio::SeqFeature::Genic->new('-object' => $object, '-subject' => $subject);
		#$interaction->warning();
	}elsif ($type == 0){

		# create an interaction intergenic interaction
		$interaction = Bio::SeqFeature::InterGenic->new('-object' => $object, '-subject' => $subject);
		#$interaction->printer();
	}else{
		croak ("Create interaction : Unrecognized TYPE. Please contact authors \n");
	}

	return $interaction;
}

=cut
	NOTICE : Interactions are added to an interaction collection
			 Refer to Bio::SeqFeature::InteractionCollection
=cut

=cut
=head1 DoItForMe
	return Bio::DB::InteractionsCollection
=head1 PARAMETERS
	- user
	- host
	- database
	- type (mRNA , transcript etc.)
	- size of window (pb)
        - biotype: if the user want the biotype in the output
	optional but going together:
	two files :
				 - a gtf file from bedtools merged containing your LncRNAs ;
				 - a gff file containing your mRNAs.
	Then will classified each of them in relation to each other
=cut

sub DoItForMe{

## first step database import
	my $pkg = shift;
	my $window = shift;
	my $lncrna_file = shift;
	my $mrna_file = shift;
	my $max_window = shift;
	my $biotype = shift;
	my $verbosity = shift;

	print STDERR "window : $window - max window : $max_window - lncrna : $lncrna_file - mrna : $mrna_file - biotype: $biotype\n";
	my $db_lnc = Bio::SeqFeature::database_part->new();
	my $db_mrna = Bio::SeqFeature::database_part->new();


	# test input lnc format
	my $lncformat;
	if (Utils::guess_format($lncrna_file) eq 'gtf'){
		$lncformat = 2;
	} elsif (Utils::guess_format($lncrna_file) eq 'gff'){
		$lncformat = 3;
	} else {
		die "Cannot recognize lncRNA input format '$lncrna_file' \n(supported extensions are : gtf, gff, gff2, gff3)\n";
	}

	# test input mrna format
	my $mrnaformat;
	if (Utils::guess_format($mrna_file) eq 'gtf'){
		$mrnaformat = 2;
	} elsif (Utils::guess_format($mrna_file) eq 'gff'){
		$mrnaformat = 3;
	} else {
		die "Cannot recognize mrna input format '$mrna_file' \n(supported extensions are : gtf, gff, gff2, gff3)\n";
	}


	# load files
	my $nblnc = $db_lnc->load_merge_gtf_into_db($lncrna_file,$lncformat, 'lncRNA');
	my $nb_mrna = $db_mrna->load_merge_gtf_into_db($mrna_file, $mrnaformat,'mRNA');


	my $prec = 0;

	my $collection =  Bio::SeqFeature::InteractionCollection->new();

	my $nombre = 0;
	my $numberlnc = 0;
	my @no_interaction=();

	my $lnc=$db_lnc->get_seq_stream(-type => "lncrna");


	while (my $object = $lnc->next_seq) {
		my $step=1;
		my $interaction_per_lncRNA = 0;
		$numberlnc++;

		while (($interaction_per_lncRNA) == 0 && ($step*$window <= $max_window)) {
		    if ($step >1) {
			if($verbosity >= 2) {
			    print  STDERR "No interaction found for ", $object->get_tag_values("transcript_id"), " within a window  of size  ", $window*($step-1), ". Expanded to : ", $window*$step,"\n";
			}
		    }

			my $autour_s = $object->start() - $step*$window; # upstream
			my $autour_e = $object->end() + $step*$window; #downstream

 			 my @features = $db_mrna->get_features_by_location(-seq_id=>$object->seq_id(),-start=>$autour_s,-end=>$autour_e);

			foreach my $subject (@features){
				next unless ($subject->primary_tag eq 'mRNA');
				$nombre ++;
				$interaction_per_lncRNA++;
				my $interaction = Bio::SeqFeature::LncRNAs_Factory->create_interaction($object, $subject);

				#$interaction->printer();

				## sixth step : add to the collection
				$collection->add_interaction($interaction);
			}

 		 	if ($interaction_per_lncRNA == 0) {
			 	# We try to expand the environment
				$step++;
 			}


		    if ($nombre%100 ==0 && $prec !=$nombre ){
			if($verbosity >= 1){
			    print STDERR " ... Looking for interactions, currently eq to $nombre \n" ;
			}
				$prec = $nombre;
			}

		}
		if ($interaction_per_lncRNA == 0) {
			# We try to expand the environment
			push @no_interaction,  $object->get_tag_values("transcript_id");
 		}

	}


	# Make a log file
	my $logName = Utils::renamefile($lncrna_file, ".feelncclassifier.log");

	# Open the file
	open FILE, "> $logName" or die "Error! Cannot access log file '". $logName . "': ".$!;

	print FILE "#FEELnc Classification\n";
	print FILE "#lncRNA file :  lncrna : $lncrna_file \n";
	print FILE "#mRNA file : $mrna_file\n";
	print FILE "#Minimal window size : $window\n";
	print FILE "#Maximal window size : $max_window\n";
	print FILE "#Number of lncRNA : $numberlnc \n";
	print FILE "#Number of mRNA : $nb_mrna\n";
	print FILE "#Number of interaction : $nombre \n";
	print FILE "#Number of lncRNA without interaction : ",scalar(@no_interaction), "\n";
	print FILE "#List of lncRNA without interaction : ", join (" ", @no_interaction), "\n";

	close FILE;

	$db_lnc->destroy();
	$db_mrna->destroy();
	## seventh step: return the collection

	return $collection;

}

=cut
=head1 DoIt_with2gtf
	return Bio::DB::InteractionsCollection
=head1 DESCRIPTION

the same than DoItfor me but with 2 gtf
=cut
sub DoIt_with2gtf{


## first step database import

	if (scalar(@_) > 2 ) {
		shift; #pkg
		my $gtf_file = shift;
		my $gtf_file2 = shift;
		my $db = Bio::SeqFeature::database_part->initialize_db();

		$db = Bio::SeqFeature::database_part->load_complete_gtf_into_db($db, $gtf_file, 'lncRNA');
		$db = Bio::SeqFeature::database_part->load_complete_gtf_into_db($db, $gtf_file2, 'mRNA');

	}


	my $db = Bio::SeqFeature::database_part->open_db();

## second step : collection initialization

	my $collection =  Bio::SeqFeature::InteractionCollection->new();

	my $nombre = 0;

## third step : extract for the DB the lncRNAs - the tag primary is lncRNA

	my $lnc=$db->get_seq_stream(-primary  => "lncRNA");

	while (my $object = $lnc->next_seq) {
		my $autour_s = $object->start() - 10000; # upstream
		my $autour_e = $object->end() + 10000; #downstream

## forth step : for each give me the mRNA wich are around : $autour_s and $autour end

		my $aphid=$db->get_seq_stream(-seq_id => $object->seq_id() , -type => "mRNA", -start => $autour_s ,-end => $autour_e );

## fifth step : create the interaction
		while (my $subject = $aphid->next_seq()){
			$nombre ++;
			my $interaction = Bio::SeqFeature::LncRNAs_Factory->create_interaction($object, $subject);

## sixth step : add to the collection

			$collection->add_interaction($interaction);
		}

		#sleep(5);
	}
## seventh step: return the collection
		return $collection;

}


# open a file ... return an file handle or crash your program if the file can't be opened
# be sure that your close your handle after usage
sub _openFile {

}

# # private method : _arethey_BioSeqFeatureI
# control if 2 objects are Bio::SeqFeature
# parameter(s) : $features1, $features2
# return : 1 or 0

sub _arethey_BioSeqFeatureI{
	my $object = shift;
	my $subject = shift;

	#control if both are SeqFeature object
	if (_isDefined ($object) == 1 && $object->isa('Bio::SeqFeatureI')){
		if (_isDefined ($object) == 1 && $subject->isa('Bio::SeqFeatureI')){
			return 1;
		}
	}
	return 0;
}

# # private method : _isDefined
# function : control if something is defined
# parameter(s) : variable
# return : 1 or 0

sub _isDefined{
	my $val = shift;
	if (defined $val) { return 1;} else { return 0;}
}

# # private method : _isDefined
# function : print exons information
# parameter(s) : feature
# return : /

sub _print_exons{
	my $feat = shift;

	my @feats = $feat->get_SeqFeatures();

	foreach my $f (@feats){
		print "\t\t  Exon  \n";
		_afficher_fvalues($f);
	}
}

# # private method : _ligne_carre
# function : a pretty design line ...

sub _ligne_carre(){
	print "\t ------\t------\t------\n";
}
# # private method : _afficher_fvalues
# function : print features information
# parameter(s) : @features
# return : /

sub _afficher_fvalues{
	my $feat = shift @_;
	print " \t ______________________________ \n ";
	print " \t     feature information \n ";
	print " \t ______________________________ \n ";

		print "\t scaffold / seq_id:",  $feat->seq_id() , "\n";
   		print "\t start :",  $feat->start() , "\n";
   		print "\t end :" , $feat->end() , "\n";
 		print "\t source_tag:" , $feat->source_tag() , "\n";
   		print "\t \t";
   		print " hash tag  -tag    => {clé => val}\n";
   		print "\t \t";print " class code : ", $feat->get_tag_values("class_code"), "\n";
   		print "\t \t";print " gene id : ", $feat->get_tag_values("gene_id"), "\n";
    	print "\t \t";print " oId : ", $feat->get_tag_values("oId"), "\n";
   		print "\t \t";print " tss_id: ", $feat->get_tag_values("tss_id"), "\n";
    	print "\t \t";print " transcript_id: ", $feat->get_tag_values("transcript_id"), "\n";
    		print " \t ______________________________ \n";
}

1;
