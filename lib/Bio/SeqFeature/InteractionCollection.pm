package Bio::SeqFeature::InteractionCollection;
use Bio::SeqFeature::Generic;
use Carp;
use strict;
use warnings;
=head2 CLASS: InteractionCollection

SUPER CLASS : None

Bio::SeqFeature::InteractionCollection

Use : Interaction ; Iterator

=head1 ABOUT
Authors: Audrey DAVID - M2 Informatique opt Bioinfo Nantes
		 Fabrice LEGEAI - IRISA/INRIA Rennes / INRA Rheu
		 Thomas	 DERRIEN - CNRS Rennes
If you have questions contact authors.

=head1 DESCRIPTION

Implement several method to manipulate Collection of interaction
Implementation of classification rules of Derrien et al 2011.
= head1 Others
Usually many of these method as Method has to be used like :

my $res = $col->Method();

some of method can not return result than it's a printer and has to be use like
$col->Method();

Please check the description

=cut

use Bio::SeqFeature::Interaction;
use Bio::SeqFeature::InteractionIterator;

=cut
=head1 new()

Bio::SeqFeature::InteractionCollection

=head1 DESCRIPTION

my $coll = Bio::SeqFeature::InteractionCollection->new(@interactions)
contains a list of subjects index and  a list of objects index

=cut


sub new {
	my $pkg = shift;
  	my $coll = bless {
	 _subjects => undef,
	 _objects => undef,
	}, $pkg;


	foreach my $interaction (@_) {
		#print $interaction," j'ajoute une interaction\n";

		if ($interaction->isa('Bio::SeqFeature::Interaction')) {
			$coll->add_interaction($interaction);# add the interaction to a pointer to the object seqfeature object
			# add the interaction to a pointer to the subject seqfeature object
		}
		else {
			carp("you have to give a list of Bio::SeqFeature::Interaction to InteractionCollection");
		}
	}

  return $coll;
}





=head1 add_interaction

Bio::SeqFeature::Interaction

=head1 DESCRIPTION

add an Interaction or a list of interaction in a collection of interactions

=cut

sub add_interaction{
	my $self = shift;
	my $interaction = shift;
	if (_isComplete($interaction) == 0){
		carp ("To be add to a collection of interactions, the interaction has to be complete \n");
	}else{

		$self->_add_interaction_subject($interaction);
		$self->_add_interaction_object($interaction);
	}
}


### add interaction as object

sub _add_interaction_object{
	my $self = shift;
	my $interaction = shift;
	my $object_indx = $interaction->object(); #index


	push (@{$self->{'_objects'}->{$object_indx}}, $interaction); # add the interaction to a pointer to the object seqfeature object
}

### add the subject of the interaction

sub _add_interaction_subject{
	my $self = shift;
	my $interaction = shift;
	my $subject_indx = $interaction->subject();

	push (@{$self->{'_subjects'}->{$subject_indx}}, $interaction);
}



=head1 get_subjects

returns a list of Bio::SeqFeature

=head1 DESCRIPTION

this function returns the list of subjects of interactions

=cut

sub get_subjects {
	my $self = shift;

	return  keys %{$self->{'_subjects'}};

}


=head1 get_objects

returns a list of Bio::SeqFeature

=head1 DESCRIPTION

this function returns the list of objects of interactions

=cut

sub get_objects {
	my $self = shift;

	return  keys %{$self->{'_objects'}};

}

=head1 get_interaction

Bio::SeqFeature

=head1 DESCRIPTION

get all Interaction share by an seqfeature

the seqfeature can be
	 object => (1),
	 subject => (2),
	 both => (3).

by defeault the seqfeature can be all the second parameter is 3

if the seqfeature is the object fill 1 for the second parameter
if the seqfeature is the subject fill 2 for the second parameter


$col->get_interaction($feature, 2 ); #look for all interaction shared by $feature and in which $feature is the subject

note:  if you pass an unrecognized second argument, than the default one will be called (means : all )

=cut

sub get_interactions{
	my $col = shift;
	my $index = shift; # adress of the subject/object
	my $option = 3; # by default all
	my $iterator = Bio::SeqFeature::InteractionIterator->new();
	if (@_)	{
		$option = shift; # maybe you want to change default argument
	}
	if ( $option == 2 ){ #want to get all interaction in which the index is the subject
		$iterator->add( $col->_get_interactions_subject($index ) );
	}elsif ($option == 1){ #want to get all interaction in which the index is the object


		$iterator->add( $col->_get_interactions_object($index ) );
	}else{ #both
		$iterator->add( $col->_get_interactions_object($index) ,$col->_get_interactions_subject($index ) );
	}

	return $iterator;
}


# private method : _getArray_from_hash
# function : get an array value from a hash with the key
# hash looks like : %hash{$key}=@array

sub _getArray_from_hash {
	my %hash = %{$_[0]};
	my $cle = $_[1];


	if ( _isExist_hash (\%hash, $cle) == 0 ){
		print " Key not found \n";
		return undef; #object not found
	}

	return (@{$hash{$cle}}); #return the list/array
}
# private method : _isExist_hash
# function : say if a key hasg exist

sub _isExist_hash {
	my %hash = %{$_[0]};
	my $cle = $_[1];

	if (exists $hash{$cle}){
		return 1;
	}

	return 0;
}


# private method : _get_interactions_object
# function : look for all interaction having object $object

sub _get_interactions_object{
	my $self = shift;
	my $index = shift;


	return ( _getArray_from_hash(\%{$self->{'_objects'}},$index));
}

# private method : _get_interactions_subject
# function : look for all interaction having subject $subject (from a collection)

sub _get_interactions_subject{
	my $self = shift;
	my $index = shift;

	return ( _getArray_from_hash(\%{$self->{'_subjects'}},$index));
}

## NOTICE : all the methods get_accordingTO are a little obsolete, you should use the generic clever method
# get_interactions_with_tags_values (), with the tag and the value
#											 -----------------

# private method : _get_accordingTO_type
# function : get interaction of a specif type (intergenic or genic) from an ARRAY of interaction !!!
# not a collection

sub _get_accordingTO_type {
	my @tab = @{$_[0]}; # array  of Interactions object
	my $type = $_[1];
	my @array_type = ();
	for ( my $i = 0; $i <= $#tab; $i++){
		if( $tab[$i]->type() == $type){
			push(@array_type, $tab[$i]);
		}
	}
	return @array_type;
}

# private method : _get_accordingTO_type
# function : get interaction of a specif distance (integer) from an ARRAY of interaction !!!
# not a collection

sub _get_accordingTO_distance {
	my @tab = @{$_[0]}; # array  of Interactions object
	my $distance = $_[1];
	my @array_distance = ();
	for ( my $i = 0; $i <= $#tab; $i++){
		if( $tab[$i]->distance() == $distance){
			push(@array_distance, $tab[$i]);
		}
	}
	return @array_distance;
}

# get interaction of a specif direction (-1, 0, 1) from an ARRAY of interaction !!!

sub _get_accordingTO_direction {
	my @tab = @{$_[0]}; # array  of Interactions object
	my $direction = $_[1];
	my @array_direction = ();

	#print " __________________________________dans le get AC TO direction  ______ \n";
	for ( my $i = 0; $i <= $#tab; $i++){
		#print " dans le for $i \n";
		#print "dois je ajouter : ", $tab[$i]->direction() ," == ", $direction ,"\n";
		if( $tab[$i]->direction() == $direction){
			push(@array_direction, $tab[$i]);
		}
	}
	return @array_direction;
}

# get interaction of a spefic subtype AND type (0 (1,0), 1(1,2,3)) from an ARRAY of interaction !!!

sub _get_accordingTO_subtype_and_type {
	my @tab = @{$_[0]}; # array  of Interactions object
	my $type = $_[1];
	my $subtype = $_[2];
	my @array_subtype = ();
	#print " __________________________________dans le get AC TO subtype  ______ \n";
	for ( my $i = 0; $i <= $#tab; $i++){
		#print " dans le for $i \n";
		#print "dois je ajouter : ", $tab[$i]->type ," == ", $type ," && " , $tab[$i]->subtype() ," == ", $subtype," \n";

		if( $tab[$i]->type == $type && $tab[$i]->subtype() == $subtype){
		#	print "j ajoute car : ", $tab[$i]->type ," == ", $type ," && " , $tab[$i]->subtype() ," == ", $subtype," \n";
			push(@array_subtype, $tab[$i]);
		}
	}
	return @array_subtype;
}

# get interaction according to their nested statut (means that's you're looking for genic ones) from an ARRAY of interaction !!!

sub _get_accordingTO_nested {
	my @tab = @{$_[0]}; # array  of Interactions object
	my $nested = $_[1];
	my @array_nested = ();

	#	print " __________________________________dans le get AC TO nested  ______ \n";
	for ( my $i = 0; $i <= $#tab; $i++){
	#	print " dans le for $i \n";
	#	print "dois je regarder : ", $tab[$i]->type ," == 1\n";
		if( $tab[$i]->type == 1 && $tab[$i]->nested() == $nested){
	#		print "__ et l'autre attribut : ",$tab[$i]->nested() ," == ",$nested,"\n";
			push(@array_nested, $tab[$i]);
		}
	}
	return @array_nested;
}

# get interaction according to their divergent statut (means that's you're looking for intergenic ones) from an ARRAY of interaction !!!

sub _get_accordingTO_divergent {
	my @tab = @{$_[0]}; # array  of Interactions object
	my $divergent = $_[1];
	my @array_divergent = ();

	#	print " __________________________________dans le get AC TO divergent  ______ \n";
	for ( my $i = 0; $i <= $#tab; $i++){
	#	print " dans le for $i \n";
	#	print "dois je regarder : ", $tab[$i]->type ," == 0\n";
		if( $tab[$i]->type() == 0 && $tab[$i]->isDivergent() == $divergent){
	#			print "et l'autre attribut : ",$tab[$i]->isDivergent() ," == ",$divergent,"\n";
			push(@array_divergent, $tab[$i]);
		}
	}
	return @array_divergent;
}

=head1 get_all_interaction

Bio::SeqFeature::CollectionInteraction

=head1 DESCRIPTION

return an BIo::SeqFeature::Iterator object containing interactions

=cut
sub get_all_interactions {
	my $collection = shift;
 	my $iterator = Bio::SeqFeature::InteractionIterator->new();

 	foreach my $k ( keys (%{$collection->{'_objects'}}) ) {
 		my @array = _getArray_from_hash(\%{$collection->{'_objects'}}, $k);

 		$iterator->add(@array);
 	}
	return $iterator;
}

=head1 print_all_interaction

Bio::SeqFeature::CollectionInteraction

=head1 DESCRIPTION

print all the collection

=cut

sub print_all_interactions {
	my $collection = shift;
	my $biotype    = shift;
	my @objects=$collection->get_objects();

	#VW: print the header
	if($biotype)
	{
	    print "isBest\tlncRNA_gene\tlncRNA_transcript\tlncRNA_biotype\tpartnerRNA_gene\tpartnerRNA_transcript\tpartnerRNA_biotype\tdirection\ttype\tdistance\tsubtype\tlocation\n";
	}
	else
	{
	    print "isBest\tlncRNA_gene\tlncRNA_transcript\tpartnerRNA_gene\tpartnerRNA_transcript\tdirection\ttype\tdistance\tsubtype\tlocation\n";
	}

	foreach my $object (@objects) {
		my $best =$collection->get_the_best_interaction($object);
		$best->printer_mini('best',$biotype);
		foreach my $interaction ($collection->_get_interactions_object($object)) {
			unless ($interaction->subject() == $best->subject()) {
			    #VW: set the first arg as undef, i.e. not the best interaction
			    $interaction->printer_mini(undef,$biotype);
			}
		}
	}
}



=head1 get the_best_interaction

Bio::SeqFeature::CollectionInteraction

=head1 DESCRIPTION

my $interaction = $coll->get_the_best_interaction($feature);

return the best interaction for a particular Feature. This is to say :
	- if there is a genic interaction it returns those with the maximal overlap in bp, otherwise it returns the closest
=cut

sub get_the_best_interaction {
	my $collection = shift;
	my $object= shift;

	my $best;

	#print STDERR "obj : $object\n";
	foreach my $interaction  ($collection->_get_interactions_object($object)) {
	#	print STDERR "inteearction : $interaction\n";
			unless (defined $best) {
				$best= $interaction;
				next;
			}
			if ($interaction->isGenic()) {
				if ($best->isGenic) {
					if ($interaction->overlap() > $best->overlap()) {
						$best=$interaction;
					}
					if ($interaction->overlap() == $best->overlap() && $interaction->subject()->get_tag_values("transcript_id") < $best->subject()->get_tag_values("transcript_id")) {
						$best=$interaction;
					}
				}
				else {
					$best=$interaction;
				}
			}

			if ($interaction->isInterGenic) {
				if ($best->isGenic) {
					next;
				}
				else {
					if ($interaction->distance < $best->distance) {
						$best=$interaction;
					}
					if ($interaction->distance == $best->distance && $interaction->subject()->get_tag_values("transcript_id") < $best->subject()->get_tag_values("transcript_id")) {
						$best=$interaction;
					}
				}
			}
	}

	return $best;
}


=head1 sayme_type

Bio::SeqFeature::CollectionInteraction

=head1 DESCRIPTION

 say if the string type is correct
 return the associated number
 or the opposite
=cut
sub sayme_type{

my $type = shift;

	my %TYPES = ( 'intergenic' =>  '0' , 'genic' => '1' , '0' => 'intergenic' , '1' => 'genic');

	if (exists $TYPES{$type}){
		#print " type: ", $type , "\n";
		return ( $TYPES{$type} );
	}
	croak (" Unknow Type : \" $type \". Maybe the syntax is incorrect. \n Please read the manual to know the subtype we used \n") ; # error
}


=head1 get_interactions_with_tags OBSOLETE use the one with tags values

 return Bio::SeqFeature::InteractionIterator

=head1 DESCRIPTION

check if interactions corresponding to a list of fixed tags
better to use get_interactions_with_tags_values

=cut


sub get_interactions_with_tags{
	my $col = shift; #collection
	my $iterator = Bio::SeqFeature::InteractionIterator->new();

	if(@_){
		#print "les tags sont : \n";
		my %args = @_;

		#while (my ($cle, $v) = each(%args))  {
		#	print "cle : $cle , val : $v, _____ $args{$cle} -------\n";
		#	}


		# on wich interaction do we working ?

		### add objects interactions wich have $object as object

		if (exists $args{'-object'}) {
			if ( !defined ($args{'-object'})){
				carp(" Warning: Value object not defined \n ");
			}else{
				$iterator = $col->get_interactions($args{'-object'},1)  ;
				print " tag object set \n";
			}
		}elsif (exists $args{'-subject'}) {
			if ( !defined ($args{'-subject'})){
				carp(" Warning: Value subject not defined \n ");
			}else{
				$iterator = $col->get_interactions($args{'-subject'},2)  ;
				print " tag subject set \n";
			}
		}elsif ( exists $args{'-both'} ) {
			if ( !defined ($args{'-both'}) ){
				carp(" Warning: Value both not defined \n ");
			}else{
				$iterator = $col->get_interactions( $args{'-both'} ,3)  ;
				print " tag both set \n";
			}
		}else { # you want to work on the complete collection

			$iterator = $col->get_all_interactions();
			print " Default : work on all collection  \n";
		}

		if (_isDefined(@{$iterator->{'_array'}}) == 0){
			print " Number of interactions : 0 \n";
			_nointeraction_found();
			return $iterator;
		}
		# we have an iterator on an array of interactions


### ### add the interactions for tag subtype

		if ( _isExist_hash(\%args,'-subtype') ) {
			print " subtype  existe! \n";
			if ( !defined ($args{'-subtype'})){
				carp(" Warning: Value '-subtype' not defined \n ");
			}else{
				print " tag subtype set \n";
				my ($subtype,$type) = _sayme_subtype_type($args{'-subtype'});
				my @array = @{$iterator->{'_array'}};
				$iterator=Bio::SeqFeature::InteractionIterator->new( _get_accordingTO_subtype_and_type(\@array,$type,$subtype) ) ;
			}
		}else{
			print " no tag subtype \n";
		}



### ### add the interactions for tag direction

		if (exists $args{'-direction'}) {
			if ( !defined ($args{'-direction'})){
				carp(" Warning: Value '-direction' not defined \n ");
			}else{
				print " tag direction set \n";
				#print " \n looking for specific direction interaction  \n";
				my @array = @{$iterator->{'_array'}};
				$iterator=Bio::SeqFeature::InteractionIterator->new(_get_accordingTO_direction(\@array,$args{'-direction'}) ) ;

			}
		}else{
			print " no tag direction \n";
		}
### ### add the interactions for tag type

		if (exists $args{'-type'}) {
			if ( !defined ($args{'-type'})){
				carp(" Warning: Value '-type' not defined \n ");
			}else{
				print " tag type set \n";
				#print " \nlooking for specific type interaction  \n";
				my @array = @{$iterator->{'_array'}};
				my $type = sayme_type($args{'-type'});
				#print " le type est : ", $type, " ", $args{'-type'},"\n";
				$iterator=Bio::SeqFeature::InteractionIterator->new( _get_accordingTO_type(\@array, $type ))  ;
			}
		}else {
			print " no tag type \n";
		}
### ### add the interactions for tag divergent

		if (exists $args{'-divergent'}) {
			if ( !defined ($args{'-divergent'})){
				carp(" Warning: Value '-divergent' not defined \n ");
			}else{
				print " tag divergent set \n";
				#print " \nlooking for specific divergent interaction  \n";
				my @array = @{$iterator->{'_array'}};

				$iterator=Bio::SeqFeature::InteractionIterator->new( _get_accordingTO_divergent(\@array, $args{'-divergent'} ))  ;
			}
		}else{
			print " no tag divergent \n";
		}

### ### add the interactions for tag nested

		if (exists $args{'-nested'}) {
			if ( !defined ($args{'-nested'})){
				carp(" Warning: Value '-nested' not defined \n ");
			}else{
				print " tag nested set \n";
				#print " \nlooking for specific divergent interaction  \n";
				my @array = @{$iterator->{'_array'}};
				$iterator=Bio::SeqFeature::InteractionIterator->new( _get_accordingTO_nested(\@array, $args{'-nested'} ))  ;
			}
		}else{
			print " no tag nested \n";
		}

	}	#if of %args
	print " \n";
	return $iterator;
}


=head1 get_interactions_with_tags_values

 return Bio::SeqFeature::InteractionIterator

=head1 DESCRIPTION

whatever are the tags and values, check if interactions corresponding
see the list of tags in the manual , or using tags_list method

note:  if you pass an unexisting tag ...

=cut

sub get_interactions_with_tags_values{
	my $collection = shift; #collection
	my $iterator = $collection->get_all_interactions();

	if (_isDefined(@{$iterator->{'_array'}}) == 0){
			print " Number of interactions : 0 \n";
			_nointeraction_found();
			return $iterator;
	}
	my $iterator2 = Bio::SeqFeature::InteractionIterator->new();


	my %args = @_;
	while (my ($cle, $v) = each(%args))  {
		print "cle : $cle , val : $v \n";
		## inside the  interactions list
		while (my $interaction = $iterator->next()){
			 print " ----",$interaction->get_tag_value($cle) , " == ", $v , "\n";
			# if ($interaction->type()==0){$interaction->printer();}
			if ( $interaction->get_tag_value($cle) eq $v ){
				print  "jadd \n";
				$iterator2->add($interaction);
			}
		}
		#aucune gestion de la mémoire
		$iterator = $iterator2;
		$iterator2 =  Bio::SeqFeature::InteractionIterator->new();

	}

	return $iterator;
}

=cut
=head1 get_with_distance_value

 return Bio::SeqFeature::InteractionIterator

=head1 DESCRIPTION

You might want to have a list of  interaction wich have a specific distance.
Example : less or equal to 400 bases
	$collection->get_with_distance_value(400,'le')
	eq : equal
	lt : less than
	le : less or equal to
	gt : greater than
	ge : greater or equal to
=cut

sub get_with_distance_value{
	my $collection = shift;

	my $distance = shift;
	my $operator = shift;
	my $iterator = $collection->get_all_interactions();
	my $iterator_output = Bio::SeqFeature::InteractionIterator->new();



	if (_isDefined($operator) == 0 or $operator eq 'eq' ){
		$operator = "1";
	}elsif($operator eq 'le'){
		$operator = "2";
	}elsif($operator eq 'ge'){
		$operator = "3";
	}elsif($operator eq 'lt'){
		$operator = "4";
	}elsif($operator eq 'gt'){
		$operator = "5";
	}else{
		croak " get_with_distance_value: Error, undefined operator ( $operator ) \n"
	}

	while (my $interaction = $iterator->next()){ # je parcours mon itérateur dans le but de regarder si distance interaction correspond

		if($operator == 1 ){
			if( $interaction->distance() == $distance ) {
				$iterator_output->add($interaction);
			}
		}elsif($operator == 2 ){
			if( $interaction->distance() <= $distance ) {
				$iterator_output->add($interaction);
			}
		}elsif($operator == 3 ){
			if( $interaction->distance() >= $distance ) {
				$iterator_output->add($interaction);
			}
		}elsif($operator == 4 ){
			if( $interaction->distance() < $distance ) {
				$iterator_output->add($interaction);
			}
		}elsif($operator == 5 ){
			if( $interaction->distance() > $distance ) {
				$iterator_output->add($interaction);
			}
		}
	}
	return $iterator_output;
}

###
=cut
=head1 get_with_intervalle

 return Bio::SeqFeature::InteractionIterator

=head1 DESCRIPTION

You might want to have a list of interaction having distance within an intervalle

	$collection->get_interactions_with_intervalle(400,800)

=cut

sub get_interactions_with_intervalle {
	my $collection = shift;
	my $iterator = $collection->get_all_interactions();
	my $iterator_output = Bio::SeqFeature::InteractionIterator->new();
	my $inf = shift;
	my $sup = shift;

	my $rejet = 0;
	my $all = 0;
	if (_isDefined($inf) == 1 && _isDefined($sup) == 1){
		while (my $interaction = $iterator->next()){
			if( $interaction->distance() >= $inf && $interaction->distance() <= $sup ) {

				$iterator_output->add($interaction);

			}else{

				$rejet++;
			}
			$all++;
		}
	}

	print " Total number of interactions :  ", $all;

	print " Total number of rejected interactions :  ", $rejet;
	return $iterator_output;

}



=cut
=head1 objects_list

 return a list of Bio::SeqFeature

=head1 DESCRIPTION

	You might want to know the list of object containing in the collection

=cut

## hash keys 3rd level
sub objects_list {
	my $collection =shift;
	my @object = keys (%{$collection->{'_objects'}});
	my @objects;

		foreach my $k (@object){
			push (@objects, $collection->{'_objects'}->{$k}->[0]->object()); # then we have the object object

		}
	return @objects;
}
# private method : _sayme_subtype_type
# function : give me the subtype string I will give you the corresponding type integer and subtype integer
sub _sayme_subtype_type{

	my $subtype = shift;
	my @up = (1,0);
	my @down = (0,0);
	my @ex = (1,1);
	my @in = (2,1);
	my @ov = (3,1);
	my %SUBS_and_TYPES = ( 'upstream' =>  \@up , 'downstream' => \@down, 'exonic' => \@ex, 'intronic' => \@in , 'overlapping' => \@ov );


	if (exists $SUBS_and_TYPES{$subtype}){
		return ( @{$SUBS_and_TYPES{$subtype}} );
	}
	croak (" Unknow subtype : \" $subtype \". Maybe the syntax is incorrect. \n Please read the manual to know the subtype we used \n") ; # error
}
# private method : _nointeraction_found
# function : message where no interaction are founded
sub _nointeraction_found{
	print " Error : I am sorry. We did not detect any interaction. \n";
	print "         Please be sure the lncRNAs and mRNAs have tag type set to respectively lncRNAs and mRNAs \n";
	print "         Otherwise please contact authors \n";
}

# theses method has been copy paste from an other package :)
sub _isDefined{
	my $val = shift;
	if (defined $val) { return 1;} else { return 0;}
}

# control if the object and the subject are defined

sub _isComplete{
	my $self = shift;

	my $object = $self->object();
	my $subject = $self->subject();
	if ( _isDefined($object)  == 1 && _isDefined($subject) == 1){
		return 1;
	}else{
		return 0;
	}
}

1;
