package Bio::SeqFeature::InteractionIterator;

=head2 CLASS: Interator

 
Bio::SeqFeature::Interator

Use: Bio::SeqFeature::Interaction
=head1 ABOUT
Authors: Audrey DAVID - M2 Informatique opt Bioinfo Nantes
		 Fabrice LEGEAI - IRISA/INRIA Rennes / INRA Rheu
		 Thomas	 DERRIEN - CNRS Rennes
If you have questions please contact authors.

=head1 DESCRIPTION 

Implement several method to manipulate SImple list of Interactions
Implementation of classification rules of Derrien et al 2011. 
=head1 Others

Usually many of these method as Method has to be used like : 

my $res = $iterator->Method();

some of method can not return result than it's a printer and has to be use like
$iterator->Method();

Please check the description
=cut

use Bio::SeqFeature::Generic;
use Carp;
use strict;
use warnings;

=cut
=head1 new
	return an Bio::SeqFeature::InteractionIterator
=head1 DESCRIPTION 
	intialize an Bio::SeqFeature::InteractionIterator object or/and add interactions
=cut
sub new {
	my $pkg = shift;
  	my $iterator = bless {
	 _array => undef,
	 _index => -1,
	}, $pkg;


	if (@_) { # one array 
		#print "je te call : ____ $_[0] \n";
		
		$iterator->add(@_);# add the interactions array to the object
	}
		
  return $iterator;
}

=cut
=head1 add
=head1 DESCRIPTION 
	add a list of interaction 
=cut

sub add{
	my $iterator = shift;


	foreach my $interaction (@_ ){
		#print " ITERATOR : j'ajoute une nouvelle interaction : ", $interaction, "\n";
		
		if (_isDefined($interaction) == 1){
			if ( $interaction->isa('Bio::SeqFeature::Interaction')){
				push ( @{$iterator->{'_array'}} , $interaction );
			}else{
				croak (" Not an interaction object ! \n");
			}
		} #else : interaction is undef
	}
	
	$iterator->_init_index(); 
}
## initialize the index to 0 where the array is not empty
sub _init_index{
	my $iterator = shift;
	#print " ITERATOR : init l'index  s'il le faut \n";
	
	if ($#{$iterator->{'_array'}} >-1 && $iterator->_index() == -1 ) { # but the indice is at -1 (can be a turn back to the first element or the first declaration)
		$iterator->_index(0);
		#print " the indice was -1 \n ";
	}#else{ print " il ne le faut pas \n";}
	
}
## update according to the array size and the index
sub _update_index{
	
	my $iterator = shift;
	#print " ITERATOR : j update l'index \n";
	
	if ($#{$iterator->{'_array'}} >-1){ # the array is not empty
		#print " the array is not empty \n";
		if ($iterator->_index() == -1 ) { # but the indice is at -1 (can be a turn back to the first element or the first declaration)
			$iterator->_index(0);
			#print " the indice was -1 \n ";
		}else {
			$iterator->_index($iterator->_index() + 1); # next object
			#print " prochain indice : ", $iterator->_index(), "\n";
			#sleep(2);
		}	
	}
}
## return and/or set the index 
sub _index {
	my $self = shift;
	if (@_){
		$self->_set_index($_[0]);
	}

	return $self->{'_index'};
}
## set the index 
sub _set_index{
	my $self = shift;
	my $index = shift;
	#print (" l indice suivant doit être ", $index,"\n");
	#print (" la taille est: ", $#{$self->{'_array'}} ,"\n");

	if ($index > $#{$self->{'_array'}}) {
		$self->{'_index'} = -1; # back to undef state
	#	print "je te set à -1 \n";
	#	sleep(2);
	}else{
		$self->{'_index'} = $index;
	}
}

=cut
=head1 iterator first
	return Bio::SeqFeature::Interaction
=head1 DESCRIPTION 
	Next interaction 
=cut

sub next {
	my $iterator = shift;
	my $i = $iterator->_index();
	if (  $i == -1 ) {
		$iterator->_update_index();
		return undef; # no more interaction
	}else{
		$iterator->_update_index();
		return $iterator->{'_array'}[$i];
	}	
}


=cut
=head1 iterator first
	return Bio::SeqFeature::InteractionIterator
=head1 DESCRIPTION 
	Reinitialize the iterator (back to the beginning) 
=cut

sub first {
	my $iterator = shift;
	$iterator->_index(-1);	
	$iterator->_update_index();
	return $iterator;
}


=cut
=head1 iterator copy
	return Bio::SeqFeature::InteractionIterator
=head1 DESCRIPTION 
	 Copy an iterator.
	 Use it before use fusion or/and intersection method, if you want to save your iterator.
=cut

sub copy{
	my $iterator = shift;
	my $copy = Bio::SeqFeature::InteractionIterator->new();
	while (my $interaction = $iterator->next()){
		$copy->add($interaction);
	}
	return $copy;
}


=cut
=head1 iterator fusion
	return Bio::SeqFeature::InteractionIterator
=head1 DESCRIPTION 
	Union operation on two iterator
	$it = $it->fusion($it2)
=cut


sub fusion{

	my $iterator = shift;
	my $fusion = shift;
	my %hash = map { $_ => 1 } @{$iterator->{'_array'}};

	while (my $interaction = $fusion->next()){
		if (!exists $hash{$interaction}){
			$iterator->add($interaction);
		}
	}
	return $iterator;
}

=cut
=head1 iterator intersection
	return Bio::SeqFeature::InteractionIterator
=head1 DESCRIPTION 
	Intersection operation on two iterator
	$it = $it->intersection($it2)	
=cut

sub intersection{
	my $self = shift;
	# may be you can add the size information to accelarate a litte bit the treatment
	my %hash = map { $_ => 1 } @{$self->{'_array'}};
	my $iterator = shift;
	my $iterator_output =  Bio::SeqFeature::InteractionIterator->new();

	while (my $interaction = $iterator->next() ) {
		if (exists $hash{$interaction}){
			$iterator_output->add($interaction);
		}
	}
	return $iterator_output;
}
sub remove_redondancy{}

############### méthode défini dans  package Bio::SeqFeature::Interaction->meth

sub _isDefined{
	my $val = shift;
	if (defined $val) { return 1;} else { return 0;}
}
1;