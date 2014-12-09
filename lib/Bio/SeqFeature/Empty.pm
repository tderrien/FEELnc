package Bio::SeqFeature::Empty;

@ISA = qw(Bio::SeqFeature::Interaction);
use Carp;
use strict;
use List::Util qw(shuffle);
use Bio::SeqFeature::Generic;
sub new {
	
	my $pkg = shift;
	$pkg->SUPER::new(
		_type => 555,
		_empty => 1,
		@_,
	);
}

				# ---------------------------------------------------- #
				# 					printer							   #
				# ---------------------------------------------------- #
	
sub printer {
	my $self = shift;
	$self->SUPER::printer();
	print "\t Other info \n";
	Bio::SeqFeature::Genic->_ligne_carre();
	print "\t Empty interaction\n";

}

sub type {
	return 555;
}

1;