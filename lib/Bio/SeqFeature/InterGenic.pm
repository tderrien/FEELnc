package Bio::SeqFeature::InterGenic;

=head2 CLASS: InterGenic

SUPER CLASS : Interaction

Bio::SeqFeature::InteractionCollection

Use : Interaction

=head1 ABOUT
Authors: Audrey DAVID - M2 Informatique opt Bioinfo Nantes
		 Fabrice LEGEAI - IRISA/INRIA Rennes / INRA Rheu
		 Thomas	 DERRIEN - CNRS Rennes
If you have questions contact authors.

=head1 DESCRIPTION

Implement several method to manipulate intergenic interaction
Implementation of classification rules of Derrien et al 2011.

=cut

@ISA = qw(Bio::SeqFeature::Interaction);
use Carp;
use strict;
use List::Util qw(shuffle);
use Bio::SeqFeature::Generic;

=cut
=head1 new()

Bio::SeqFeature::InterGenic

=head1 DESCRIPTION

constructor
=cut
sub new {

	my $pkg = shift;
	#print "je supoer \n";
	$pkg->SUPER::new(
		_type => 0,
		_divergent => undef,
		_subtype => undef,
		@_
	);

}

# private method : _get_upstream
# function : get upstream parameter
# arg(s): $self
# return: upstream value

sub _get_upstream{
	my $self = shift;
	if ( _isDefined($self->{"_upstream"}) == 0 ){
		$self->_set_upstream($self->isUpstream());
	}
	return $self->{'_upstream'};

}
# private method : _get_divergent
# function : get divergent parameter
# arg(s): $self
# return: divergent value

sub _get_divergent{
	my $self = shift;
	if ( _isDefined($self->{"_divergent"}) == 0 ){
		$self->_set_divergent($self->isDivergent());
	}
	return $self->{'_divergent'};
}


# private method : _set_upstream
# function : set upstream parameter
# arg(s): $self
# return: /
sub _set_upstream{
	my $self = shift;
	if (@_){
		$self->{"_upstream"} = shift;
	}else{
		croak " Upstream value undefined \n";
	}
}

# private method : _set_divergent
# function : set divergent parameter
# arg(s): $self
# return: /
sub _set_divergent{
	my $self = shift;
	if (@_){
		$self->{"_divergent"} = shift;
	}else{
		croak " Divergent value undefined \n";
	}
}

=cut
=head1 isUpstream()

=head1 DESCRIPTION
	my $res = $interaction->isUpstream();
	return 1 if it's an upstream interaction
	return 0 if it's a downstream interaction
=cut

sub isUpstream{
	my $self = shift;
	if( _isComplete($self) == 0){
		carp "Incomplete interaction : Cannot look for an upstream interaction.\n ";
		return 555;
	}
	my $object = $self->object();
	my $subject = $self->subject();

	# Obj = lncRNA
	#print STDERR "obj : ", $object->get_tag_values("transcript_id"), "\n";

	if ($subject->strand == 1) {
		if ($object->start() > $subject->end()){
			return 0; #is downstream
		}
		return 1; #is upstream
	}

	if ($subject->strand == -1) {
		if ($object->start() > $subject->end()){
			return 1; #is upstream
		}
		return 0; #is  downstream
	}

}

# private method : _get_subtype
# function : get subtype parameter
# genic or intergenic
# arg(s): $self
# return: subtype value

sub _get_subtype{
	my $self = shift;
	if ( _isDefined ($self->{'_subtype'}) == 0) { #undef

		$self->_set_subtype($self->_calculate_subtype());
	}
	return $self->{'_subtype'};
}

# private method : _calculate_subtype
# function : calculate subtype parameter
# genic or intergenic
# arg(s): $self
# return: 1 or 0

sub _calculate_subtype(){
	my $interaction = shift;
	if ($interaction->_isComplete() == 0) {
		croak ("To make operation on interaction, you must have an complete interaction \n");
	}
	if ($interaction->isUpstream){
		return 1; # upstream
	}
	return 0; # downstream

}


# private method : _set_subtype
# function : set subtype parameter
# genic or intergenic
# arg(s): $self
# return: /
sub _set_subtype{
	my $self = shift;
	if (@_){
		$self->{"_subtype"} = shift;
	}else{
		croak " Subtype value undefined \n";
	}
}


=cut
=head1 isDownstream()

=head1 DESCRIPTION
	my $res = $interaction->isDownstream();
	return 1 if it's a downstream interaction
	return 0 if it's an upstream interaction
=cut

# oppposite :  is amont
sub isDownstream{
	my $self = shift;
	my $resp = $self->isUpstream();
	if ($resp == 1){
		return 0;
	}else{
		return 1;
	}
}


=head1 isDivergent()

Bio::SeqFeature:Interaction

=head1 DESCRIPTION
	my $res = $interaction->isDivergent();
	say if the interaction belong to Divergent class
	object direction 1 , subject direction -1
=cut

sub isDivergent{
	my $self = shift;
	my $direction = $self->direction();
	if ($direction == 1){ # same strand
		return -1;
	}elsif ($direction == 0){ #unknow direction
				return -2;
	}else{ # anti direction
		my $object = $self->object();
		my $subject = $self->subject();

		#		if ( ($object->strand() ==1 &&  $self->isUpstream()) || ($object->strand() == -1 && $self->isDownstream() == 1 ) ) {
		# VW: modification of the condition
		if ( $self->isDownstream() == 1 ) {
			return 0;	# convergent ; we're already known that it's in anti sens ( [1;-1] or [-1;1])
		}
		return 1;
	}
}


=head1 subtype()

Bio::SeqFeature:Interaction

=head1 DESCRIPTION
	my $res = $interaction->subtype();
	get subtype value of an Bio::SeqFeature:Interaction object
=cut

sub subtype{

	my $self = shift;
	return $self->_get_subtype();
}


=head1 printer()

Bio::SeqFeature:Interaction

=head1 DESCRIPTION
	$interaction->printer();
	print information of an Bio::SeqFeature:Interaction object

=cut
# this fonction use the SUPER printer and is also defined in intergenic
sub printer {
	my $self = shift;
	$self->SUPER::printer();
	print "\t Other info \n";
	Bio::SeqFeature::Genic->_ligne_carre();
	print "\t About divergent statut: ",_print_divergent($self->_get_divergent()),"\n";
	print "\t Subtype : ",_print_subtype($self->subtype()),"\n";

}


=head1 printer_mini()

Bio::SeqFeature:Interaction

=head1 DESCRIPTION
	$interaction->printer_mini();
	print information of an Bio::SeqFeature:Interaction object

=cut
# this fonction use the SUPER printer and is also defined in intergenic
sub printer_mini {
	my $self = shift;
	my $best=shift;
	my $biotype = shift;

	#VW: for the biotype printing
	$self->SUPER::printer_mini($best,$biotype);
	# print "\t Status=",_print_divergent($self->_get_divergent());
	# print "\t Subtype=",_print_subtype($self->subtype()),"\n";
	#VW modif
	print "\t",_print_divergent($self->_get_divergent());
	print "\t",_print_subtype($self->subtype()),"\n";

}
# private method : _print_divergent
# function : print divergent parameter
# arg(s) : $self
# return : /
sub  _print_divergent{
	my $div = shift;
	my %div_types = ( -2 => "unknow strand(s)", -1 => 'same_strand', 1 => 'divergent', 0 => 'convergent');
	if (exists $div_types{$div}) {
		return $div_types{$div};
	}
	return $div;
}

# private method : _print_subtype
# function : print divergent parameter
# arg(s) : $self
# return : /

sub _print_subtype{
	my $subtype = shift;
	my %subtypes = ( 1 => 'upstream', 0 => 'downstream');
	if (exists $subtypes{$subtype}) {
		return $subtypes{$subtype};
	}
	return $subtype;
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
