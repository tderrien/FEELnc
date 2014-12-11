package Bio::SeqFeature::Extended;

use Bio::SeqFeature::Generic;
@ISA=('Bio::SeqFeature::Generic');


=head1 NAME 

Bio::SeqFeature::Extended 

=head1 DESCRIPTION 

bla bla 

=cut


=head2 add_interaction
 Title   : add_interaction
 Usage   : $feature->add_interaction($interaction);
 Function: 
 Returns : the start of this range
 Args    : optionaly allows the start to be set
          : using $loc->start($start)
=cut


sub add_interaction {
	$self = shift;
	print STDERR "troto\n";
	my $interaction = shift;
	unless ($interaction->isa('Bio::SeqFeature::Interaction')) {
		carp("Error add_interaction needs a Bio::SeqFeature::Interaction\n");
		exit();
	}
	push @{$self->{'_interactions'}}, $_;	
}


=head2 get_interactions
 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 
 
=cut

sub get_interactions {
	$self = shift;
	return @{$self->{'_interactions'}};	
}	