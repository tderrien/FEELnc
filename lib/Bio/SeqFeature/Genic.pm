package Bio::SeqFeature::Genic;

=head2 CLASS: Genic

SUPER CLASS : Interaction

Bio::SeqFeature::Genic
=head1 ABOUT
Authors: Audrey DAVID - M2 Informatique opt Bioinfo Nantes
		 Fabrice LEGEAI - IRISA/INRIA Rennes / INRA Rheu
		 Thomas	 DERRIEN - CNRS Rennes
If you have questions contact authors.

=head1 DESCRIPTION

Implement several method to manipulate genic interaction
Implementation of classification rules of Derrien et al 2011.

=head1 Others

Usually many of these method as Method has to be used like :

my $res = $interaction->Method();

some of method can not return result than it's a printer and has to be use like
$interaction->Method();

Please check the description
=cut
@ISA = qw(Bio::SeqFeature::Interaction);
use Carp;
use strict;
use List::Util qw(shuffle);
use Bio::SeqFeature::Generic;



=head1 new

Bio::SeqFeature::Genic
=head1 DESCRIPTION
Constructor

=cut
sub new {

	my $pkg = shift;
	$pkg->SUPER::new(
		_type => 1,
		_distance => 0,
		_status => undef,
		_subtype => undef,
		@_,
	);

}


=head1 subtype()

Bio::SeqFeature::Genic
=head1 DESCRIPTION

$interaction->subtype();
get subtype information

=cut


sub subtype{

	my $self = shift;
	return  $self->_get_subtype() ;
}

# _set_subtype : private method
# function : set the subtype ( number => string)
sub _set_subtype(){
	my $self = shift;
	my $sub_t = shift;
	my %subtypes =('exonic' => 1 ,'intronic' => 1, 'none'=> 1);
	if (defined($subtypes{$sub_t})){
		$self->{'_subtype'} = $sub_t;
	}
	else {
		croak("Error unrecognized subtype : $sub_t\n");
	}
}

#_get_subtype() : private method
# function : getter (get the subtype); called by subtype()

sub _get_subtype{
	my $self = shift;
	if ( _isDefined ($self->{'_subtype'}) == 0) { #undef
		$self->_set_subtype($self->_calculate_subtype());
	}
	return $self->{'_subtype'};
}


=head1 calculate_subtype()

Bio::SeqFeature::Interaction

=head1 DESCRIPTION

if it's not intronic and not exonic it's overlapping
see classification rules (T.Derrien)

=cut

sub _calculate_subtype(){
	my $interaction = shift;
	if ($interaction->_isComplete() == 0) {
		croak ("To make operation on interaction, you must have an complete interaction \n");
	}

	if (isExonic($interaction)){
			return 'exonic';
		}elsif ( isIntronic($interaction)){
			return 'intronic';
		}else{
			return 'none'; #overlapping
		}
}

=head1 isExonic()

Bio::SeqFeature::Interaction

=head1 DESCRIPTION

say if there is/are at least $_[3] (number) intersection btwn exon of the subject seqfeature and exons of the object seqfeature
see classification rules from T.Derrien

=cut

sub isExonic{

	my $interaction = shift;
	if ($interaction->_isComplete() == 0) {

		croak ("To look if your interaction is Exonic, you must have a complete interaction \n");
	}
	my $object = $interaction->object();
	my $subject = $interaction->subject();

	my @exons_mRNA= _getExons($object);
	my @exons_lncRNA = _getExons($subject);
	if (_overlaps_btwn_feature_array(\@exons_lncRNA, \@exons_mRNA, 1) == 1){
		return 1;
	}else{
		return 0;
	}
}


=head1 isIntronic()

Bio::SeqFeature::Interaction

=head1 DESCRIPTION

say if one exon of the subject seqfeature overlaps $_[3] (number) intron of the object seqfeature
see classification rules from T.Derrien

=cut
sub isIntronic {

	#my $pkg = shift;
	my $interaction = shift;
	if ($interaction->_isComplete() == 0) {
		croak ("To look if your interaction is Intronic, you must have an complete interaction \n");
	}

	my $lncrna = $interaction->object();
	#print STDERR $lncrna->get_tag_values("transcript_id"), " start: ", $lncrna->start(), " end: ", $lncrna->end(), "\n";

	my @exons_lncRNA = _getExons($lncrna);

	my $mrna = $interaction->subject();

	my @introns_mRNA= _getIntrons($mrna);


	#print "intronic call \n";
	# VW: put the line in commentary because it causes a case where monoexonic mRNA in intron of lncRNA are classed as none instead of containing->intron
	#if (scalar(@introns_mRNA) == 0) {return 0}


	foreach my $lncrna_exon (@exons_lncRNA) {
		next if ($lncrna_exon->end() < $mrna->start());
		next if ($lncrna_exon->start() > $mrna->end());
		unless (_is_included($lncrna_exon, \@introns_mRNA)) {return 0}
	}

	return 1;
}


# _is_included : private method
# it returns true if the first argument is included into at least one of the following arguments

sub _is_included {
	my $first = shift;
	my $others = shift;

	foreach my $other (@{$others}) {
		if ($other->contains($first)) {return(1)}
	}

	return (0);
}


#_overlaps_btwn_feature_array : private method
# function : say if there is an overlap between the exons of 2 features
					# ---------------------------------------------------- -------------------------#
					#_overlaps_btwn_feature_array(@exons_lncRNA, @introns_mRNA, $number_of_overlaps)#
					# ---------------------------------------------------- -------------------------#
sub _overlaps_btwn_feature_array{

	my @array1 = @{$_[0]};
	my @array2 = @{$_[1]};
	my $number = $_[2];
	my $required = 0;
	my $i =0 ;
	my $j=0;
	#print "overlaps \n";
	for ( $i = 0; $i <= $#array1; $i++) {


		for ( $j =0 ; $j <= $#array2 ; $j++){
			if ($array2[$j]->overlaps($array1[$i])){
				$required++;
				if ($required == $number ){ #do we obtain the number of overlaps required?
					return 1;
				}
			}
		}
	}
}


=head1 overlap

=head1 DESCRIPTION

returns the number of bases that are overlapping between the subject and the object

=cut

sub overlap {
	my $interaction=shift;

	my $object = $interaction->object();
	my $subject = $interaction->subject();

	my @exons_mRNA= _getExons($object);
	my @exons_lncRNA = _getExons($subject);

	my $overlap=0;
	foreach my $exon_mRNA (@exons_mRNA) {

		foreach my $exon_lncRNA (@exons_lncRNA) {
				my $first=$exon_mRNA;
				my $second=$exon_lncRNA;
				if ($exon_mRNA->start() > $exon_lncRNA->start()) {
					$first=$exon_lncRNA;
					$second=$exon_mRNA;
				}


				next if ($first->end() < $second->start);

				my $begin = $first->start();
				my $end = $first->end();
				if ($second->end() < $first->end) {
					$end = $second->end();
				}

				$overlap= $overlap+$end-$begin;
		}
	}
	return $overlap;
}

=head1 _getIntrons()

protected
Bio::SeqFeature::Generic

=head1 DESCRIPTION

get the position of intron from an array of ordinate seqfeature exons

return an array of seqfeatures introns with coordinates e.g. start & end

=cut


sub _getIntrons{

	my $feature = shift;
	my @exons= _getExons($feature);
	my @introns=();
	for (my $i=0; $i < $#exons ; $i++){

		$introns[$i] = Bio::SeqFeature::Generic->new( -start => $exons[$i]->end+1 ,
													-end => $exons[$i+1]->start-1 ,
													-strand => $exons[$i]->strand(),
												    );
		$introns[$i]->seq_id($exons[$i]->seq_id());
		#print "  ", $exons[$i]->end(), " (", $introns[$i]->start(), " - ",$introns[$i]->end," ) ", $exons[$i+1]->start(), "\n";


	}
	return @introns;

}

=head1 _getExons()

protected
Bio::SeqFeature::Generic

=head1 DESCRIPTION

get the subfeature exon from a seqfeature transcripts
return an array of ordinate exons coordinates

=cut

sub _getExons {

	my $feature = shift;
	my @exons = $feature->get_SeqFeatures();
	return (_order_array_features(@exons));
}

#_order_array_features : private method
# function : for example : return an ordonnate exons array
				# ---------------------------------------------------- #
				# order an array of feature by position (start and ends)
				# ---------------------------------------------------- #

sub _order_array_features {
	my @fs_array = @_; #features array
	my @ordered_array = sort {$a->start() <=> $b->start()} @fs_array;


	return (@ordered_array);
}



=head1 printer_mini()
=head1 DESCRIPTION

$interaction->printer_mini()
print short interaction informations

=cut
# this fonction use the SUPER printer and is also defined in intergenic


sub printer_mini {
	my $self = shift;
	my $best=shift;
	my $biotype = shift;

	#VW: for the biotype printing
	$self->SUPER::printer_mini($best, $biotype);

	# print "\t Status=",$self->status();
	# print "\t Subtype=", $self->subtype(),"\n";
	#VW modif
	print "\t", $self->status();
	print "\t", $self->subtype(),"\n";

}

=head1 status()

Bio::SeqFeature::Generic

=head1 DESCRIPTION

my $status = $interaction->status()

says what is the interaction status (containing, nested or overlapping) or not (see classification rules Derrien et al 2011)
=cut

sub status {
	my $self =shift;
	my $status = $self->_get_status();
	return $status;
}


#_get_status : private method
# function : getter of status attribute

sub _get_status{
	my $self = shift;
	if ( _isDefined($self->{"_status"}) == 0 ){
		$self->_set_status($self->_calculate_status());
	}
	return $self->{'_status'};
}


#_set_status : private method
# function : setter of status attribute


sub _set_status{
	my $self = shift;
	if (@_) {
		my $status = shift;
		my %status =('nested' => 1 ,'containing' => 1, 'overlapping'=> 1);
		unless (exists $status{$status}) {
			croak("Error unrecognized status : $status\n");
		}
		$self->{'_status'} = $status;
	}
	else {
		croak ("type undefined, cannot set \n"); # i want it to die because it's a protected method; then the developper make an error
	}
}

				# ---------------------------------------------------- #
				# ---------------------------------------------------- #

# calculate_status : private method
# function : returns the status of the interaction (nested, containing or overlapping)

sub _calculate_status{
	my $self = shift;

	if ($self->_isComplete() == 0) {
		croak ("To get the status of your interaction, you must have a complete interaction \n");
	}

	my $lncrna = $self->object();
	my $mrna = $self->subject();

 	if ($lncrna->contains($mrna)) {
 		return 'containing';
 	}


	if ($mrna->contains($lncrna)){
		return ('nested');
	}

	 return 'overlapping';
}


=head1 reciprocal()
DO NOT USED IT prior to create a new interaction in wich object and subject are switch
Bio::SeqFeature::Genic

=head1 DESCRIPTION

return the reciprocal interaction e.g. subject and object are reversed

=cut
sub reciprocal {

	my $self = shift;
	my %reciprocal_hash=(3 =>2, 1 => 1, 2=>2);
	my %nested_reciprocal = (1=> 0, 0=>1);


	my $reciprocal_object = Bio::SeqFeature::Genic->new( "-object" => $self->subject(), "-subject"=> $self->object(), "-type" => $self->type(),
								"-distance"=>$self->distance(), "_nested"=> $nested_reciprocal{$self->nested()} , "_subtype"=> $reciprocal_hash{$self->subtype()}
							);
}

				# ---------------------------------------------------- #
				# ---------------------------------------------------- #
=head1 warning()

Bio::SeqFeature::Genic

=head1 DESCRIPTION
when a genic interaction is created, this function has to be called to check if it's an exonic same strand interaction
then a error has to be solved
=cut
sub warning{ # if exonic and direction == 1 that's an error

	my $self = shift;
	if ($self->isExonic == 1){
		if ($self->direction == 1){
			print STDERR  $self->object()->get_tag_values("transcript_id"), " and ",  $self->subject()->get_tag_values("transcript_id") , " are overlapping in the same strand\n";
			return 0; # fatal error for the entry
		}
	}
	return 1;
}

# theses method has been copy paste from an other package :)

# control if a value is defined
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
		return 0; #not complete
	}
}

1;
