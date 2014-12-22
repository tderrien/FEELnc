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
		_nested => undef,
		_subtype => undef,
		@_,
	);
	#TODO n
	#is it really genic? 
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
	my %subtypes =( 1 => 'exonic', 2=>' intronic', 3 => 'overlapping');
	if (defined($subtypes{$sub_t})){
		$self->{'_subtype'} = $sub_t;

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
=head1 isOverlapping() 

Bio::SeqFeature::Interaction

=head1 DESCRIPTION 

if it's not intronic and not exonic it's overlapping
see classification rules (T.Derrien)

=cut

sub isOverlapping{
	
	my $interaction = shift;
	if ($interaction->_isComplete() == 0) {
		
		croak ("To look if your interaction is Overlapping, you must have an complete interaction \n");
	}	
	if(_isDefined($interaction->{'_subtype'})==0) {
		$interaction->_set_subtype($interaction->_calculate_subtype()); # what is the subtype?
	}
	if($interaction->_get_subtype()==3){
		return 1;
	}else {return 0;}
	
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
 
			return 1;
		}elsif ( isIntronic($interaction)){
			return 2;
		}else{
			return 3; #overlapping
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
		
		croak ("To look if your interaction is Exonic, you must have an complete interaction \n");
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
	my $object = $interaction->object();
	my $subject = $interaction->subject();

	my @introns_mRNA= _getIntrons($object);
	my @exons_lncRNA = _getExons($subject);

	#print "intronic call \n";
	
	if (_overlaps_btwn_feature_array(\@exons_lncRNA, \@introns_mRNA, 1) == 1){
		return 1;
	}else{
		return 0;
	}
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
			#print "exons", $array1[$i]->start, " - ", $array1[$i]->end() , "\n"; #partie censé se chevaucher
			#print " \t\t introns ( ",$array2[$j]->start, " - " , $array2[$j]->end(), ")\n"; 
		
			if ($array2[$j]->overlaps($array1[$i])){
				$required++;
				#print "oui \n";
				if ($required == $number ){ #do we obtain the number of overlaps required?
					return 1;
				}
			}
		}
	}
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
	my @f_array = shuffle @fs_array;
	my $i = 0;
	my %hash=();
	#print "taille array : ", $#f_array, "\n";
	for ($i =0 ; $i<=$#f_array; $i++){
		$hash{$i} = $f_array[$i]->start();
		
	}
 
	my @keys = sort { $hash{$a} <=> $hash{$b} } keys %hash; # sort the hash table on values 
 	$i = 0;
	foreach my $cle (@keys) {		
		$fs_array[$i] = $f_array[$cle]; # cle : indice of the array 
		
		$i++;
	}
	return (@fs_array);
}
				# ---------------------------------------------------- #
				# 					printer							   #

				# ---------------------------------------------------- #
=head1 printer() 

=head1 DESCRIPTION 

$interaction->printer()
print interaction informations

=cut	
# this fonction use the SUPER printer and is also defined in intergenic
# sub printer {
# 	my $self = shift;
# 	$self->SUPER::printer();
# 	print "\t Other info \n";
# 	Bio::SeqFeature::Genic->_ligne_carre();
# 	print "\t Nested : ",$self->nested(),"\n";
# 	print "\t Subtype : ",_print_subtype($self->subtype()),"\n";
# 
# }


=head1 printer_mini() 
=head1 DESCRIPTION 

$interaction->printer_mini()
print short interaction informations

=cut	
# this fonction use the SUPER printer and is also defined in intergenic


sub printer_mini {
	my $self = shift;
	$self->SUPER::printer_mini();
#	print "\t Other info \n";
#	Bio::SeqFeature::Genic->_ligne_carre();
	print "\t Status=",$self->nested();
	print "\t Subtype=",_print_subtype($self->subtype()),"\n";

}

=head1 DESCRIPTION 

$interaction->printer()
print interaction informations

=cut	
# this fonction use the SUPER printer and is also defined in intergenic
sub printer {
	my $self = shift;
	$self->SUPER::printer();
	print "\t Other info \n";
	Bio::SeqFeature::Genic->_ligne_carre();
	print "\t Nested : ",$self->nested(),"\n";
	print "\t Subtype : ",_print_subtype($self->subtype()),"\n";

}
				# ---------------------------------------------------- #
				# ---------------------------------------------------- #
#_print_subtype : private method
# function : do the correspondance between the subtype number and the string corresponding

sub _print_subtype{
	my $subtype = shift;
	my %subtypes =( 1 => 'exonic', 2=>' intronic', 3 => 'overlapping');
	if (exists $subtypes{$subtype}) {
		return $subtypes{$subtype};
	}
	return $subtype;
}
				# ---------------------------------------------------- #
				# ---------------------------------------------------- #
=head1 nested() 

Bio::SeqFeature::Generic

=head1 DESCRIPTION 

my $nested = $interaction->nested()

say if the interaction is nested or not (see classification rules Derrien et al 2011)
=cut	

sub nested {
	my $self =shift;
	my $nested = $self->_get_nested();
	return $nested;
}
				# ---------------------------------------------------- #
				# ---------------------------------------------------- #
#_get_nested : private method
# function : getter of nested attribute
sub _get_nested{

	my $self = shift;
	if ( _isDefined($self->{"_nested"}) == 0 ){
		$self->_set_nested($self->_isNested());
	}
	return $self->{'_nested'};
}
				# ---------------------------------------------------- #
				# ---------------------------------------------------- #

#_set_nested : private method
# function : setter of nested attribute
sub _set_nested{
	my $self = shift;
	if (@_) {
		$self->{'_nested'} = shift;
	}else {
		croak ("Nested undefined, cannot set nested \n"); # i want it to die because it's a protected method; then the developper make an error 
	}
}
	
				# ---------------------------------------------------- #
				# ---------------------------------------------------- #
# isNested : private method
# function : say if nested or not

sub _isNested{
	
	my $self = shift;
	
	if ($self->_isComplete() == 0) {	
		croak ("To look if your interaction is Nested, you must have an complete interaction \n");
	}	
	my $object = $self->object();
	my $subject = $self->subject();
 
	if ($object->contains($subject)){
	 
		return 1;
	}	else{
	 
		 return 0;
	}
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
=head1 gestion_error() 
IMPORTANT

Bio::SeqFeature::Genic

=head1 DESCRIPTION 
when a genic interaction is created, this function has to be called to check if it's an exonic same strand interaction
then a error has to be solved
=cut
sub gestion_error{ # if exonic and direction == 1 that's an error
	
	my $self = shift;
	if ($self->isExonic == 1){
		if ($self->direction == 1){
			print "The subject and the object are in the same direction and there is at least one overlap btwn exon \n";
			print "If you are trying to class lncRNAs in relation to mRNA, that's an error of classification \n";
			print "Please have a look of your entry sets \n";
			print "If you use our method ... please contact us \n";
			sleep(1);
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