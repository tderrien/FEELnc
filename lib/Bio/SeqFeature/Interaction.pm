package Bio::SeqFeature::Interaction;

=head2 CLASS: Interaction

SUPER CLASS : NONE
Bio::SeqFeature::Interaction

Use Bio::SeqFeature
=head1 ABOUT
Authors: Audrey DAVID - M2 Informatique opt Bioinfo Nantes
		 Fabrice LEGEAI - IRISA/INRIA Rennes / INRA Rheu
		 Thomas	 DERRIEN - CNRS Rennes
If you have questions please contact authors.

=head1 DESCRIPTION

Implement several method to manipulate Interaction
Implementation of classification rules of Derrien et al 2011.
=head1 Others

Usually many of these method as Method has to be used like :

my $res = $interaction->Method();

some of method can not return result than it's a printer and has to be use like
$interaction->Method();

Please check the description
=cut

use strict;
use Carp;

=head1 new

Bio::SeqFeature::Interaction

=head1 DESCRIPTION

	_subject: Bio::SeqFeatureI
	_object : Bio::SeqFeatureI
	_overlaps: Boolean
	_direction 	 : 1, 0 or -1
	_distance: integer value
notice :
	only the object and the subject can be set using new method
	the others attributes have to be calculate using the appropriate method
=cut


our $VERSION='1.01';
our $ABSTRACT='';

sub new {
my ($pkg, %args) = @_;

my $interaction;
	if( ($pkg eq 'Bio::SeqFeature::InterGenic') or ($pkg eq 'Bio::SeqFeature::Genic')) {

	  	$interaction = bless {
		 _subject => undef,
		 _object => undef,
		 _overlap => undef,
		 _direction => undef,
		 _distance => undef,
		}, $pkg;
	}elsif ($pkg eq 'Bio::SeqFeature::Empty') {
		$interaction = bless {
		 _subject => undef,
		 _object => undef,
		 _empty => 1,
		}, $pkg;
	}else {
		croak (" Abstract class : Error \n Report to the manual or contact the authors \n");
	 	## make it looks like an abstract class
	}


	if (exists $args{'-subject'}) {
		if ( !defined ($args{'-subject'})){
			carp(" Warning: Value subject not defined \n ");
		}else{
			$interaction->subject($args{'-subject'});
		}
	}
	if (exists $args{'-object'}) {
		if (!defined ($args{'-object'}) ){
			carp(" Warning: Value object not defined \n ");
		}else{
			$interaction->object($args{'-object'});
		}
	}

  return $interaction;
}

=head1 object()

Bio::SeqFeature:Interaction

=head1 DESCRIPTION
	set or get object values of an Bio::SeqFeature:Interaction object
	must be a Bio::SeqFeatureI
=cut

sub subject {
	my $self = shift;
	if (@_) {
		my $subject = shift;
		if (_isSeqFeatureI ($subject) ) { #correct type?
			$self->{'_subject'} = $subject;
		}else{ print " Subject \n";}
	}else {
		return $self->{'_subject'};
	}
}

=head1 subject()

Bio::SeqFeature:Interaction

=head1 DESCRIPTION
	set or get subject values of an Bio::SeqFeature:Interaction object
	must be a Bio::SeqFeatureI
=cut

sub object {
	my $self = shift;
	if (@_) {
		my $object = shift;
		if (_isSeqFeatureI ($object) ) { #correct type?
			$self->{'_object'} = $object;
		}else{ print " Object \n";}
	}else {
		return $self->{'_object'};
	}
}


# control if is a SeqFeatureI

sub _isSeqFeatureI {
	my $subject = shift;
	unless ($subject->isa('Bio::SeqFeatureI')) {
		carp("Error : object must be a Bio::SeqFeatureI");
		return 0;
	}
	return 1;
}

# control if a value is define

sub _isDefined{
	my $val = shift;
	if (defined $val) { print return 1;} else {return 0;}
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

=head1 direction()

Bio::SeqFeature:Interaction

=head1 DESCRIPTION
	get direction values of an Bio::SeqFeature:Interaction object

=cut

sub direction {

	my $self = shift;
	my $direction = $self->_get_direction();
	return $direction;
}

# protected method set_direction()

sub _set_direction {
	my $self = shift;
	if (@_) {
		$self->{'_direction'} = shift;
	}else {
		croak ("direction undefined, cannot set direction \n"); # i want it to die because it's a protected method; then the developper make an error
	}
}

# protected method get_direction()

# Bio::SeqFeature:Interaction

sub _get_direction {
	my $self = shift;
	if (_isDefined($self->{'_direction'}) == 0 ) { #1st call is not defined
		$self->_set_direction($self->_calculate_direction);
		my $direction = $self->_get_direction();
		$self->_set_direction($direction);
	}
	return $self->{'_direction'};
}

# private calculate the direction of the interaction

sub _calculate_direction {
	my $self = shift;
	my $object = $self->object();
	my $subject = $self->subject();
	if ( _isComplete($self)== 1){
		if ($subject->strand() == 0 || $object->strand() == 0 ){ # unknow
			return(0);
		}elsif( $subject->strand() != $object->strand() ) { # anti direction
			return(-1);
		}else { #direction
			return(1);
		}
	} else { croak (" Interaction : calculate_direction : You must have an object and a subject to have calculate a direction \n ");}
}
=head1 type()

Bio::SeqFeature:Interaction

=head1 DESCRIPTION
	get type of an Bio::SeqFeature:Interaction object
=cut
sub type {

	my $self = shift;
	my $type = $self->_get_type();
	return $type;
}
# protected method set_type()
sub _set_type {
	my $self = shift;
	if (@_) {
		$self->{'_type'} = shift;
	}else {
		croak ("Type undefined, cannot set type value \n"); # i want it to die because it's a protected method; then the developper make an error
	}
}

# protected method get_direction()

# Bio::SeqFeature:Interaction

sub _get_type {
	my $self = shift;
	if (! defined $self->{'_type'}) { #1st call is not defined
		$self->_set_type($self->_calculate_type);
		my $type = $self->_get_type();
		$self->_set_type($type); #set the type
	}
	return $self->{'_type'};
}

# private calculate the type of the interaction

sub _calculate_type {
	my $self = shift;
	my $object = $self->object();
	my $subject = $self->subject();
	if (_isComplete($self) == 1 ){
		if ($subject->overlaps($object) == 1){
			return(1); # genic
		}else {
			return(0); #intergenic
		}
	} else { croak (" Interaction : calculate_direction : You must have an object and a subject to have calculate a direction \n ");}
}

=head1 distance()

Bio::SeqFeature:Interaction

=head1 DESCRIPTION
	get distance inside a complete interaction : object and subject are store
	if type is not set it will be calculate and store
=cut

sub distance {
	my $self = shift;
	my $distance = $self->_get_distance();
	return $distance;
}
# protected

sub _get_distance {
	my $self = shift;
	my $object = $self->object();
	my $subject = $self->subject();

	if(_isDefined($self->{'_distance'}) == 0) { #not defined
		$self->_set_distance($self->_calculate_distance());
	}
	return ($self->{'_distance'});
}

sub _calculate_distance{
	my $self = shift;
	if (_isComplete($self) == 0 ){
		croak (" Interaction : get_distance : calculate_distance : You must have an object and a subject to have calculate a distance \n ");
	}else{
		my $type = $self->type();
		my $distance = 0;
		if ( $type == 0 ) {
			$distance = $self->{'_object'}->start() - $self->{'_subject'}->end();
			if ( $distance < 0 ) {
				$distance = $self->{'_subject'}->start() - $self->{'_object'}->end();
			}
		}
		return $distance;
	}
}
# protected

sub _set_distance {
	my $self = shift;
	if (@_) {
		$self->{'_distance'} = shift;
	}else{
		croak (" Interaction : set_distance : There must have an value \n ");
	}
}


=head1 overlap()

returns the number of bases that are overlapping for a genic interaction

=head1 DESCRIPTION
this fucntions calculates the number of bases that overalps between the subject and object

=cut

sub interaction_overlap {
	my $self = shift;
	my $overlap = $self->_overlap();
	return $overlap;
}
# protected

sub _get_overlap {
	my $self = shift;
	my $object = $self->object();
	my $subject = $self->subject();

	if(_isDefined($self->{'_overlap'}) == 0) { #not defined
		$self->_set_overlap($self->_calculate_overlap());
	}
	return ($self->{'_overlap'});
}


sub _calculate_overlap{
	my $self = shift;

	if (_isComplete($self) == 0 ){
		croak (" Interaction : get_distance : calculate_distance : You must have an object and a subject to have calculate a distance \n ");
	}
	else{
		my $type = $self->type();
		my $overlap =0;
		if ( $type == 1) {
			$overlap=$self->overlap();
		}
		return $overlap;
	}
}
# protected

sub _set_overlap {
	my $self = shift;
	if (@_) {
		$self->{'_overlap'} = shift;
	}else{
		croak (" Interaction : set_overlap needs a value \n ");
	}
}


=head1 printer_mini()

Bio::SeqFeature:Interaction

=head1 DESCRIPTION
	$interaction->printer_mini()
	print minimal interaction informations
=cut

sub printer_mini {
	my $self=shift;
	my $best = shift;
	my $biotype = shift;

	# VW: variable best
	my $bestVal = 0;
	# VW: modif
	if (defined $best) {$bestVal = 1}

	# VW: if the transcripts get a biotype print it
	if($biotype)
	{
	    print join ("\t", $bestVal, $self->object()->get_tag_values("gene_id"), $self->object()->get_tag_values("transcript_id"), $self->object()->get_tag_values("transcript_biotype"), $self->subject()->get_tag_values("gene_id"), $self->subject()->get_tag_values("transcript_id"), $self->subject()->get_tag_values("transcript_biotype"), _conversion_direction($self->direction()), _conversion_type($self->type()), $self->distance());
	}
	else
	{
	    print join ("\t", $bestVal, $self->object()->get_tag_values("gene_id"), $self->object()->get_tag_values("transcript_id"), $self->subject()->get_tag_values("gene_id"), $self->subject()->get_tag_values("transcript_id"), _conversion_direction($self->direction()), _conversion_type($self->type()), $self->distance());
	}
}


=head1 printer()

Bio::SeqFeature:Interaction

=head1 DESCRIPTION
	$interaction->printer()
	print interaction informations
=cut

sub printer {
	my $self = shift;
	my $inversion = _sayme_switching($self->{'_object'});

#	print "\t Object \n";
#	_ligne_carre();
	_print_finfo($self->{'_object'});
	# print object specific informations according to his primary tag
		if ($inversion == 1){
			_lncRNA_informations($self->{'_object'});
		}else{
			_mRNA_informations($self->{'_object'});
		}


	print "\t Subject \n";
	_ligne_carre();
	_print_finfo($self->{'_subject'});

	# print subject specific informations according to his primary tag
		if ($inversion == 1){
			_mRNA_informations($self->{'_subject'});
		}else{
			_lncRNA_informations($self->{'_subject'});
		}

	if ($self->{'_subject'}->seq_id() eq "-555") {

		print "\t Interaction information \n";
		_ligne_carre();
		print "\t Direction : ", _conversion_direction($self->direction()),"\n";
		print "\t Type : ", _conversion_type($self->type()),"\n";
		print "\t Distance : ";
		if( _isComplete($self) == 1) { print ("",$self->distance(), "\n");} else { print "NULL \n";}
	}
	print "Yahohooh\n";
	_ligne_carre();
}
sub _mRNA_informations{
	my $elem = shift;
	print "\t || Transcript id : ",$elem->get_tag_values("transcript_id"), "\n";

}
sub _lncRNA_informations{
	my $elem = shift;
	print "\t || Transcript id : ",$elem->get_tag_values("transcript_id"), "\n";
    #print "\t || Gene id : ",$elem->get_tag_values("gene_id"), "\n";

}
# when you switch the object and subject , then, theirs attributes are not the same, the bioPerl class
# raise fatal error
sub _sayme_switching {
	my $objet = shift;
	if ($objet->primary_tag() eq "lncRNA" ) {return 1;} # object is lncRNA
	return 0; # object is mRNA
}
sub _ligne_carre(){
	print "\t ------\t------\t------\n";
}

# print a feature info, test if the feature is defined

sub _print_finfo{
	#print "il y a ", scalar(@_), "element passé \n";
	my $f = shift;
	if (_isDefined($f) == 0 ){
		print " You want to print informations about an undefined feature \n";
		exit;
	}
	print " \t    # Info feature \n   ";
#	_ligne_carre();
	print	"\t | Primary tag : ",
    		$f->primary_tag(),"\n \t | source tag : ",
#     		$f->source_tag,"\n \t | seq id : ",
    		$f->seq_id() , "\n \t | start : ",
    		$f->start, "\n \t |  end : ",
    		$f->end, "\n \t | strand : ",
    		$f->strand ,"\n  ";
  #  _ligne_carre();
}

# transform the direction interger value in the appropriate direction string value <<< lawl
# 0 <=> direction
# 1 <=> anti direction
sub _conversion_direction {
	my $i= shift;
	my %direction = ("1","sense","-1","antisense","0","strand_unknow");
	if (exists $direction{$i} ) {
		return ($direction{$i});
	}else{ return "$i unknow direction";}
}


# transform the type interger value in the appropriate type string value <<< lawl
# 0 <=> intergenic
# 1 <=> genic

sub _conversion_type {
	my $i = shift;
	my @TYPES= ("intergenic", 'genic');
	return ($TYPES[$i]);

}

sub isGenic {
	my $self = shift;
	my $type = $self->type();
	if ($type == 1){ return 1;} else { return 0;}

}
sub isInterGenic{
	my $self = shift;
	my $type = $self->type();
	if ($type == 1){ return 0;} else { return 1;}
}

## print the object tags list
sub tags_list {
	my $self = shift;
	my @ks = keys(%{$self});
	foreach my $k (@ks) {
		print " $k \n";
	}
}
### give me an tag value

sub get_tag_value{
	my $self  = shift; #interaction
	my $tag   = shift; # the tag

	if (_isExist_hash($self,$tag) == 1){
		return $self->{$tag};
	}
	#$self->tags_list();
	return -555; #"tag does not exist";

}

# look for the existance of a key in a hash
sub _isExist_hash {
	my %hash = %{$_[0]};
	my $cle = $_[1];

	if (exists $hash{$cle}){
		return 1;
	}

	return 0;
}
1;
