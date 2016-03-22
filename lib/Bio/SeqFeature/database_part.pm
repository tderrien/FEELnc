package Bio::SeqFeature::database_part;
use strict;
use warnings;

use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;

use Bio::Tools::GFF;
use Bio::SeqIO;

use Data::Dumper;

my $BASEDIR="/tmp";


=head1  new
	 return : Bio::DB::SeqFeature::Store
=head1 DESCRIPTION
	load a complete GTF file in a parameterizable DATABASE
=cut

sub new {
	my $pkg = shift;

	my $rand_number =  int(rand(10000));
	my $dsn =  $BASEDIR."/interact.$rand_number";

	my $db = Bio::DB::SeqFeature::Store->new( -adaptor => 'berkeleydb',
                                            -dsn     => "$BASEDIR/interact.$rand_number",
                                            -create => 1,
                                            -compress => 1 );
	 my $db_part = bless {
		 _db => $db,
		 _rand => $rand_number,
		 _dsn => $dsn
		}, $pkg;

	return $db_part;
}

=cut


=head1 db

	returns the database

=cut

sub db {
	my $self=shift;

	return $self->{_db};
}

=head1 dsn

	returns disk location of the db

=cut

sub dsn {
	my $self=shift;

	return $self->{_dsn};
}


=head1 store

	stores a feature into the database

=cut

sub store {
	my $self=shift;
	my $feat = shift;

	return $self->db()->store($feat);
}


sub get_seq_stream {
	my $self=shift;
	$self->db()->get_seq_stream(@_);
}


sub get_features_by_location{
	my $self=shift;
	$self->db()->get_features_by_location(@_);
}

=head1 destroy
	destroy the database
=cut

sub destroy {
	my $self=shift;

	system ("rm -r " . $self->dsn());
}


=head1 load_complete_gff_into_DB
	args:
		$gff_file
	return : Bio::DB::SeqFeature::Store
=head1 DESCRIPTION
	load a complete GFF files in a DATABASE
=cut

sub load_complete_gff_into_db {

	my $self = shift;
	my $file = shift;


	print STDERR " the file is : ", $file, "\n";
	my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(-store    => $self->db(),
                                                           -verbose  => 1,
                                                           -fast     => 1);


	$loader->load($file); # chargement annotations db

}
=cut

=head1 load_complete_gtf_into_DB
	arg :
		- database
		- file
	 return:  Bio::DB::SeqFeature::Store
=head1 DESCRIPTION
	load a complete GTF file in a DATABASE
=cut

sub load_complete_gtf_into_db {
	my $self = shift;
	my $gtf = shift;
	my $format = shift;
	my $primary = shift;


	print STDERR  " the file is : ", $gtf, " tagged as ", $primary, "\n";
	my $input = Bio::Tools::GFF->new(-file => $gtf, -gff_version => $format);
	my $line = 0;
	while ( my $feature = $input->next_feature() ) {
		$line++;
		$feature->set_attributes(-primary => $primary); #it's an mRNA and LncRNAs
		$self->store($feature) or die "Importation complete GTF failed line : $line \n";
	}
	print STDERR " Importation  complete GTF : done \n";
}

=head1 load_merge_gtf_into_DB
	Bio::DB::SeqFeature::Store::GFF3Loader
	args:
		$file : from merge step
	return : Bio::DB::SeqFeature::Store
=head1 DESCRIPTION
	load a GTF file come from the bedtools merge in a DATABASE
=cut

## in such file transcript lines are missing and need to be correctly reconstruct
## reconstruct transcript - from theirs exons - and attach them theirs exons subsequences

sub load_merge_gtf_into_db {


	my $self = shift;
	my $gtf = shift;
	my $format = shift;
	my $type= shift;

	my $input= Bio::Tools::GFF->new(-file => $gtf, -gff_version => $format);  #read gtf file - lncRNAs
	my $line = 0;
	my %transcrit_id=();

	my $trash=0;
	my %tx2gene_id; # correspondence gene_id to tx_id for gff

	# first step : link together all exon of a transcript_id


	while ( my $feature = $input->next_feature() ) {

		my $cle;
		my $gene_id;
		if ($format == 2) {
			next unless ($feature->primary_tag eq 'exon');
 			($cle) 		= $feature->get_tag_values("transcript_id"); # nom du transcrit

 			# if some line does not have gene_id (bedtools converter tx VW:)
 			# put gene_id == transcript_id
 			if ( ! $feature->has_tag("gene_id") ){
 				$feature->add_tag_value('gene_id', $cle);
 			}

			# VW: if the transcript didn't have a transcript_biotype, set it to 'none'
			if( ! $feature->has_tag("transcript_biotype") )
			{
			    $feature->add_tag_value('transcript_biotype', 'none');
			}
		}
		if ($format == 3) {

			# store gene_id transcript_id in hash to add tag gene_id for gff format
			if ($feature->primary_tag eq 'mRNA'){
				my ($gene_id)			= 	$feature->get_tag_values("Parent");
 				my ($tx_id)			= 	$feature->get_tag_values("ID");
 				$tx2gene_id{$tx_id}	=	$gene_id;
				next;


			} elsif($feature->primary_tag eq 'exon'){
				($cle) = $feature->get_tag_values("Parent"); # nom du transcrit
				$feature->add_tag_value('transcript_id', $cle);


			} else {
				next;
			}
		}
 		push(@{$transcrit_id{$cle}},$feature); # add to array
	}




	# ad gene_id tag for gff for consistency in print
	if ($format == 3){
		foreach  my $k (keys(%transcrit_id))  {
			$transcrit_id{$k}->[0]->set_attributes(-tag =>{"gene_id" => $tx2gene_id{$k}});
		}
# 		print Dumper \%transcrit_id;
	}

	# second step : foreach exons array make the corresponding transcript with correct start and ends
	foreach  my $k (keys(%transcrit_id))  {

	   	my $feat = _copy_feature($type, $transcrit_id{$k}->[0]); #copy the first exon

 		# not sure it is important since it is in the tag is cp in primary
	   	$feat->set_attributes(-type => $type); #it's an lncRNA or mRNA

	   	my ($start,$end) = _position_from_farray (@{$transcrit_id{$k}}); # give me the position
	   	$feat->start($start)  ;
	   	$feat->end($end)   ;

		# third step : attach to each transcript SeqFeature his exons SeqFeature
		foreach my $f (@{$transcrit_id{$k}}){

	   		#might check if it's correctly reconstruct
	   		if ($feat->start() > $f->start () or $feat->end() < $f->end() ){
	   			croak (" Fatal error in  start/end features checking...\n");
	   		}

	   	   	$feat->add_SeqFeature($f);
	   	}

	# fourth step : add to database
	   	$line ++;
	   	#print "Lets add this\n";
	   	#print STDERR "$line\n";
	   	$self->store($feat) or die " GTF FILE from merge IMPORTATION : Error entry $line \n";
	 	#sleep(2);
		}

		print STDERR  " Importation  merge GTF file : ", $line, " $type done \n";
		return $line;
	#last : return the $db
}

# private method : _copy_feature
# function : deep copy a feature object
# arg : a feature object
# return a Bio::SeqFeature::Generic object

sub _copy_feature{
	# imperfect copy: source tag modified, primary modified
	my $type = shift;
	my $feat = shift (@_);


	my $cfeat= new Bio::SeqFeature::Generic (
	    -seq_id  => $feat->seq_id(),
	    -start   => $feat->start(), -end => $feat->end(),
	    -strand  => $feat->strand(),
	    -primary => $type,
	    -score   => $feat->score(),
	    -tag     => {
		"transcript_id" => $feat->get_tag_values("transcript_id"),
		"gene_id"       => $feat->get_tag_values("gene_id"),
		"transcript_biotype" => $feat->get_tag_values("transcript_biotype"), #VW: adding the transcript_biotype into the deep copy
	    });;
	# "class_code" => $feat->get_tag_values("class_code"),
	#"gene_id" => $feat->get_tag_values("gene_id"),
	#"oId" => $feat->get_tag_values("oId"),
	#"transcript_id" => $feat->get_tag_values("transcript_id"),
	#"tss_id" => $feat->get_tag_values("tss_id"),
	#} );;

   	#_afficher_tag_names($cfeat);
  	#	print " je l'ai recopié en : \n";

	#   _afficher_fvalues($cfeat);

	return $cfeat;
}

# private method : _position_from_farray
# function : Get the greater end and the smaller start from feature array
# arg : array of features
# return: $start $end

sub _position_from_farray {
	#print "___ position from feature aray ___ \n";
	my @t = @_;

	my $start= $t[0]->start();
	my $end = $t[0]->end();
	foreach my $p (@t){ #parcours tableau
		if ( $start > $p->start() ) {
			$start = $p->start;
		}

	# end position
		if ( $end < $p->end() ) {
			$end = $p->end;
		}
	}

	return ($start,$end);
}

# private method :  _afficher_fvalues
# function : Print feature informations
# arg : $feature
# return: /

sub _afficher_fvalues{
    my $feat = shift @_;
    print " **** affichage ----- \n";
    print " scaffold / seq_id:",  $feat->seq_id() , "\n";
    print " start :",  $feat->start() , "\n";
    print " end :" , $feat->end() , "\n";
    print " source_tag:" , $feat->source_tag() , "\n";
    print "\t \t";
    print " hash tag  -tag    => {clé => val}\n";
    #print "\t \t";print " class code : ", $feat->get_tag_values("class_code"), "\n";
    print "\t \t";print " gene id : ", $feat->get_tag_values("gene_id"), "\n";
    #print "\t \t";print " oId : ", $feat->get_tag_values("oId"), "\n";
    #print "\t \t";print " tss_id: ", $feat->get_tag_values("tss_id"), "\n";
    print "\t \t";print " transcript_id: ", $feat->get_tag_values("transcript_id"), "\n";
    print "\t \t";print " transcript_biotype: ", $feat->get_tag_values("transcript_biotype"), "\n";
    print " **** fin affichage ----- \n";
}





# private method :  _afficher_fvalues
# function :  print tags names
# arg : $feature
# return: /

sub _afficher_tag_names{
	print "------------------------- affichage tags  --------------------- \n";
	my $feature= shift (@_);
	my @tags = $feature->get_all_tags();

	foreach my $t (@tags){
		print "tags : ", $t ,"\n";
	}
		print "------------------------- FIN --------------------- \n";
}




1;
