package Filter;

# 
$VERSION = v0.0.1;

use warnings;
use strict;
use ExtractFromHash;
use ExtractFromFeature;
use Data::Dumper;

$| = 1;


# if we compare all features from transcripts
sub filter{

	my ($refh1, $refh2, $uniquexa, $uniquexb, $common, $fraction, $verbosity)	= @_;
	my %h1			=	%{$refh1}; # parsing a hash in sub need dereference in shift
	my %h2			=	%{$refh2}; # parsing a hash in sub need dereference in shift
	$fraction 		||= 1;
	$verbosity 		||= 0;

# 	print STDERR "$refh1, $refh2,  $keepduptx, $keepdupint, $verbosity\n";

	# variable that should be options
	my $stranded 		= 1; #require  same strandness for overlap	
    
    # nb tx1 and
	my $nb_tx1 = keys( %h1 );
	my $i = 0;

	my $nb_tx2 = keys( %h2 );
	my $j = 0;	    
    
   	my %commonA=();  	
   	my %commonB=();  	   	
   	
    # Sort both hash by chr and transcript start ie startt
    foreach my $tr1 (sort {$h1{$a}->{'startt'} <=> $h1{$b}->{'startt'} }keys %h1) {

		# infos
		$i++;
		my $s1 =	$h1{$tr1}->{"strand"};
		
		# verbose on file1
	    if ($verbosity > 10){ Utils::showProgress( $nb_tx1, $i++, "Analyse fileA: ");}


	    # Sort both hash by chr and transcript start ie startt
        foreach my $tr2 (sort { $h2{$a}->{'startt'} <=> $h2{$b}->{'startt'} }keys %h2) {
			
			my $s2 =	$h2{$tr2}->{"strand"};

			# trick to speed (?) loop
			last if ($h2{$tr2}->{"startt"} > $h1{$tr1}->{"endt"});
			next if ($h2{$tr2}->{"endt"} < $h1{$tr1}->{"startt"});            

			##### Overlap window  ########
			
			# ############## Exon overlap
			my $tr1_exon_size		= ExtractFromHash::cumulSize($h1{$tr1}->{"feature"});
    		my $tr2_exon_size		= ExtractFromHash::cumulSize($h2{$tr2}->{"feature"});
			
			# Exon overlap
			my $cumul_overlap_size	=	ExtractFromFeature::intersectFeatures($h1{$tr1}->{'feature'}, $h2{$tr2}->{'feature'}, $s1, $s2, $stranded, $verbosity);
			
			# proportion
			my $fractionoverexon1	=	$cumul_overlap_size/$tr1_exon_size;
			my $fractionoverexon2	=	$cumul_overlap_size/$tr2_exon_size;
							
			# remove same tx exons in 2
			if ( $fractionoverexon1 >= $fraction  && $fractionoverexon2 >= $fraction ){
				$commonA{$tr1}++;
				$commonB{$tr2}++;
				
				print STDERR "SameExons:	$tr1 - $tr2...\n";
				
			}
        }
	}
	################
	
	my %hfinal;
	# return 
	if ($uniquexa){
		for my $tx1 (keys %h1){
		
			if (!exists($commonA{$tx1})){
		
				$hfinal{$tx1}	=	$h1{$tx1};
			}
		
		}
	} elsif ($uniquexb){
		for my $tx2 (keys %h2){
		
			if (!exists($commonB{$tx2})){
		
				$hfinal{$tx2}	=	$h2{$tx2};
			}
		
		}
	}
	return %hfinal;
}


1;
