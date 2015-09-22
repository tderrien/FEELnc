package RandomForest;

$VERSION = v0.0.1;

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use Utils;
use ExtractFromFeature; # getKeyFromFeature line 932
use Bio::DB::Fasta;
use Bio::SeqIO;

# Perl libs
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Data::Dumper;

# lib directory : ~tderrien/bin/perl/lib/
use Parser;
use ExtractFromHash;
use ExtractFromFeature;
use Intersect;
use Utils;
use Orf;


=encoding UTF-8

=head1 RandomForest.pm

=over

=item .

getKmerRatio: Get one kmer ratio from fasta training file

=item .

getKmerRatioSep: Get one kmer ratio between coding and non coding fasta training file

=item .

scoreORF: Scores ORF test file

=item .

getSizeFastaFile: Get the size of each element of a multi fasta and write it in a file

=item .

mergeKmerScoreSize: Merge a list of kmer score file, a file with the size of ORF and a file with the size of mRNA

=item .

getRunModel: Get a threshold using 10-fold cross-validation on training data if it is not defined and crate and apply the model to predict biotypes for test sequences using an R script


=item .

runRF: Run all the process using ORF from coding sequences, fasta from coding and non coding sequences, ORF and fasta from test sequences, a list of kmer size and a threshold (if defined)

=item .

rfPredToOut: With a .gtf and a result file from a random forest, write 2 .gtf file, each one respectively for coding and non coding genes

=back

=cut


# Run KmerInShort (kis) on training set and generate output model (logRatio values for each kmer)
#	$codFile   = fasta file with the ORF of coding genes
#	$nonFile   = fasta file with sequence of non coding genes
#	$outFile   = file where the logRatio of kmer frequency need to be written
#	$kmerSize  = size of the kmer
#	$codStep   = step for the counting of kmer for the ORF coding genes
#	$nonStep   = step for the counting of kmer for the non coding genes
#	$proc      = number of proc to be use for KmerInShort
#	$verbosity = value to define the verbosity
#       $nameTmp   = absolute path and prefix for the temporary files
#	$keepTmp   = keeping or not the temporary files
#	Return value:
#		Return the array of the log ratio value
sub getKmerRatio
{
    my($codFile, $nonFile, $outFile, $kmerSize, $codStep, $nonStep, $proc, $verbosity, $nameTmp, $keepTmp) = @_;
    $codFile ||= undef;
    $nonFile ||= undef;
    $outFile ||= undef;
    $codStep ||= 3;
    $nonStep ||= 1;
    $proc    ||= 2;
    $keepTmp ||= 0;

    if($kmerSize < $codStep)
    {
	$codStep = 1;
	$nonStep = 1;
    }

    # Check if mendatory arguments have been given
    die "Bulding model: ORF coding genes training file not defined... exiting\n" if (!defined $codFile);
    die "Bulding model: lncRNA training file not defined... exiting\n"           if (!defined $nonFile);
    die "Bulding model: output file for kmer model is missing... exiting\n"      if (!defined $outFile);

    # emtpy
    die "Bulding model: ORF coding genes training file '$codFile' is empty... exiting\n" unless (-s $codFile);
    die "Bulding model: lncRNA training file '$nonFile' is empty... exiting\n"           unless (-s $nonFile);

    # Path to kis
    my $kisPath = Utils::pathProg("KmerInShort");
    my $cmd = "";
    # Temporary files to put kmer counting for ORF coding genes and non coding genes
    my $codOut = $nameTmp.".coding_size".$kmerSize."_kmerCounting.tmp";
    my $nonOut = $nameTmp.".noncoding_size".$kmerSize."_kmerCounting.tmp";

    # Run kis on ORF for coding genes
    print "\tRunning KmerInShort on '$codFile'.\n" if($verbosity >= 5);
    $cmd = "$kisPath -file $codFile -nb-cores $proc -kmer-size $kmerSize -out $codOut -dont-reverse -step $codStep 1>/dev/null 2>/dev/null";
    system($cmd);
    #warn "$cmd"
    # Run kis on non coding genes
    print "\tRunning KmerInShort on '$nonFile'.\n" if($verbosity >= 5);
    $cmd = "$kisPath -file $nonFile -nb-cores $proc -kmer-size $kmerSize -out $nonOut -dont-reverse -step $nonStep 1>/dev/null 2>/dev/null";
    system($cmd);

    # Read the two kmer files and put value in a table and get the total number of kmer to comput frequency
    my @kmerTab;
    my @codVal;
    my $codTot = 0;
    my @nonVal;
    my $nonTot = 0;
    my $i;

    # Read the kmer counting for ORF on coding genes
    print "\tRead kmer counting for coding genes of size '$kmerSize' in '$codOut'.\n" if($verbosity >= 5);
    # Open file
    open FILE, "$codOut" or die "Error! Cannot open kmerFile '". $codOut . "': ".$!;
    $i = 0;
    while(<FILE>)
    {
	chop;
	my ($kmer, $val) = split(/\t/);
	#$kmerTab[$i]     = $kmer;
	#$codVal[$i]      = $val;
	push(@kmerTab,    $kmer);
	push(@codVal, int($val));
	$codTot = $codTot + int($val);
	$i = $i+1;
    }
    close FILE;

    # Read the kmer counting for non coding genes
    print "\tRead kmer counting for non coding genes of size '$kmerSize' in '$nonOut'.\n" if($verbosity >= 5);
    # Open file
    open FILE, "$nonOut" or die "Error! Cannot open kmerFile '". $nonOut . "': ".$!;
    $i = 0;
    while(<FILE>)
    {
	chop;
	my ($kmer, $val) = split(/\t/);

	die "Error: kmer order is not the same between '$codOut' and '$nonOut'.\nExit." if($kmerTab[$i] ne $kmer);

	#$nonVal[$i] = $val;
	push(@nonVal, int($val));
	$nonTot = $nonTot + int($val);
	$i = $i+1;
    }
    close FILE;

    # Write the output file with the log ratio directly
    print "\tWrite the log ratio of kmer between coding and non coding genes frequency for a size of kmer of '$kmerSize' in '$outFile'.\n" if($verbosity >= 5);
    open FILE, "> $outFile" or die "Error! Cannot access output file '". $outFile . "': ".$!;
    my $nbKmer = @kmerTab;
    my $log    = 0;

    # VW No header for KmerInShort
    # Print the file header
    #print FILE "kmer\tkmerSize_logRatio\n";
    for($i=0; $i<$nbKmer; $i++)
    {
	## VW MODIF
	# # logratio = 0                      -- if the kmer is not found in any gene files
	# if($codVal[$i]==0 && $nonVal[$i]==0)
	# {
	#     $log = 0;
	# }
	# # logratio = 1                      -- if the kmer is found in coding ORFs but not in non coding genes
	# elsif($codVal[$i]>0 && $nonVal[$i]==0)
	# {
	#     $log = 1;
	# }
	# # logratio = -1                     -- if the kmer is not found in coding ORFs but found in non coding genes
	# elsif($codVal[$i]==0 && $nonVal[$i]>0)
	# {
	#     $log = -1;
	# }
	# # logratio = log( codFreq/nonFreq ) -- if the kmer is found on the two files
	# else
	# {
	#     $log = log( ($codVal[$i]/$codTot) / ($nonVal[$i]/$nonTot) );
	# }

	$log = log( (($codVal[$i]/$codTot)+1) / (($nonVal[$i]/$nonTot)+1) );

	#$logTab[$i] = $log;
	#print FILE "$kmerTab[$i]\t$log\n";
	# VW only the log for KmerInShort
	print FILE "$log\n";
    }
    close FILE;

    # Delete the temporary files if keepTmp != 0
    if($keepTmp == 0)
    {
	unlink $codOut;
	unlink $nonOut;
    }

    return(1);
}



# Run KmerInShort (kis) on training set and generate output model (logRatio values for each kmer depending only on the kmer)
#	$codFile   = fasta file with the ORF of coding genes
#	$nonFile   = fasta file with sequence of non coding genes
#	$outFile   = file where the logRatio of kmer frequency need to be written
#	$kmerSize  = size of the kmer
#	$codStep   = step for the counting of kmer for the ORF coding genes
#	$nonStep   = step for the counting of kmer for the non coding genes
#	$proc      = number of proc to be use for KmerInShort
#	$verbosity = value to define the verbosity
#       $nameTmp   = absolute path and prefix for the temporary files
#	$keepTmp   = keeping or not the temporary files
#	Return value:
#		Return the array of the log ratio value
sub getKmerRatioSep
{
    my($codFile, $nonFile, $outFile, $kmerSize, $codStep, $nonStep, $proc, $verbosity, $nameTmp, $keepTmp) = @_;
    $codFile ||= undef;
    $nonFile ||= undef;
    $outFile ||= undef;
    $codStep ||= 3;
    $nonStep ||= 3;
    $proc    ||= 1;
    $keepTmp ||= 0;

    if($kmerSize < $codStep)
    {
	$codStep = 1;
	$nonStep = 1;
    }

    # Check if mendatory arguments have been given
    die "Bulding model: ORF coding genes training file not defined... exiting\n" if (!defined $codFile);
    die "Bulding model: lncRNA training file not defined... exiting\n"           if (!defined $nonFile);
    die "Bulding model: output file for kmer model is missing... exiting\n"      if (!defined $outFile);

    # emtpy
    die "Bulding model: ORF coding genes training file '$codFile' is empty... exiting\n" unless (-s $codFile);
    die "Bulding model: lncRNA training file '$nonFile' is empty... exiting\n"           unless (-s $nonFile);

    # Path to kis
    my $kisPath = Utils::pathProg("KmerInShort");
    my $cmd = "";
    # Temporary files to put kmer counting for ORF coding genes and non coding genes
    my $codOut = $nameTmp.".coding_size".$kmerSize."_kmerCounting.tmp";
    my $nonOut = $nameTmp.".noncoding_size".$kmerSize."_kmerCounting.tmp";

    # Run kis on ORF for coding genes
    print "\tRunning KmerInShort on '$codFile'.\n" if($verbosity >= 5);
    $cmd = "$kisPath -file $codFile -nb-cores $proc -kmer-size $kmerSize -out $codOut -dont-reverse -step $codStep 1>/dev/null 2>/dev/null";
    system($cmd);
    if ($? != 0)
    {
        die "\nFailed to run KmerInShort: check your KmerInShort PATH (LINUX/MAC). Command line:\n$cmd\n";
    }

    # Run kis on ORF non coding genes
    print "\tRunning KmerInShort on '$nonFile'.\n" if($verbosity >= 5);
    $cmd = "$kisPath -file $nonFile -nb-cores $proc -kmer-size $kmerSize -out $nonOut -dont-reverse -step $nonStep 1>/dev/null 2>/dev/null";
    system($cmd);
    if ($? != 0)
    {
        die "\nFailed to run KmerInShort: check your KmerInShort PATH (LINUX/MAC). Command line:\n$cmd\n";
    }

    # Read the two kmer files to get the kmer frequency and write the ratio in the output file
    my @codTab;
    my @nonTab;
    my $codKmer = "";
    my $nonKmer = "";
    my $codVal  = 0;
    my $nonVal  = 0;
    my $codTot  = 0;
    my $nonTot  = 0;
    my $ratio   = 0;

    open FILECOD,   "$codOut"  or die "Error! Cannot open kmerFile '". $codOut . "': ".$!;
    open FILENON,   "$nonOut"  or die "Error! Cannot open kmerFile '". $nonOut . "': ".$!;
    open FILEOUT, "> $outFile" or die "Error! Cannot access output file '". $outFile . "': ".$!;
    print "\tWrite the ratio for each kmer between coding and non coding kmer counting for a size of kmer of '$kmerSize' in '$outFile'.\n" if($verbosity >= 5);

    # Read the two kmer files to get the kmer frequency
    while(not eof FILECOD and not eof FILENON)
    {
	($codKmer, $codVal) = split(/\t/, <FILECOD>);
	($nonKmer, $nonVal) = split(/\t/, <FILENON>);

	die "Error: kmer order is not the same between '$codOut' and '$nonOut' ($codKmer, $nonKmer).\nExit." if($codKmer ne $nonKmer);

	push(@codTab, int($codVal));
	push(@nonTab, int($nonVal));
	$codTot = $codTot + $codVal;
	$nonTot = $nonTot + $nonVal;
    }

    # Write the output file
    my $nbKmer = @codTab;
    my $codR   = 0;
    my $nonR   = 0;

    for(my $i=0; $i<$nbKmer; $i++)
    {
	$codR = $codTab[$i]/$codTot;
	$nonR = $nonTab[$i]/$nonTot;

	if(($codR+$nonR) == 0)
	{
	    $ratio = 0.5;
	}
	else
	{
	    $ratio = ($codR) / ($codR + $nonR);
	}
	print FILEOUT "$ratio\n";

    }
    close FILECOD;
    close FILENON;
    close FILEOUT;

    # Delete the temporary files if keepTmp != 0
    if($keepTmp == 0)
    {
	unlink $codOut;
	unlink $nonOut;
    }

    return(1);
}


# Function to compute the scoring function of a multifasta file of ORF
#	$orfFile   = multi fasta file with the ORF of the genes to be scored
#	$modFile   = file with the log ratio score compute on learning files
#	$outFile   = file where the logRatio of kmer frequency need to be written
#	$kmerSize  = size of the kmer
#	$step      = step for the counting of kmer
#	$proc      = number of proc to be use for kis
#	$verbosity = value to define the verbosity
sub scoreORF
{
    my($orfFile, $modFile, $outFile, $kmerSize, $step, $proc, $verbosity) = @_;
    $orfFile  ||= undef;
    $modFile  ||= undef;
    $outFile  ||= undef;
    $kmerSize ||= 6;
    $step     ||= 3;
    $proc     ||= 1;

    if($kmerSize < $step)
    {
	$step = 1;
    }

    # Check if mendatory arguments have been given
    die "Scoring ORF file: ORF file for scoring sequences is not defined... exiting\n"                     if(!defined $orfFile);
    die "Scoring ORF file: model of log ratio score file is not defined... exiting\n"                      if(!defined $modFile);
    die "Scoring ORF file: output file to write kmer score for test sequences is not defined... exiting\n" if(!defined $outFile);

    # emtpy
    die "Scoring ORF file: ORF file for scoring sequences '$orfFile' is empty... exiting\n" unless (-s $orfFile);
    die "Scoring ORF file: model of log ratio score file '$modFile' is empty... exiting\n"  unless (-s $modFile);

    # Path to kis
    my $kisPath = Utils::pathProg("KmerInShort");

    print "\tRun KmerInShort on '$orfFile' and '$modFile' to '$outFile'" if($verbosity >= 5);
    # Print header
    open FILE, "> $outFile";
    print FILE "name\tkmerScore_".$kmerSize."mer\n";
    close FILE;

    # Run KmerInShort
    my $cmd = "$kisPath -file $orfFile -kval $modFile -nb-cores $proc -kmer-size $kmerSize -dont-reverse -step $step 1>> $outFile 2>/dev/null";
    system($cmd);

    # To add the header to the $outFile
    # open FILE, "$outFile";
    # my @cp = <FILE>;
    # close FILE;
    # open FILE, "> $outFile";
    # print FILE "name\tkmerScore_".$kmerSize."mer\n";
    # print FILE @cp;
    # close FILE;
}


# Fusion of the kmerScore files for a list of kmer size
#	@kmerFileList = list of path to each kmer score files
#	$orfSizeFile  = file with the size of the ORF as name\tORFSize
#	$rnaSizeFile  = file with the size of the mRNA as name\tmRNASize
#	$outFile      = file to write the merge of all kmer scores and ORF and mRNA size
#       $nameTmp      = absolute path and prefix for the temporary files
#	$keepTmp      = keep or not temporary files
sub mergeKmerScoreSize
{
    my ($RefKmerFileList, $orfSizeFile, $rnaSizeFile, $outFile, $keepTmp) = @_;
    $RefKmerFileList ||= undef;
    $orfSizeFile     ||= undef;
    $rnaSizeFile     ||= undef;
    $outFile         ||= undef;

    # Check if mendatory arguments have been given
    die "Merging kmer scores and size files: list of kmer scores files is not defined... exiting\n"   if(!defined $RefKmerFileList);
    die "Merging kmer scores and size files: ORF size file is not defined... exiting\n"               if(!defined $orfSizeFile);
    die "Merging kmer scores and size files: mRNA size file is not defined... exiting\n"              if(!defined $rnaSizeFile);
    die "Merging kmer scores and size files: output file for the merging is not defined... exiting\n" if(!defined $outFile);

    # empty
    die "Merging kmer scores and size files: ORF size file '$orfSizeFile' is empty... exiting\n"  unless(-s $orfSizeFile);
    die "Merging kmer scores and size files: mRNA size file '$rnaSizeFile' is empty... exiting\n" unless(-s $rnaSizeFile);


    # Get a hash to stock score values for each sequences and a tab for the header
    my %seq;
    my @head;
    my $flag = 0;
    my $file;

    # Read kmer score files
    foreach $file (@{$RefKmerFileList})
    {
	$flag = 0;
	open FILE, "$file" or die "Error! Cannot access kmer score file '". $file . "': ".$!;
	while(<FILE>)
	{
	    chop;
	    my ($name, $val) = split(/\t/);

	    if($flag != 0)
	    {
		if(!exists $seq{$name})
		{
		    $seq{$name} = [$val];
		}
		else
		{
		    push(@{$seq{$name}}, $val);
		}
	    }
	    else
	    {
		push(@head, $val);
		$flag = 1;
	    }
	}
	close FILE;
	unlink $file unless($keepTmp != 0);
    }

    # Read ORF size
    $flag = 0;
    open FILE, "$orfSizeFile" or die "Error! Cannot access ORF size file '". $orfSizeFile . "': ".$!;
    while(<FILE>)
    {
	chop;
	my($name, $orfSize) = split(/\t/);

	if($flag != 0)
	{
	    if(!exists $seq{$name})
	    {
		die "Error at merge ORF size step! '$name' is not in the kmer score files: ".$!;
	    }

	    push(@{$seq{$name}}, $orfSize);
	}
	else
	{
	    push(@head, $orfSize);
	    $flag = 1;
	}
    }
    close FILE;
    unlink $orfSizeFile unless($keepTmp != 0);

    # Read mRNA size
    $flag = 0;
    open FILE, "$rnaSizeFile" or die "Error! Cannot access mRNA size file '". $rnaSizeFile . "': ".$!;

    while(<FILE>)
    {
	chop;
	my($name, $rnaSize) = split(/\t/);

	if($flag != 0)
	{
	    if(!exists $seq{$name})
	    {
		warn "The sequence id $name didn't have an ORF with a start and stop codon. Skip this sequence...\n";
		next;
		#die "Error at merge mRNA size step! $name is not in the kmer score files: ".$!;
	    }

	    push(@{$seq{$name}}, $rnaSize);
	}
	else
	{
	    push(@head, $rnaSize);
	    $flag = 1;
	}
    }
    close FILE;
    unlink $rnaSizeFile unless($keepTmp != 0);

    # Write the output file
    open FILE, "> $outFile" or die "Error! Cannot access to the output file '". $outFile . "': ".$!;
    my $seqId;

    # Write the header
    print FILE "name\t";
    print FILE join("\t", @head), "\n";

    # Write the values
    foreach $seqId (sort(keys(%seq)))
    {
	print FILE "$seqId\t";
	print FILE join("\t", @{$seq{$seqId}}), "\n";
    }
    close FILE;
}


# Get the size of each element of a multi fasta and write it in a file
#	$inFile  = multifasta file from which the size need to be write
#	$outFile = output file
#	$header  = header to print after the name
sub getSizeFastaFile
{
    my ($inFile, $outFile, $header) = @_;
    $inFile  ||= undef;
    $outFile ||= undef;
    $header  ||= "size";

    # Check if mendatory arguments have been given
    die "Get size: the input file is not defined... exiting\n" if(!defined $inFile);
    die "Get size: the ouput file is not defined... exiting\n" if(!defined $outFile);

    # empty
    die "Get size: input file '$inFile' is empty... exiting\n" unless(-s $inFile);

    # Put all input sequence in array
    my $multiFasta = new Bio::SeqIO(-file  => $inFile);
    my @seq_array;
    my $seq;
    while( $seq = $multiFasta->next_seq() )
    {
	push(@seq_array, $seq);
    }

    # Get length for each sequence and print it
    open FILE, "> $outFile" or die "Error! Cannot access to the output file '". $outFile . "': ".$!;
    my $length;
    my $id;

    print FILE "name\t$header\n";
    foreach my $seq ( @seq_array )
    {
        $length = $seq->length;
        $id     = $seq->id();
	# print "$id\t$length\n";
	print FILE "$id\t$length\n";
    }
    close FILE;
}

# Run the R script RSCRIPT_RF.R to comput the random forest model on learning data and apply it the test data
#	$codLearnFile = file where the kmerscores, mRNA and ORF size are put for learning coding sequences
#	$nonLearnFile = file where the kmerscores, mRNA and ORF size are put for learning non coding sequences
#	$testFile     = file where the kmerscores, mRNA and ORF size are put for test sequences
#	$outFile      = file to write the result of the random forest
#	$thres        = if a valid value is given ([0,1]) then it would be the threshold for the random forest as val>=$thres => coding, if undef then the threshold is obtain by 10-fold cross validation
#	$nTree           = number of trees used in random forest
sub getRunModel
{
    my ($codLearnFile, $nonLearnFile, $testFile, $outFile, $thres, $nTree, $seed) = @_;
    $codLearnFile ||= undef;
    $nonLearnFile ||= undef;
    $testFile     ||= undef;
    $outFile      ||= undef;
    $thres        ||= undef;
    $nTree        ||= 500;

    # Check if mendatory arguments have been given
    die "Running random forest: predictor file for learning coding sequences is not defined... exiting\n"     if(!defined $codLearnFile);
    die "Running random forest: predictor file for learning non coding sequences is not defined... exiting\n" if(!defined $nonLearnFile);
    die "Running random forest: predictor file for testing sequences is not defined... exiting\n"             if(!defined $testFile);
    die "Running random forest: output file is not defined... exiting\n"                                      if(!defined $outFile);

    # emtpy
    die "Running random forest: predictor file for learning coding sequences '$codLearnFile' is empty... exiting\n"     unless (-s $codLearnFile);
    die "Running random forest: predictor file for learning non coding sequences '$nonLearnFile' is empty... exiting\n" unless (-s $nonLearnFile);
    die "Running random forest: predictor file for testing sequences '$testFile' is empty... exiting\n"              unless (-s $testFile);


    my $rprogpath = "";
    my $cmd       = "";
    # First check if $thres is defined
    if(!defined $thres)
    {
	# Get the path to codpot_randomforest.r to run the learning and assignment of the sequences
	$rprogpath = $ENV{'FEELNCPATH'}."/utils/codpot_randomforest.r";
	$cmd       = "";
	# If no threshold given, run codpot_randomforest.r without threshold (6 arguments)
	print "\tThe threshold for the voting in random forest is not defined. Use 10-fold cross-validation to determine the best threshold.\n";
	$cmd = "$rprogpath $codLearnFile $nonLearnFile $testFile $outFile $nTree $seed";
	#warn $cmd;
    }
    # Second check if there is one threshold or two
    else
    {
	my @thresList = split(/,/, $thres);
	if(@thresList==1)
	{
	    # Get the path to codpot_randomforest.r to run the learning and assignment of the sequences
	    $rprogpath = $ENV{'FEELNCPATH'}."/utils/codpot_randomforest.r";
	    $cmd       = "";
	    # If a unique threshold is given, run codpot_randomforest.r with this threshold (7 arguments)
	    print "\tThe threshold for the voting in random forest is '$thres'.\n";
	    $cmd = "$rprogpath $codLearnFile $nonLearnFile $testFile $outFile $nTree $seed $thres";
	}
	else
	{
	    # Get the path to codpot_randomforest.r to run the learning and assignment of the sequences
	    $rprogpath = $ENV{'FEELNCPATH'}."/utils/codpot_randomforest_2thres.r";
	    $cmd       = "";
	    # If two thresholds are given, run codpot_randomforest_2thres.r
	    print "\tTwo specificity thresholds based on mRNA and lncRNA and 10-fold cross-validation: '$thresList[0]' and '$thresList[1]'.\n";
	    $cmd = "$rprogpath $codLearnFile $nonLearnFile $testFile $outFile $nTree $seed $thresList[0] $thresList[1]";
	}
    }
    system($cmd);
}


# Run all the process using ORF from coding sequences, fasta from coding and non coding sequences, ORF and fasta from test sequences, a list of kmer size and a threshold (if defined)
#	$codLearnFile    = fasta file of the learning coding sequences
#	$orfCodLearnFile = fasta file of the ORF learning coding sequences
#	$nonLearnFile    = fasta file of the learning non coding sequences
#	$orfNonLearnFile = fasta file of the ORF learning non coding sequences
#	$testFile        = fasta file of the test sequences
#	$orfTestFile     = fasta file of the ORF test sequences
#	$outFile         = file to write the result of the random forest
#	$kmerListString  = list of size of kmer as '3,6,9' (default value) as string
#	$thres           = the threshold for the random forest, if it is not defined then it is set using a 10-fold cross-validation on learning data
#	$nTree           = number of trees used in random forest
#	$verbosity       = value for level of verbosity
#       $nameTmp         = absolute path and prefix for the temporary files
#	$keepTmp         = keep or not the temporary files
sub runRF
{
    my ($REFcodLearnFile, $REForfCodLearnFile, $REFnonLearnFile, $REForfNonLearnFile, $testFile, $orfTestFile, $outFile, $kmerListString, $thres, $nTree, $outDir, $verbosity, $nameTmp, $keepTmp, $seed, $proc) = @_;
    $kmerListString ||= "3,6,9";
    $verbosity      ||= 0;
    $proc           ||= 1;

    # Get the name of the training files for the random forest
    my $codLearnFile    = $REFcodLearnFile->[1];
    my $orfCodLearnFile = $REForfCodLearnFile->[1];
    my $nonLearnFile    = $REFnonLearnFile->[1];
    my $orfNonLearnFile = $REForfNonLearnFile->[1];

    # 1. Compute the size of each sequence and ORF
    print "1. Compute the size of each sequence and ORF\n";
    my $sizeCodLearnFile    = $nameTmp.".coding_rnaSize.tmp";
    my $sizeOrfCodLearnFile = $nameTmp.".coding_orfSize.tmp";
    my $sizeNonLearnFile    = $nameTmp.".noncoding_rnaSize.tmp";
    my $sizeOrfNonLearnFile = $nameTmp.".noncoding_orfSize.tmp";
    my $sizeTestFile        = $nameTmp.".test_rnaSize.tmp";
    my $sizeOrfTestFile     = $nameTmp.".test_orfSize.tmp";



    # Learning
    ## Coding
    &getSizeFastaFile($codLearnFile,    $sizeCodLearnFile,    "mRNA_size");
    &getSizeFastaFile($orfCodLearnFile, $sizeOrfCodLearnFile, "ORF_size");
    ## Non coding
    &getSizeFastaFile($nonLearnFile,    $sizeNonLearnFile,    "mRNA_size");
    &getSizeFastaFile($orfNonLearnFile, $sizeOrfNonLearnFile, "ORF_size");
    # Test
    &getSizeFastaFile($testFile,    $sizeTestFile,    "mRNA_size");
    &getSizeFastaFile($orfTestFile, $sizeOrfTestFile, "ORF_size");

    # 2. Compute the kmer ratio for each kmer and put the output file name in a list
    print "2. Compute the kmer ratio for each kmer and put the output file name in a list\n";
    my @kmerList = split(/,/, $kmerListString);
    my @kmerRatioFileList;
    my @kmerScoreCodLearnFileList;
    my @kmerScoreNonLearnFileList;
    my @kmerScoreTestFileList;
    my $kmerFile;
    my $kmerSize;
    my $codStep       = 3;
    my $nonStep       = 3;
    my $lenKmerList   = @kmerList;
    my $i             = 0;

    foreach $kmerSize ( @kmerList )
    {
	$kmerFile = $nameTmp.".kmerScoreValues_size$kmerSize.tmp";
	push(@kmerRatioFileList, $kmerFile);

	## VW: modification of the score
	&getKmerRatioSep($REForfCodLearnFile->[0], $REForfNonLearnFile->[0], $kmerFile, $kmerSize, $codStep, $nonStep, $proc, $verbosity, $nameTmp, $keepTmp);
    }

    # 3. Compute the kmer score for each kmer size on learning and test ORF and for each type
    print "3. Compute the kmer score for each kmer size on learning and test ORF\n";
    $kmerFile = "";
    for($i=0; $i<$lenKmerList; $i++)
    {
	print "\t- kmer size: $kmerList[$i]\n";
	# Learning
	## Coding
	$kmerFile = $nameTmp.".coding_sequencesKmer_".$kmerList[$i]."_ScoreValues.tmp";
	push(@kmerScoreCodLearnFileList, $kmerFile);
	scoreORF($orfCodLearnFile, $kmerRatioFileList[$i], $kmerFile, $kmerList[$i], $codStep, $proc, $verbosity);

	## Non coding
	$kmerFile = $nameTmp.".noncoding_sequencesKmer_".$kmerList[$i]."_ScoreValues.tmp";
	push(@kmerScoreNonLearnFileList, $kmerFile);
	scoreORF($orfNonLearnFile, $kmerRatioFileList[$i], $kmerFile, $kmerList[$i], $codStep, $proc, $verbosity);

	# Test
	$kmerFile = $nameTmp.".test_sequencesKmer_".$kmerList[$i]."_ScoreValues.tmp";
	push(@kmerScoreTestFileList, $kmerFile);
	scoreORF($orfTestFile, $kmerRatioFileList[$i], $kmerFile, $kmerList[$i], $codStep, $proc, $verbosity);
    }


    # 4. Merge the score and size files into one file for each type (learning coding and non coding and test)
    print "4. Merge the score and size files into one file for each type\n";
    my $outModCodLearn = $nameTmp.".modelCoding.out";
    my $outModNonLearn = $nameTmp.".modelNonCoding.out";
    my $outModTest     = $nameTmp.".modelTest.out";

    # Learning
    ## Coding
    &mergeKmerScoreSize(\@kmerScoreCodLearnFileList, $sizeOrfCodLearnFile, $sizeCodLearnFile, $outModCodLearn, $keepTmp);
    ## Non coding
    &mergeKmerScoreSize(\@kmerScoreNonLearnFileList, $sizeOrfNonLearnFile, $sizeNonLearnFile, $outModNonLearn, $keepTmp);
    # Test
    &mergeKmerScoreSize(\@kmerScoreTestFileList,     $sizeOrfTestFile,     $sizeTestFile,     $outModTest,     $keepTmp);


    # 5. Make the model on learning sequences and apply it on test sequences
    print "5. Make the model on learning sequences and apply it on test sequences\n";

    &getRunModel($outModCodLearn, $outModNonLearn, $outModTest, $outFile, $thres, $nTree, $seed);

    # Delete temporary files
    if($keepTmp == 0)
    {
	my $file = "";
	unlink $sizeCodLearnFile;
	unlink $sizeOrfCodLearnFile;
	unlink $sizeNonLearnFile;
	unlink $sizeOrfNonLearnFile;
	unlink $sizeTestFile;
	unlink $sizeOrfTestFile;
	foreach $file (@kmerRatioFileList)
	{
	    unlink $file;
	}
	foreach $file (@kmerScoreCodLearnFileList)
	{
	    unlink $file;
	}
	foreach $file (@kmerScoreNonLearnFileList)
	{
	    unlink $file;
	}
	foreach $file (@kmerScoreTestFileList)
	{
	    unlink $file;
	}
	unlink $outModCodLearn;
	unlink $outModNonLearn;
	unlink $outModTest;
    }
}





# With a .GTF/FASTA and a result file from a random forest, write 2 or 3 (if TUCps) .GTF/.FASTA file, each one respectively for coding and non coding genes plus one for sequences with no ORF. Unlink $outTuc and $noOrf if they are empty
#	$testFile = the GTF/FASTA file of the unknown transcripts
#	$reFile   = the output file from the random forest giving transcripts class (-1: TUCp; 0: non coding; 1: coding)
#	$outDir   = ouptut directory
sub rfPredToOut
{
    my($testFile, $rfFile, $outDir, $outName) = @_;
    $testFile ||= undef;
    $rfFile   ||= undef;
    $outDir   ||= undef;
    $outName  ||= basename($testFile);


    # Check if mendatory arguments have been given
    die "Parsing random forest output: GTF file with the new transcripts is not defined... exiting\n" if(!defined $testFile);
    die "Parsing random forest output: random forest output file is not defined... exiting\n"         if(!defined $rfFile);

    # emtpy
    die "Parsing random forest output: GTF file with the new transcripts '$testFile' is empty... exiting\n" unless (-s $testFile);
    die "Parsing random forest output: random forest output file '$rfFile' is empty... exiting\n"           unless (-s $rfFile);


    # Start by reading the result file from the random forest
    open FILE, "$rfFile" or die "Error! Cannot access to the random forest output file '". $rfFile . "': ".$!;
    # Put transcript names on 3 tab, one for coding genes, one for non coding genes and the other for TUCp
    my %nonHas;
    my %codHas;
    my %tucHas;
    my @info;
    my $flag = 0;

    while(<FILE>)
    {
	chop;
	@info = split(/\t/);

	# Check for the first line...
	if($flag==0)
	{
	    $flag = 1;
	    next;
	}

	# Check if the transcript is annotated as non coding (0) or coding (1)
	if($info[-1] == 0)
	{
	    $nonHas{$info[0]} = 0;
	}
	elsif($info[-1] == 1)
	{
	    $codHas{$info[0]} = 1;
	}
	elsif($info[-1] == -1)
	{
	    $tucHas{$info[0]} = -1;
	}
	else
	{
	    die "Not a valid value for coding label for the transcript '$info[0]'. Exit.\n";
	}
    }
    close FILE;


    if(Utils::guess_format($testFile) eq "gtf")
    {
	# Read the GTF file and put the line in the right file depending on which tab the transcript is
	my $outNon = $outDir.$outName.".lncRNA.gtf";
	my $outCod = $outDir.$outName.".mRNA.gtf";
	my $outTuc = $outDir.$outName.".TUCp.gtf";
	my $noOrf  = $outDir.$outName.".noORF.gtf";
	my $line   = "";
	my $name   = "";

	print "Writing the GTF output files\n";
	open FILE,  "$testFile" or die "Error! Cannot access to the GTF input for new transcripts '". $testFile . "': ".$!;
	open LNC, "> $outNon"   or die "Error! Cannot access to the lncRNA GTF output file '". $outNon . "': ".$!;
	open RNA, "> $outCod"   or die "Error! Cannot access to the mRNA GTF output file '". $outCod . "': ".$!;
	open TUC, "> $outTuc"   or die "Error! Cannot access to the TUCp GTF output file '". $outTuc . "': ".$!;
	open NO,  "> $noOrf"    or die "Error! Cannot access to the no ORF GTF output file '". $noOrf . "': ".$!;

	while(<FILE>)
	{
	    $line = $_;
	    $name = ($line =~/.*transcript_id "([^ ]*)";/);
	    $name = $1;

	    # Check if the name is in the non coding or coding array
	    if(exists $nonHas{$name})
	    {
		print LNC $line;
	    }
	    elsif(exists $codHas{$name})
	    {
		print RNA $line;
	    }
	    elsif(exists $tucHas{$name})
	    {
		print TUC $line;
	    }
	    else
	    {
		print NO $line;
	    }
	}
	close FILE;
	close LNC;
	close RNA;
	close TUC;
	close NO;

	# Delete $outTuc and $noOrf if they are empty
	unlink $outTuc if( ! (-s $outTuc) );
	unlink $noOrf  if( ! (-s $noOrf) );
    }
    else # if FASTA format
    {
	# Read the FASTA file and put the line in the right file depending on which tab the transcript is
	my $outNon = $outDir.$outName.".lncRNA.fa";
	my $outCod = $outDir.$outName.".mRNA.fa";
	my $outTuc = $outDir.$outName.".TUCp.fa";
	my $noOrf  = $outDir.$outName.".noORF.fa";

	print "Writing the FASTA output files\n";
	my $multiFasta = new Bio::SeqIO(-file => "$testFile", '-format' => 'Fasta');
	my $lnc        = new Bio::SeqIO(-file => "> $outNon", '-format' => 'Fasta');
	my $rna        = new Bio::SeqIO(-file => "> $outCod", '-format' => 'Fasta');
	my $tuc        = new Bio::SeqIO(-file => "> $outTuc", '-format' => 'Fasta');
	my $noorf      = new Bio::SeqIO(-file => "> $noOrf",  '-format' => 'Fasta');

	my @seqTab;
	my $seq;
	while( $seq = $multiFasta->next_seq() )
	{
	    # Check if the name is in the non coding or coding array
	    if(exists $nonHas{$seq->id()})
	    {
		$lnc->write_seq($seq);
	    }
	    elsif(exists $codHas{$seq->id()})
	    {
		$rna->write_seq($seq);
	    }
	    elsif(exists $tucHas{$seq->id()})
	    {
		$tuc->write_seq($seq);
	    }
	    else
	    {
		$noorf->write_seq($seq);
	    }
	}

	# Delete $outTuc and $noOrf if they are empty
	unlink $outTuc if( ! (-s $outTuc) );
	unlink $noOrf  if( ! (-s $noOrf) );
    }
}


1;



__END__

=encoding UTF-8

=head2 getKmerRatio

Run KmerInShort on learning set and generate output model (logRatio values for each kmer)

=over

=item $codFile

fasta file with the ORF of coding genes

=item $nonFile

fasta file with sequence of non coding genes

=item $outFile

file where the logRatio of kmer frequency need to be written

=item $kmerSize

size of the kmer

=item $codStep

step for the counting of kmer for the ORF coding genes

=item $nonStep

step for the counting of kmer for the non coding genes

=item $proc

number of proc to be use for KmerInShort

=item $verbosity

define the verbosity of the program

=item $nameTmp

absolute path and prefix for the temporary files

=item $keepTmp

if the temporary files need to be kept (0: delete; !0: keep)

=back

Return value: Return the array of the log ratio value

##############################################################################

=head2 getKmerRatioSep

Run KmerInShort on learning set and generate output model (ratio values for each kmer depending on the occurancy of each kmer between coding and non coding files)

=over

=item $codFile

fasta file with the ORF of coding genes

=item $nonFile

fasta file with sequence of non coding genes

=item $outFile

file where the logRatio of kmer frequency need to be written

=item $kmerSize

size of the kmer

=item $codStep

step for the counting of kmer for the ORF coding genes

=item $nonStep

step for the counting of kmer for the non coding genes

=item $proc

number of proc to be use for KmerInShort

=item $verbosity

define the verbosity of the program

=item $nameTmp

absolute path and prefix for the temporary files

=item $keepTmp

if the temporary files need to be kept (0: delete; !0: keep)

=back

Return value: Return the array of the log ratio value

##############################################################################

=head2 scoreORF

Temporary function to compute the scoring function of a multifasta file of ORF

=over

=item $orfFile

multi fasta file with the ORF of the genes to be scored

=item $modFile

file with the log ratio score compute on learning files

=item $outFile

file where the logRatio of kmer frequency need to be written

=item $kmerSize

size of the kmer

=item $step

step for the counting of kmer

=item $proc

number of proc to be use for KmerInShort

=item $verbosity

define the verbosity of the program

=back

##############################################################################

=head2 getSizeFastaFile

Get the size of each element of a multi fasta and write it in a file

=over

=item $inFile

multifasta file from which the size need to be write

=item $outFile

output file

=item $header

header to print after the name

=back

##############################################################################

=head2 mergeKmerScoreSize

Fusion of the kmerScore files for a list of kmer size

=over

=item \@kmerFileList

list of path to each kmer score files

=item $orfSizeFile

file with the size of the ORF as name\tORFSize

=item $rnaSizeFile

file with the mRNA size as name\tmRNASize

=item $outFile

file to write the merge of all kmer scores and ORF/mRNA size

=item $nameTmp

absolute path and prefix for the temporary files

=item $keepTmp

if the temporary files need to be kept (0: delete; !0: keep)

=back

#############################################################################

=head2 getRunModel

Run the R script BLABLA to comput the random forest model on learning data and apply it the test data

=over

=item $codLearnFile

file where the kmerscores, mRNA and ORF size are put for learning coding sequences

=item $nonLearnFile

file where the kmerscores, mRNA and ORF size are put for learning non coding sequences

=item $testFile

file where the kmerscores, mRNA and ORF size are put for test sequences

=item $thres

if it is given, the threshold for random forest (>$thres = coding), if undef then the threshold is obtain by 10-fold cross validation and if it is a string with two float separated by a ',', then use two specificity thresholds


=item $nTree

the number of trees used in random forest

=item $seed

seed to fixe random fonction during the shuffling of the input matrix and for the random forest

=back

##############################################################################

=head2 runRF

Run all the process using ORF from coding sequences, fasta from coding and non coding sequences, ORF and fasta from test sequences, a list of kmer size and a threshold (if defined)

=over

=item $codLearnFile

fasta file of the learning coding sequences

=item $orfCodLearnFile

fasta file of the ORF learning coding sequences

=item $nonLearnFile

fasta file of the learning non coding sequences

=item $orfNonLearnFile

fasta file of the ORF learning non coding sequences

=item $testFile

fasta file of the test sequences

=item $orfTestFile

fasta file of the ORF test sequences

=item $outFile

file to write the result of the random forest

=item $kmerListString

list of size of kmer as '3,6,9' (default value) as a string

=item $thres

the threshold for the random forest, if it is not defined then it is set using a 10-fold cross-validation on learning data and if it is a string with two float separated by a ',', then use two specificity thresholds

=item $nTree

the number of trees used in random forest

=item $outDir

output directory

=item $verbosity

define the verbosity of the program

=item $nameTmp

absolute path and prefix for the temporary files

=item $keepTmp

if the temporary files need to be kept (0: delete; !0: keep)

=item $seed

seed to fixe random fonction during the shuffling of the input matrix and for the random forest

=back

##############################################################################

=head2 rfPredToOut

With a .gtf and a result file from a random forest, write 2 .gtf file, each one respectively for coding and non coding genes

=over

=item $testGTF

the GTF file of the unknown transcripts

=item $reFile

the output file from the random forest giving ranscripts class (0: non coding; 1: coding)

=item $outDir

output directory

=item $outName

output name

=back

=cut
