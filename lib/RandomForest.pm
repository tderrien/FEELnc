package RandomForest;

$VERSION = v0.0.1;

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;

use Utils;
use ExtractFromFeature; # getKeyFromFeature line 932
use Bio::DB::Fasta;

=head1 RandomForest.pm

=over

=item .

getKmerRatio: Get one kmer ratio from fasta training file

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

=back

=cut


# Run minidsk on training set and generate output model (logRatio values for each kmer)
#	$codFile = fasta file with the ORF of coding genes
#	$nonFile = fasta file with sequence of non coding genes
#	$outFile = file where the logRatio of kmer frequency need to be written
#	$kmerSize = size of the kmer
#	$codStep = step for the counting of kmer for the ORF coding genes
#	$nonStep = step for the counting of kmer for the non coding genes
#	$proc = number of proc to be use for minidsk
#	Return value:
#		Return the array of the log ratio value
sub getKmerRatio
{
    my($codFile, $nonFile, $outFile, $kmerSize, $codStep, $nonStep, $proc) = @_;
    $codFile ||= undef;
    $nonFile ||= undef;
    $outFile ||= undef;
    $codStep ||= 3;
    $nonStep ||= 1;
    $proc    ||= 1;    
   
    # Check if mendatory arguments have been given
    die "Bulding model: ORF coding genes training file not defined... exiting\n" if (!defined $codFile);
    die "Bulding model: lncRNA training file not defined... exiting\n"           if (!defined $nonFile);
    die "Bulding model: output file for kmer model is missing... exiting\n"      if (!defined $outFile);
    
    # emtpy
    die "Bulding model: ORF coding genes training file '$codFile' is empty... exiting\n" unless (-s $codFile);
    die "Bulding model: lncRNA training file '$nonFile' is empty... exiting\n"           unless (-s $nonFile);
    
    # Path to minidsk
    my $minidskPath = Utils::pathProg("minidsk");
    # Temporary files to put kmer counting for ORF coding genes and non coding genes
    my $codOut = "/tmp/cod_".int(rand(1000)).".tmp"
    my $nonOut = "/tmp/non_".int(rand(1000)).".tmp"
    
    # Run minidsk on ORF for coding genes
    print "Running minidsk on $codFile:\n";
    my $cmd = "$minidskPath -file $codFile -nb-cores $proc -kmer-size $kmerSize -out $codOut -dont-reverse -step $codStep 1>/dev/null 2>/dev/null";
    system($cmd);

    # Run minidsk on non coding genes
    print "Running minidsk on $nonFile:\n";
    my $cmd = "$minidskPath -file $nonFile -nb-cores $proc -kmer-size $kmerSize -out $nonOut -dont-reverse -step $nonStep 1>/dev/null 2>/dev/null";
    system($cmd);

    # Read the two kmer files and put value in a table and get the total number of kmer to comput frequency
    my @kmerTab;
    my @codVal;
    my $codTot = 0;
    my @nonVal;
    my $nonTot = 0;
    my $i;
    
    # Read the kmer counting for ORF on coding genes
    print "Read kmer counting for coding genes of size $kmerSize in $codOut\n";
    # Open file
    open FILE, "$codOut" or die "Error! Cannot open kmerFile ". $codOut . ": ".$!;
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
    print "Read kmer counting for non coding genes of size $kmerSize in $nonOut\n";
    # Open file
    open FILE, "$nonOut" or die "Error! Cannot open kmerFile ". $nonOut . ": ".$!;
    $i = 0;
    while(<FILE>)
    {
	chop;
	my ($kmer, $val) = split(/\t/);

	die "Error: kmer order is not the same between $codOut and $nonOut.\nExit." if($kmerTab[$i] ne $kmer);

	#$nonCod[$i] = $val;
	push(@nonCod, int($val));
	$nonTot = $nonTot + int($val);
	$i = $i+1;
    }
    close FILE;


    # Write the output file with the log ratio directly
    print "Write the log ratio of kmer bewteen coding and non coding genes frequency for a size of kmer of $kmerSize\n";
    open FILE, "> $outFile" or die "Error! Cannot access output file ". $outFile . ": ".$!;
    my $nbKmer = @kmerTab;
    my $log    = 0;
    my @logTab;
    
    # Print the file header
    print FILE "kmer\t$kmerSize_logRatio\n";
    for($i=0; $i<$nbKmer; $i++)
    {
	# logratio = 0                      -- if the kmer is not found in any gene files
	if($codVal[$i]==0 && $nonVal[$i]==0)
	{
	    $log = 0;
	}
	# logratio = 1                      -- if the kmer is found in coding ORFs but not in non coding genes
	elsif($codVal[$i]>0 && $nonVal[$i]==0)
	{
	    $log = 1;
	}
	# logratio = -1                     -- if the kmer is not found in coding ORFs but found in non coding genes
	elsif($codVal[$i]==0 && $nonVal[$i]>0)
	{
	    $log = -1;
	}
	# logratio = log( codFreq/nonFreq ) -- if the kmer is found on the two files
	else
	{
	    $log = log( ($codVal[$i]/$codTot) / ($nonVal[$i]/$nonTot) );
	}

	$logTab[$i] = $log;
	print FILE "$kmerTab[i]\t$log\n";
    }
    close FILE;

    # Delete the temporary files
    unlink $codOut;
    unlink $nonOut;

    return(@logTab);
}



# Temporary function to compute the scoring function of a multifasta file of ORF
#	$orfFile  = multi fasta file with the ORF of the genes to be scored
#	$modFile  = file with the log ratio score compute on learning files
#	$outFile  = file where the logRatio of kmer frequency need to be written 
#	$kmerSize = size of the kmer
#	$step     = step for the counting of kmer
#	$proc     = number of proc to be use for minidsk
sub scoreORF
{
    my($orfFile, $modFile, $outFile, $kmerSize, $step, $proc) = @_;
    $orfFile  ||= undef;
    $modFile  ||= undef;
    $outFile  ||= undef;
    $kmerSize ||= 6;
    $step     ||= 3;
    $proc     ||= 1;

    # Check if mendatory arguments have been given
    die "Scoring ORF file: ORF file for test sequences is not defined... exiting\n"                        if(!defined $orfFile);
    die "Scoring ORF file: model of log ratio score file is not defined... exiting\n"                      if(!defined $modFile);
    die "Scoring ORF file: output file to write kmer score for test sequences is not defined... exiting\n" if(!defined $outFile);

    # emtpy
    die "Scoring ORF file: ORF file for test sequences '$orfFile' is empty... exiting\n"   unless (-s $orfFile);
    die "Scoring ORF file: model of log ratio score file '$modFile' is empty... exiting\n" unless (-s $modFile);

    # Path to minidsk
    my $minidskPath = Utils::pathProg("minidsk");

    # Parse ORF file and put ORF sequence in an array
    print "Read ORF file $orfFile\n";
    my $multiFasta = new Bio::SeqIO(-file  => $orfFile);
    my @seqTab;
    my $seq;
    while( $seq = $multiFasta->next_seq() )
    {
	push(@seq_array,$seq);
    }
    # Fermer le fichier ou autres ?
    

    # Read the model file
    print "Read the log ratio file $modFile\n";
    my @kmerTab;
    my @logTab;
    open FILE, "$modFile" or die "Error! Cannot open the model file ". $modFile . ": ".$!;
    while(<FILE>)
    {
	chop;
	my ($kmer, $log) = split(/\t/);
	push(@kmerTab, $kmer);
	push(@logTab, float($log));
    }
    close FILE;

    my $tmpFile    = "/tmp/orf_".int(rand(1000)).".tmp";
    my $tmpFileOut = "/tmp/out_".int(rand(1000)).".tmp";
    my $id         = "";
    my $length     = 0;
    my $logSum     = 0;
    my $totKmer    = 0;
    my $i          = 0;
    open FILEOUT, "> $outFile" or die "Error! Cannot access output file ". $outFile . ": ".$!;

    # Print the header of the output file
    print FILEOUT "name\t$kmerSize_score\n";
    # For each ORF, print the sequence on a temporary file $tmpFile and run minidsk on this sequence and print the result in $tmpFileOut
    foreach $seq (@seqTab)
    {
	$id     = $seq->id();
	$length = $seq->length();

	# Check if there is a sequence
	if($length == 0)
	{
                warn "Warning: no sequence for ", $id, " (certainly because the scaffold of the transcript does not exists and the sequence has not been extracted correctly)";
                next;
	}

	# Write the ORF sequence on $tmpFile
	my $seqout  = Bio::SeqIO->new(-format => 'fasta', -file => '> '.$tmpFile, -alphabet =>'dna', ); # la ',' à la fin chelou
	my $new_seq = Bio::Seq->new(-id => $id, -seq => $seq);

	# Run minidsk on this sequence
	my $cmd = "$minidskPath -file $tmpFile -nb-cores $proc -kmer-size $kmerSize -out $tmpFileOut -dont-reverse -step $step 1>/dev/null 2>/dev/null";
	system($cmd);

	# Read the minidsk output (FILEMINI) and write it on the output file (FILEOUT)
	open FILEMINI, "$tmpFileOut" or die "Error! Cannot access the temporary output minidsk file ". $tmpFileOut . ": ".$!;
	while(<FILEMINI>)
	{
	    chop;
	    my ($kmer, $nbr) = split(/\t/);

	    die "Error: kmer order is not the same between $modFile and $tmpFileOut.\nExit." if($kmerTab[$i] ne $kmer);
	    
	    $logSum  = $logSum + ( int($nbr)*$logTab[$i] );
	    $totKmer = $totKmer + int($nbr);
	}
	close FILEMINI;

	print FILEOUT "$id\t$logSum";
    }
    close FILEOUT;


    # Delete temporary files
    unlink $tmpFile;
    unlink $tmpFileOut;
}




# Fusion of the kmerScore files for a list of kmer size
#	@kmerFileList = list of path to each kmer score files
#	$orfSizeFile  = file with the size of the ORF as name\tORFSize
#	$rnaSizeFile  = file with the size of the mRNA as name\tmRNASize
#	$outFile      = file to write the merge of all kmer scores and ORF and mRNA size
sub mergeKmerScoreSize
{
    my ($RefKmerFileList, $orfSizeFile, $rnaSizeFile, $outFile) = @_;
    $RefKmerFileList ||= undef;
    $orfSizeFile     ||= undef;
    $rnaSizeFile     ||= undef;
    $outFile         ||= undef;

    # Check if mendatory arguments have been given
    die "Merging kmer scores and size files: list of kmer scores files is not defined... exiting\n"   if(!defined $RefKmerFileList);
    die "Merging kmer scores and size files: ORF size file is not defined... exiting\n"               if(!defined $orfSizeFile);
    die "Merging kmer scores and size files: mRNA size file is not defined... exiting\n"              if(!defined $rnaSizeFile);
    die "Merging kmer scores and size files: output file for the merging is not defined... exiting\n" if(!defined $orfRnaSizeFile);

    # empty
    die "Merging kmer scores and size files: ORF size file is empty... exiting\n" unless(-s $orfSizeFile);
    die "Merging kmer scores and size files: mRNA size file is empty... exiting\n" unless(-s $rnaSizeFile);


    # Get a hash to stock score values for each sequences and a tab for the header
    my %seq;
    my @head;
    my $flag = 0;
    
    # Read kmer score files
    foreach $file in (@{$RefKmerFileList})
    {
	$flag = 0;
	open FILE, "$file" or die "Error! Cannot access kmer score file ". $file . ": ".$!;
	while(<FILE>)
	{
	    chop;
	    my ($name ,$val) = split(/\t/);

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
	unlink $file;
    }

    # Read ORF size
    $flag = 0;
    open FILE, "$orfSizeFile" or die "Error! Cannot access ORF size file ". $orfSizeFile . ": ".$!;
    while(<FILE>)
    {
	chop;
	my($name, $orfSize) = split(/\t/);

	if($flag != 0)
	{
	    if(!exists $seq{$name})
	    {
		die "Error! $name is not in the kmer score files: ".$!;
	    }

	    push(@{$seq{$name}}, $orfSize);
	}
	else
	{
	    push(@head, $val);
	    $flag = 0;
	}
    }
    close FILE;

    # Read mRNA size
    $flag = 0;
    open FILE, "$rnaSizeFile" or die "Error! Cannot access mRNA size file ". $rnaSizeFile . ": ".$!;
    while(<FILE>)
    {
	chop;
	my($name, $rnaSize) = split(/\t/);

	if($flag != 0)
	{
	    if(!exists $seq{$name})
	    {
		die "Error! $name is not in the kmer score files: ".$!;
	    }

	    push(@{$seq{$name}}, $rnaSize);
	}
	else
	{
	    push(@head, $val);
	    $flag = 0;
	}
    }
    close FILE;

    # Write the output file
    open FILE, "> $outFile" or die "Error! Cannot access to the output file ". $outFile . ": ".$!;
    my $seqId;

    # Write the header
    print FILE join("\t", @head), "\n";

    # Write the values
    foreach $seqId (keys(%seq))
    {
	print FILE "$seqId\t";
	print join("\t", @seq[$seqId]), "\n";
    }
    close FILE;
}


# Get the size of each element of a multi fasta and write it in a file
#	$inFile  = multifasta file from which the size need to be write
#	$outFile = output file
#	$header  = header to print after the name
sub getSizeFastaFile
{
    my ($inFile, $outFile) = @_;
    $inFile  ||= undef;
    $outFile ||= undef;
    $header  ||= "size";

    # Check if mendatory arguments have been given
    die "Get size: the input file is not defined... exiting\n" if(!defined $inFile);
    die "Get size: the ouput file is not defined... exiting\n" if(!defined $outFile);

    # empty
    die "Get size: input file $inFile is empty... exiting\n" unless(-s $inFile);
 
    # Put all input sequence in array
    my @seq_array;
    my $seq;
    while( $seq = $in->next_seq() )
    {
	push(@seq_array,$seq);
    }

    # Get length for each sequence and print it
    open FILE, "> $outFile" or die "Error! Cannot access to the output file ". $outFile . ": ".$!;
    my $length;
    my $id;

    print FILE "name\t$header";
    foreach my $seq ( @seq_array ) 
    {
        $length = $seq->length;
        $id     = $seq->id();
	print FILE "$id\t$length";
    }
    close FILE;
}

# Run the R script RSCRIPT_RF.R to comput the random forest model on learning data and apply it the test data
#	$codLearnFile = file where the kmerscores, mRNA and ORF size are put for learning coding sequences
#	$nonLearnFile = file where the kmerscores, mRNA and ORF size are put for learning non coding sequences
#	$testFile     = file where the kmerscores, mRNA and ORF size are put for test sequences
#	$outFile      = file to write the result of the random forest
#	$thres        = if a valid value is given ([0,1]) then it would be the threshold for the random forest as val>=$thres => coding, if undef then the threshold is obtain by 10-fold cross validation
sub getRunModel
{
    my ($codLearnFile, $nonLearnFile, $testFile, $outFile, $thres) = @_;
    $codLearnFile ||= undef;
    $nonLearnFile ||= undef;
    $testFile     ||= undef;
    $outFile      ||= undef;
    $thres        ||= undef;

    # Check if mendatory arguments have been given
    die "Running random forest: predictor file for learning coding sequences is not defined... exiting\n"     if(!defined $codLearnFile);
    die "Running random forest: predictor file for learning non coding sequences is not defined... exiting\n" if(!defined $nonLearnFile);
    die "Running random forest: predictor file for testing sequences is not defined... exiting\n"             if(!defined $testFile);
    die "Running random forest: output file is not defined... exiting\n"                                      if(!defined $outFile);

    # emtpy
    die "Running random forest: predictor file for learning coding sequences '$learnFile' is empty... exiting\n"     unless (-s $codLearnFile);
    die "Running random forest: predictor file for learning non coding sequences '$learnFile' is empty... exiting\n" unless (-s $nonLearnFile);
    die "Running random forest: predictor file for testing sequences '$testFile' is empty... exiting\n"              unless (-s $testFile);


    # Get the path to RSCRIPT_RF.R to run the learning and assignment of the sequences
    print "Running RSCRIPT_RF.R\n";
    my $rprogpath = $ENV{'FEELNCPATH'}."/bin/RSCRIPT_RF.R";
    # Not a valid value for the threshold
    # It is done in FEELnc_codpot.pl
    # No threshold given
    if(!defined $thres)
    {
	print "The threshold for the voting in random forest is not defined. Use 10-fold cross-validation to determine the best threshold.\n";
	my $cmd = "$rprogpath $codLearnFile $nonLearnFile $testFile $outFile";
    }
    else
    {
	print "The threshold for the voting in random forest is $thres.\n";
	my $cmd = "$rprogpath $codLearnFile $nonLearnFile $testFile $outFile $thres";
    }
    #system($cmd);
    print "LA COMMANDE : $cmd\n";
}


# Run all the process using ORF from coding sequences, fasta from coding and non coding sequences, ORF and fasta from test sequences, a list of kmer size and a threshold (if defined)
#	$codLearnFile    = fasta file of the learning coding sequences
#	$orfCodLearnFile = fasta file of the ORF learning coding sequences
#	$nonLearnFile    = fasta file of the learning non coding sequences
#	$orfNonLearnFile = fasta file of the ORF learning non coding sequences
#	$testFile        = fasta file of the test sequences
#	$orfTestFile     = fasta file of the ORF test sequences
#	$outFile         = file to write the result of the random forest
#	$kmerListString  = list of size of kmer as '2,3,4,5,6' (default value) as string
#	$thres           = the threshold for the random forest, if it is not defined then it is set using a 10-fold cross-validation on learning data
sub runRF
{
    my ($codLearnFile, $orfCodLearnFile, $nonLearnFile, $orfNonLearnFile, $testFile, $orfTestFile, $outFile, $kmerListString, $thres) = @_;
    $kmerListString ||= "2,3,4,5,6";
    
    # 1. Compute the size of each sequence and ORF
    print "1. Compute the size of each sequence and ORF\n";
    my $sizeCodLearnFile    = "/tmp/".basename($codLearnFile)."_".int(rand(1000)).".tmp";
    my $sizeOrfCodLearnFile = "/tmp/".basename($orfCodLearnFile)."_".int(rand(1000)).".tmp";
    my $sizeNonLearnFile    = "/tmp/".basename($nonLearnFile)."_".int(rand(1000)).".tmp";
    my $sizeOrfNonLearnFile = "/tmp/".basename($orfNonLearnFile)."_".int(rand(1000)).".tmp";
    my $sizeTestFile        = "/tmp/".basename($testFile)."_".int(rand(1000)).".tmp";
    my $sizeOrfTestFile     = "/tmp/".basename($orfTestFile)."_".int(rand(1000)).".tmp";

    # Learning
    ## Coding
    &getSizeFastaFile($codLearnFile,    $sizeCodLearnFile);
    &getSizeFastaFile($orfCodLearnFile, $sizeOrfCodLearnFile);
    ## Non coding
    &getSizeFastaFile($nonLearnFile,    $sizeNonLearnFile);
    &getSizeFastaFile($orfNonLearnFile, $sizeOrfNonLearnFile);
    # Test
    &getSizeFastaFile($testFile,    $sizeTestFile);
    &getSizeFastaFile($orfTestFile, $sizeOrfTestFile);

    # 2. Compute the kmer ratio for each kmer and put the output file name in a list
    print "2. Compute the kmer ratio for each kmer and put the output file name in a list\n";
    my @kmerList = split(/,/, $kmerListString);
    my @kmerRatioFileList;
    my @kmerScoreCodLearnFileList;
    my @kmerScoreNonLearnFileList;
    my @kmerScoreTestFileList;
    my $kmerFile;
    my $kmerSize;
    my $codStep     = 3;
    my $nonStep     = 1;
    my $proc        = 1;
    my $lenKmerList = @kmerList;
    my $i           = 0;
    
    foreach $kmerSize ( @kmerList )
    {
	$kmerFile = "/tmp/".basename($orfCodLearnFile)."_".int(rand(10000))."_kmerRatio.tmp";
	push(@kmerRatioFileList, $kmerFile);

	&getKmerRatio($orfCodLearnFile, $nonLearnFile, $kmerFile, $kmerSize, $codStep, $nonStep, $proc);
    }
    
    # 3. Compute the kmer score for each kmer size on learning and test ORF and for each type
    print "3. Compute the kmer score for each kmer size on learning and test ORF\n";
    $kmerFile = "";
    for($i=0; i<$lenKmerList, $i++)
    {
	# Learning
	## Coding
	$kmerFile = "/tmp/".basename($orfCodLearnFile)."_".int(rand(10000))."_kmerScoreCodLearn.tmp";
	push(@kmerScoreCodLearnFileList, $kmerFile);
	scoreORF($orfCodLearnFile, $kmerRatioFileList[$i], $kmerFile, $kmerList[$i], $codStep, $proc);

	## Non coding
	$kmerFile = "/tmp/".basename($orfNonLearnFile)."_".int(rand(10000))."_kmerScoreNonLearn.tmp";
	push(@kmerScoreNonLearnFileList, $kmerFile);
	scoreORF($orfNonLearnFile, $kmerRatioFileList[$i], $kmerFile, $kmerList[$i], $codStep, $proc);

	# Test
	$kmerFile = "/tmp/".basename($orfTestFile)."_".int(rand(10000))."_kmerScoreTest.tmp";
	push(@kmerScoreTestFileList, $kmerFile);
	scoreORF($orfTestFile, $kmerRatioFileList[$i], $kmerFile, $kmerList[$i], $codStep, $proc);
    }    


    # 4. Merge the score and size files into one file for each type (learning coding and non coding and test)
    print "4. Merge the score and size files into one file for each type\n";
    my $outModCodLearn = basename($outFile).".modelCoding.out";
    my $outModNonLearn = basename($outFile).".modelNonCoding.out";
    my $outModTest     = basename($outFile).".modelTest.out";
    
    # Learning
    ## Coding
    &mergeKmerScoreSize(\@kmerScoreCodLearnFileList, $sizeOrfCodLearnFile, $sizeCodLearnFile, $outModCodLearn);
    ## Non coding
    &mergeKmerScoreSize(\@kmerScoreNonLearnFileList, $sizeOrfNonLearnFile, $sizeNonLearnFile, $outModNonLearn);
    # Test
    &mergeKmerScoreSize(\@kmerScoreTestFileList,     $sizeOrfTestFile,     $sizeTestFile,     $outModTest);


    # 5. Make the model on learning sequences and apply it on test sequences
    print "5. Make the model on learning sequences and apply it on test sequences\n";

    &getRunModel($outModCodLearn, $outModNonLearn, $outModTest, $outFile);
}



1;



__END__

=head2 getKmerRatio

Run minidsk on learning set and generate output model (logRatio values for each kmer)

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

number of proc to be use for minidsk

=back

Return value: Return the array of the log ratio value

=cut

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

number of proc to be use for minidsk

=back

=cut

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

=back

=cut

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

if it is given, the threshold for random forest (>$thres = coding), if undef then the threshold is obtain by 10-fold cross validation

=back

=cut

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

=item $testFile

fasta file of the test sequences

=item $orfTestFile

fasta file of the ORF test sequences

=item $outFile

file to write the result of the random forest

=item $kmerListString

list of size of kmer as '2,3,4,5,6' (default value) as a string

=item $thres

the threshold for the random forest, if it is not defined then it is set using a 10-fold cross-validation on learning data

=back

=cut
