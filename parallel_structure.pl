#! /usr/bin/perl

# the next line gives the location of my personal modules.  I did this on the cluster to use the File::Copy::Recursive module
use lib qw(/home/mussmann/perl5/lib/perl5);

use warnings;
use strict;
use Cwd;
use Getopt::Std;
use File::Path;
use File::Copy qw( move );
use File::Copy::Recursive qw( dirmove );

# kill program and print help if no command line arguments were given
if( scalar( @ARGV ) == 0 ){
  &help;
  die "Exiting program because no command line options were used.\n\n";
}

# take command line arguments
my %opts;
getopts( 'aA:b:B:c:CdDeEf:Fg:G:hHi:Ijk:K:l:Lm:MnNoOpPqQr:RsStuUvwW:x:X:yYzZ:', \%opts );

# if -h flag is used, or if no command line arguments were specified, kill program and print help
if( $opts{h} ){
  &help;
  die "Exiting program because help flag was used.\n\n";
}

# parse the command line options and return values
my( $data, $startk, $endk, $burn, $gen, $krun, $numind, $numloci, $miss, $excols, $ploidy, $alpha, $lambda, $admix, $row, $uid, $pop, $flag, $loc, $pheno, $rec, $map, $phaseinfo, $phased, $mname, $markov, $link, $usepop, $locprior, $fst, $infalpha, $popalphas, $inflam, $popspeclam, $randomize, $evanno, $structure, $pops, $method, $clumpp, $greed, $distruct, $bottomlabels, $blab, $toplabels, $tlab, $clumppall, $sortq ) = &parsecom( \%opts );

# test if parallel is installed and executable by user if p flag is used
# also set randomize to 0 so random number seeds can be fed to Structure
# if GNU Parallel is not installed, kill program and print error message
if( $opts{p} ){
  my $program = "parallel";
  my( $parpres ) = &progtest( $program );
  if( $parpres == 0 ){
    die "ERROR:\tInstall ", $program, " before proceeding\n\tOr run without -p option\n\n";
  }
  $randomize = 0;
}

# Names of directories to be used and output files to be written by the program
my $wd = cwd(); # current working directory
my $harvestdir = $wd . "/harvester"; # directory that contains StructureHarvester.py, CLUMPP, and DISTRUCT output
my $logdir = $wd . "/log"; # directory that contains log files from Structure
my $resultsdir = $wd . "/results"; # directory that contains results output from Structure
my @dirpresent = ( $logdir, $resultsdir, $harvestdir );
my $mainout = "mainparams"; # mainparams file for Structure
my $extraout = "extraparams"; # extraparams file for Structure
my $strun = "structurerun.pl"; # perl script to run Structure
my $com = "structurecom.txt"; # text file to hold Structure commands if Structure is run using GNU Parallel
my @filepresent = ($mainout, $extraout, $strun, $com );
my $clusters; # will hold the number of clusters found by Evanno method (StructureHarvester.py)
my( $oldrun, $oldrundir ) = &getdaytime( $wd ); # will hold the name ($oldrun) and path ($oldrundir) if needed to move old Structure runs to new folder.
my $structuredatafile = $wd . "/" . $data; #path to the structure data file


# if -S flag is used prepare parameter files and scripts, then run structure
if( $structure == 1 ){
  
  # see if structure is installed
  my $program = "structure";
  my( $structpres ) = &progtest( $program );
  if( $structpres == 0 ){
    die "ERROR:\tInstall ", $program, " before proceeding\n\tOr run without -S option\n\n";
  }

  # If the -w flag is used, overwrite current output files
  # If the -j flag is used, preserve old output in a new folder
  # If neither flag is used, prompt what to do
  if( $opts{w} ){
    print "\nAll existing files are being overwritten\n";
  }elsif( $opts{j} ){
    mkpath( [$oldrun], 1, 0755 ) or die "Couldn't make directory to hold previous run in $wd: $!";
    my $targetdir = &copyoldfolder( \@dirpresent, $oldrundir );
    &copyoldfiles( \@filepresent );
  }else{
    mkpath( [$oldrun], 1, 0755 ) or die "Couldn't make directory to hold previous run in $wd: $!";
    &promptoverwrite( \@filepresent, \@dirpresent, $wd, $oldrundir );
    # remove $oldrundir directory if no files were chosen to be saved
    rmdir $oldrundir;
  }
  
  # Write the mainparams file
  &main( $burn, $data, $endk, $excols, $flag, $gen, $loc, $mainout, $map, $markov, $miss, $mname, $numind, $numloci, $phased, $phaseinfo, $pheno, $ploidy, $pop, $rec, $row, $uid );
  
  # Write the extraparams file
  &extra( $admix, $alpha, $extraout, $fst, $infalpha, $inflam, $lambda, $link, $locprior, $popalphas, $popspeclam, $randomize, $usepop );
  
  # Write file to run structure
  if( $opts{p} ){
    &parstruct( $data, $com, $endk, $startk, $strun );
    &strcommand( $com, $data, $endk, $krun, $startk );
  }else{
    &struct( $data, $endk, $krun, $startk, $strun );
  }

  # Make structurerun.pl executable
  chmod( 0755, $strun ) or die "ERROR: Couldn't change permissions for $strun: $!\n";

  # get directory path to structurerun.pl script for the next system call
  my $structurerun = $wd . "/" . $strun;

  # run structure
  system( $structurerun ) == 0 or die "ERROR: System call to $strun failed: $?";
}

# run structure harvester if -E flag is used
if( $evanno == 1 ){
  # see if structureHarvester.py is installed
  my( $program ) = "structureHarvester.py";
  my( $harvpres ) = &progtest( $program );
  if( $harvpres == 0 ){
    die "ERROR:\tInstall ", $program, " before proceeding\n\tOr run without -E option\n\tFile ", $data, "_harvester-upload.zip is in the harvester directory.\n\n";
  }
  # run structure harvester
  system( "structureHarvester.py --dir results/ --out harvester/ --evanno --clumpp" ) == 0 or die "ERROR: System call to structureHarvester.py failed: $?";
}

# run clumpp if -C flag is used
# get best number of clusters found by structureHarvester.py
if( $clumpp == 1){
  # see if clumpp is installed
  my $program = "clumpp";
  my( $clumpppres ) = &progtest( $program );
  if( $clumpppres == 0 ){
    die "ERROR:\tInstall ", $program, " before proceeding\n\tOr run without -C option\n\n";
  }
  
  my $efile = $wd . "/harvester/evanno.txt";
  $clusters = &evannoparse( $efile );
  
  # Variables for creating files
  my $paramfile = "paramfile";
  my $weight = 1;
  my $similarity = 2;
  my( $inds, $outfile, $miscfile, $indfile, $popfile);
  
  # change directories to where clumpp files are stored
  chdir $harvestdir;
  
  if( $clumppall == 1 ){
    my $clumppcom = "clumppcommands.txt";
    open( CLUMPPCOM, '>', $clumppcom );
    # write paramfile and run clumpp for indfile and popfile
    for( my $i = $startk; $i < $endk+1; $i++ ){
      if( $i < 5 ){
	$method = 2;
      }else{
	$method = 3;
      }
      $indfile = "K$i.indfile";
      $popfile = "K$i.popfile";
      for( my $datatype = 0; $datatype < 2; $datatype++ ){
	print $datatype, "\n";
	if( $datatype == 1 ){
	  $inds = $pops;
	  $outfile = "K$i.POPQ";
	  $miscfile = "K$i.POPQ.miscfile";
	  $paramfile = "paramfile.k$i.pop";
	}else{
	  $inds = $numind;
	  $outfile = "K$i.INDIVQ";
	  $miscfile = "K$i.INDIVQ.miscfile";
	  $paramfile = "paramfile.k$i.ind";
	}
	&clumppwrite( $paramfile, $datatype, $indfile, $popfile, $outfile, $miscfile, $i, $inds, $krun, $method, $weight, $similarity, $greed );
	print CLUMPPCOM "clumpp $paramfile\n";
      }
    }
    close CLUMPPCOM;
    system( "cat $clumppcom | parallel" );# == 0 or die "ERROR: System call to clumpp failed: $?\n\n";
  }else{
    $indfile = "K$clusters.indfile";
    $popfile = "K$clusters.popfile";
    for( my $datatype = 0; $datatype < 2; $datatype++ ){
      #print $datatype, "\n";
      if( $datatype == 1 ){
	$inds = $pops;
	$outfile = "K$clusters.POPQ";
	$miscfile = "K$clusters.POPQ.miscfile";
	$paramfile = "paramfile.k$clusters.pop";
      }else{
	$inds = $numind;
	$outfile = "K$clusters.INDIVQ";
	$miscfile = "K$clusters.INDIVQ.miscfile";
	$paramfile = "paramfile.k$clusters.ind";
      }
      &clumppwrite( $paramfile, $datatype, $indfile, $popfile, $outfile, $miscfile, $clusters, $inds, $krun, $method, $weight, $similarity, $greed );
      system( "clumpp $paramfile" );
    }
  }
  
  
  # change directories back to the original working directory
  chdir $wd;
}

# if the user desires the files to be sorted by Q value, similar to the functionality of the Structure GUI, then this can be done here using the -Q flag
if( $sortq == 1 ){
	print "clumpall = $clumppall\n";
	chdir $harvestdir;
	if( $clumppall == 1 ){
		for( my $i = $startk; $i < $endk+1; $i++ ){
			my $indfile = "K$i.INDIVQ";
			#$indfile = $harvestdir . "/" . $indfile;
			my $popfile = "K$i.POPQ";
			#$popfile = $harvestdir . "/" . $popfile;
			&qsort( $indfile, $popfile, $i );
		}
	}else{
		my $indfile = "K$clusters.INDIVQ";
		print "clusters is $clusters\n";
		print "indfile is $indfile\n";
		#$indfile = $harvestdir . "/" . $indfile;
		my $popfile = "K$clusters.POPQ";
		#$popfile = $harvestdir . "/" . $popfile;
		&qsort( $indfile, $popfile, $clusters );
	}
	chdir $wd;
}

# run distruct if -D flag is used
if( $distruct == 1){
  chdir $harvestdir;
  
  # write labels file for distruct plot
  &labels( $structuredatafile, $harvestdir );

  # see if distruct is installed
  my( $program ) = "distruct";
  my( $distructpres ) = &progtest( $program );
  if( $distructpres == 0 ){
    die "ERROR:\tInstall ", $program, " before proceeding\n\tOr run without -C option\n\n";
  }
  my( $drawparams, $indivq, $popq, $outfile );
  
  if( $clumppall == 1 ){
    my $distructcom = "distructcommands.txt";
    open( DISTRUCTCOM, '>', $distructcom );
    for( my $i = $startk; $i < $endk+1; $i++ ){
      $drawparams = "drawparams.k$i";
      $indivq = "K$i.INDIVQ";
      $popq = "K$i.POPQ";
      $outfile = "K$i.ps";
      if( $sortq == 1 ){
		$pops = $i;
	  }
      &distructwrite( $drawparams, $indivq, $popq, $outfile, $i, $numind, $pops, $bottomlabels, $blab, $toplabels, $tlab );
      
      print DISTRUCTCOM "distruct -d $drawparams\n";
    }
    close DISTRUCTCOM;
    
    system( "cat $distructcom | parallel" ); # == 0 or die "ERROR: System call to distruct failed: $?\n\n";
  }else{
    $drawparams = "drawparams.k$clusters";
    $indivq = "K$clusters.INDIVQ";
    $popq = "K$clusters.POPQ";
    $outfile = "K$clusters.ps";
	if( $sortq == 1 ){
		$pops = $clusters;
	}
    &distructwrite( $drawparams, $indivq, $popq, $outfile, $clusters, $numind, $pops, $bottomlabels, $blab, $toplabels, $tlab );
    system( "distruct -d $drawparams" );
  }
  chdir $wd;
}

exit;

#####################################################################################################
############################################ Subroutines ############################################
#####################################################################################################

# subroutine to print help
sub help{
  
  print "\nParallel Structure is a perl script developed by Steven Michael Mussmann\n\n";
  print "To report bugs send an email to mussmann\@email.uark.edu\n";
  print "When submitting bugs please include all input files, options used for the program, and all error messages that were printed to the screen\n\n";
  
  print "parallel_structure.pl\n";
  print "\tThe following flags require no additional input.  They are switches to turn options on/off.\n";
  print "\tAny of these options can be combined.  For example:\n";
  print "\t\tparallel_structure.pl -apPtuw\n\n";
  
  print "\tCommands that are general options for the parallel_structure.pl script\n";
  print "\t\t[ -C | -D | -E | -h | -p | -S | -w ]\n";
  print "\tFlag\tExplanation\n";
  print "\t-C:\tThis flag will run the program CLUMPP if used.\n";
  print "\t-D:\tThis flag will run the program DISTRUCT if used.\n";
  print "\t-E:\tThis flag will run the Evanno method if the structureharvester python script is present.\n";
  print "\t-h:\tThis flag displays this help message.  If invoked, all other options will be ignored and only the help message will be displayed.\n";
  print "\t-p:\tThis flag will invoke parallelization of structure using the program parallel.\n";
  print "\t\tThis flag also sets [RANDOMIZE] in the EXTRAPARAMS file to 0.\n";
  print "\t-S:\tThis flag will run the program Structure using the files generated by this script.\n";
  print "\t-w:\tThis flag will overwrite all previous output from the parallel_structure.pl script without prompting the user for permission.\n\n";
  
  print "\tCommands that are entered into the MAINPARAMS file:\n";
  print "\t\t[ -d | -e | -F | -L | -M | -n | -P | -R | -t | -u | -v | -y ]\n";
  print "\tFlag\tStructure Equivalent\tExplanation\n";
  print "\t-d:\t[PHASEINFO]\t\tUse this flag if phase information is provided in the data file.\n";
  print "\t-e:\t[PHASED]\t\tUse this flag if the data is in the correct phase.\n";
  print "\t-F:\t[POPFLAG]\t\tUse this flag to turn off the POPFLAG option in Structure if the popflag is used in your input file.\n";
  print "\t-L:\t[LOCDATA]\t\tUse this flag to turn on the LOCDATA option in Structure if populations are identified in your input file.\n";
  print "\t-M:\t[MAPDISTANCES]\t\tUse this flag if your input file contains information on may locations for individual markers.\n";
  print "\t-P:\t[POPDATA]\t\tUse this flag when populations are not identified in the input file.\n";
  print "\t-n:\t[MARKERNAMES]\t\tUse this flag when your input file does not contain a row of marker names.\n";
  print "\t-R:\t[RECESSIVEALLELES]\tUse this flag when you are using dominant markers such as AFLP.\n";
  print "\t-t:\t[ONEROWPERIND]\t\tUse this flag when the input file has more than one row per individual.\n";
  print "\t-u:\t[LABEL]\t\t\tUse this flag when the input file does not contain unique identifiers for each individual.\n";
  print "\t-v:\t[MARKOVPHASE]\t\tUse this flag when the phase follows a markov chain.\n";
  print "\t-y:\t[PHENOTYPE]\t\tUse this flag if your input file contains phenotypic information.\n\n";
  
  print "\tCommands that are entered into the EXTRAPARAMS file:\n";
  print "\t\t[ -a | -H | -I | -N | -o | -O | -s | -U | -z ]\n";
  print "\tFlag\tStructure Equivalent\tExplanation\n";
  print "\t-a:\t[NOADMIX]\t\tUsing this flag turns on the admixture model\n";
  print "\t-H:\t[INFERALPHA]\t\tUse this flag if you do not want to model the parameter alpha from the data.  Command is ignored in the noadmix model.\n";
  print "\t-I:\t[LINKAGE]\t\tUsing this flag turns on the linkage model\n";
  print "\t-N:\t[INFERLAMBDA]\t\tUse this flag if you want Structure to infer a suitable value for lambda\n";
  print "\t-o:\t[LOCPRIOR]\t\tUsing this flag turns on the locprior model\n";
  print "\t-O:\t[POPSPECIFICLAMBDA]\tUse this flag if you want to infer a separate lambda for each population\n";
  print "\t-s:\t[ONEFST]\t\tUse this flag if you want to assume the same FST value for all populations\n";
  print "\t-U:\t[USEPOPINFO]\t\tIf this flag is used the -P flag must also be used\n";
  print "\t-z:\t[POPALPHAS]\t\tUse this flag if you want to infer a separate alpha for each population.\n\n";
  
  print "The following flags are required.  Each of the following requires an additional argument.\n";
  print "\tFlag\tExplanation\n";
  print "\t-f\tSpecify the input file [required]\n\n";
  print "\t-i\tSpecify the number of individuals in the input file [required]\n\n";
  print "\t-l\tSpecify the number of loci in the input file [required]\n\n";
  print "\t-X\tSpecify the number of sampled populations in the input file [required]\n\n";
  
  print "The following flags are not required.  However, if used each will require an additional argument. Default values are shown below.\n";
  print "\tFlag\tExplanation\n";
  print "\t-A\tSpecify the alpha value for the degree of admixture\n";
  print "\t\tDefault = 1.0\n\n";
  print "\t-b\tSpecify the burnin value for Structure\n";
  print "\t\tDefault = 100,000 generations\n\n";
  print "\t-B\tSpecify the lambda value for Structure\n";
  print "\t\tDefault = 1.0\n\n";
  print "\t-g\tSpecify the number of generations to be run by Structure\n";
  print "\t\tDefault = 1,000,000 generations\n\n";
  print "\t-k\tSpecify the starting K value\n";
  print "\t\tDefault = 1\n\n";
  print "\t-K\tSpecify ending maximum K value to be tested by structure\n";
  print "\t\tDefault = 10\n\n";
  print "\t-m\tSpecify the value that represents missing data\n";
  print "\t\tDefault = -9\n\n";
  print "\t-r\tSpecify the number of runs for each K value\n";
  print "\t\tDefault = 10 runs\n\n";
  print "\t-x\tSpecify the number of extra columns in the input file\n";
  print "\t\tDefault = 0\n\n";
  print "\t-Y\tSpecify the ploidy value\n";
  print "\t\tDefault = 2\n\n";
  
  # clumpp options
  print "The following commands are given to the CLUMPP program.\n";
  print "\t\t[ -c | -G | -q ]\n";
  print "\tFlag\tCLUMPP Option\tExplanation\n";
  print "\t-c:\tM\t\tSets the search method to be used.  1 = FullSearch, 2 = Greedy, 3 = LargeKGreedy (Default = 3)\n";
  print "\t-G:\tGREEDY_OPTION\tSets the input order if the Greedy or LargeKGreedy options are used.  1 = All possible input orders, 2 = random input orders, 3 = pre-specified input orders.  Default = 2\n";
  print "\t-q:\tna\t\tRuns CLUMPP for all K values\n\n";
  
  # distruct options
  print "The following commands are given to the DISTRUCT program.\n";
  print "\t\t[ -Q | -W | -Z ]\n";
  print "\tFlag\tDISTRUCT Option\t\tExplanation\n";
  print "\tna\tSorts individuals by Q value for printing in Distruct\n";
  print "\t-W:\tINFILE_LABEL_ATOP\tTakes the name of the top labels row for Distruct.\n\n";
  print "\t-Z:\tINFILE_LABEL_BELOW\tTakes the name of the bottom labels row for Distruct.\n\n";
}

#####################################################################################################

# subroutine tests if a program is both installed in user's path and executable when the name of a program is passed to it.
sub progtest{
  
  my( $prog ) = @_;
  
  my( $pres );
  
  if( grep{ -x "$_/$prog" }split /:/, $ENV{PATH} ){
    print "\n", $prog, " is installed\n\n";
    $pres = 1;
  }else{
    $pres = 0;
  }
  
  return( $pres );
}

#####################################################################################################

sub getdaytime{
  
  my( $wd ) = @_;
  
  #get information to name directory for old data based on computer time
  my( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime(time);
  $year+=1900;
  my @mabbr = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
  
  # declare directories for moving objects
  my $oldrun = "Run" . $mday . $mabbr[$mon] . $year . "_" . $hour . $min . $sec;
  my $oldrundir = $wd . "/" . $oldrun;
  
  return( $oldrun, $oldrundir );
}

#####################################################################################################
# subroutine tests if files are present and prompts user if he/she wishes to overwrite those files
# program will be terminated if user wishes to not overwrite any of the files

# If any of these files exists you may choose to kill the program without overwriting any of the files
sub promptoverwrite{
  
  # take array of files that have been passed to the subroutine
  my( $files, $dirs, $wd, $oldrundir  ) = @_;
  
  # check to see if each file in the array already exists
  foreach my $item( @$files ){
    my $filename = "$wd/$item";
    my( $answer ) = &yesno( $filename );
    if ($answer == 1 ) {
      next;
    }elsif ($answer == 0 ){
      my @save;
      push( @save, $item );
      &copyoldfiles( \@save, $oldrundir );
    }
  }
  
  foreach my $item( @$dirs ){
    my( $answer ) = &yesno( $item );
    if ($answer == 1 ) {
      next;
    }elsif ($answer == 0 ){
      my @save;
      push( @save, $item );
      &copyoldfolder( \@save, $oldrundir );
    }
  }
}

#####################################################################################################

# subroutine to prompt for yes/no input

sub yesno{
  
  my( $item ) = @_;
  
  my $answer;
  
  if (-e $item) {
    # if file exists, prompt if it should be overwritten
    print "\n$item already exists.\n";
    my $key = '';
    # Program will continue to prompt you for an answer unless your answer begins with an N or a Y
    while ($key !~ /^Y|^N/i ){
      print "\nDo you want to overwrite $item? (y/n)\n";
      $key = <STDIN>;
    }
    # If yes, move on to the next file
    if ($key =~ /^Y/i ){
      $answer = 1;
      next;
      # If no, kill the program
    }elsif ($key =~ /^N/i ){
      $answer = 0;
    }
    # if the file does not exist, move to the next file in the array
  }else{
    next;
  }
  
  return $answer;
}



#####################################################################################################

# subroutine to copy old Structure run into directory named by current year, mont, day, and hour, minute, second when perl script was executed.
sub copyoldfolder{
  
  my( $dirs, $oldrundir ) = @_;
  
  foreach my $source( @$dirs ){
    my @temp = split( /\//, $source );
    my $last = pop(@temp);
    my $target = $oldrundir . "/" . $last;
    print $target, "\n";
    dirmove( $source, $target ) or die "Couldn't copy old $last directory to $target: $!";
  }
  
}

#####################################################################################################

# subroutine to copy old scripts used to run Structure

sub copyoldfiles{
  
  my( $files, $targetdir ) = @_;
  
  foreach my $file( @$files ){
    my @temp = split( /\//, $file );
    my $last = pop( @temp );
    my $target = $targetdir . "/" . $file;
    print $target, "\n";
    move( $file, $target );
  }
  
}

#####################################################################################################

# subroutine to write the mainparams file

sub main{
  
  my( $burn, $data, $endk, $excols, $flag, $gen, $loc, $mainout, $map, $markov, $miss, $mname, $numind, $numloci, $phased, $phaseinfo, $pheno, $ploidy, $pop, $rec, $row, $uid ) = @_;
  
  open ( MAINPARAMS, '>', $mainout ) or die "Can't open $mainout: $!\n";
  
  print MAINPARAMS "Basic program parameters\n";
  print MAINPARAMS "#define MAXPOPS\t\t", $endk, "\n";
  print MAINPARAMS "#define BURNIN\t\t", $burn, "\n";
  print MAINPARAMS "#define NUMREPS\t\t", $gen, "\n\n";
  print MAINPARAMS "Input file\n";
  print MAINPARAMS "#define INFILE\t", $data, "\n\n";
  print MAINPARAMS "Data format and other parameters\n";
  print MAINPARAMS "#define EXTRACOLS\t", $excols, "\n";
  print MAINPARAMS "#define LABEL\t", $uid, "\n";
  print MAINPARAMS "#define LOCDATA\t", $loc, "\n";
  print MAINPARAMS "#define MAPDISTANCES\t", $map, "\n";
  print MAINPARAMS "#define MARKERNAMES\t", $mname, "\n";
  print MAINPARAMS "#define MARKOVPHASE\t", $markov, "\n";
  print MAINPARAMS "#define MISSING\t", $miss, "\n";
  print MAINPARAMS "#define NUMINDS\t", $numind, "\n";
  print MAINPARAMS "#define NUMLOCI\t", $numloci, "\n";
  print MAINPARAMS "#define ONEROWPERIND\t", $row, "\n";
  print MAINPARAMS "#define PHASED\t", $phased, "\n";
  print MAINPARAMS "#define PHASEINFO\t", $phaseinfo, "\n";
  print MAINPARAMS "#define PHENOTYPE\t", $pheno, "\n";
  print MAINPARAMS "#define PLOIDY\t", $ploidy, "\n";
  print MAINPARAMS "#define POPDATA\t", $pop, "\n";
  print MAINPARAMS "#define POPFLAG\t", $flag, "\n";
  print MAINPARAMS "#define RECESSIVEALLELES\t", $rec, "\n";
  
  close MAINPARAMS;
}

#####################################################################################################

# subroutine to write extraparams file

sub extra{
  
  my( $admix, $alpha, $extraout, $fst, $infalpha, $inflam, $lambda, $link, $locprior, $popalphas, $popspeclam, $randomize, $usepop ) = @_;
  
  open ( EXTRAPARAMS, '>', $extraout ) or die "Can't open $extraout: $!\n";
  
  print EXTRAPARAMS "#define ALPHA\t", $alpha, "\n";
  print EXTRAPARAMS "#define INFERALPHA\t", $infalpha, "\n";
  print EXTRAPARAMS "#define INFERLAMBDA\t", $inflam, "\n";
  print EXTRAPARAMS "#define LAMBDA\t", $lambda, "\n";
  print EXTRAPARAMS "#define LINKAGE\t", $link, "\n";
  print EXTRAPARAMS "#define LOCPRIOR\t", $locprior, "\n";
  print EXTRAPARAMS "#define NOADMIX\t", $admix, "\n";
  print EXTRAPARAMS "#define ONEFST\t", $fst, "\n";
  print EXTRAPARAMS "#define POPALPHAS\t", $popalphas, "\n";
  print EXTRAPARAMS "#define POPSPECIFICLAMBDA\t", $popspeclam, "\n";
  print EXTRAPARAMS "#define RANDOMIZE\t", $randomize, "\n";
  print EXTRAPARAMS "#define USEPOPINFO\t", $usepop, "\n";
  
  close EXTRAPARAMS;
}

#####################################################################################################

sub parstruct{
  
  my( $data, $com, $endk, $startk, $strun ) = @_;
  
  open( STRUN, '>', $strun ) or die "Can't open $strun: $!\n";
  
  print STRUN "#! /usr/bin/perl\n";
  print STRUN "use warnings;\n";
  print STRUN "use strict;\n";
  print STRUN "use Cwd;\n";
  print STRUN "use File::Path;\n";
  # get current working directory
  print STRUN "my \$dir = cwd();\n";
  # specify names of paths
  print STRUN "my \$results = \$dir . \"/results\";\n";
  print STRUN "my \$log = \$dir . \"/log\";\n";
  print STRUN "my \$harvester = \$dir . \"/harvester\";\n";
  # remove directories if they exist
  print STRUN "rmtree( [\$results, \$log, \$harvester], 1, 1 );\n";
  # make directories
  print STRUN "mkpath( [\$results, \$log, \$harvester], 1, 0755 ) or die \"Couldn't make directories in \$dir: \$!\";\n";;
  print STRUN "chdir \$log;\n";
  print STRUN "mkpath( [";
  for( my $i = $startk; $i < $endk; $i++ ){
    print STRUN "'k$i', ";
  }
  print STRUN "'k$endk'], 1, 0755 ) or die \"Couldn't make directories in \$log: \$!\";\n";
  print STRUN "chdir \$dir;\n";
  print STRUN "system( \"cat structurecom.txt | parallel\" );\n";
  print STRUN "my \$zipfile = \$dir . \"/harvester/$data.results.zip\";";
  print STRUN "my \$zipfolder = \$dir . \"/results/\";";
  print STRUN "system( \"zip -rj \$zipfile \$zipfolder\" );\n";
  print STRUN "exit;\n";
  
  close STRUN;
}

#####################################################################################################

sub struct{
  
  my( $data, $endk, $krun, $startk, $strun ) = @_;
  
  open( STRUN, '>', $strun ) or die "Can't open $strun: $!\n";
  
  print STRUN "#! /usr/bin/perl\n";
  print STRUN "use warnings;\n";
  print STRUN "use strict;\n";
  print STRUN "use Cwd;\n";
  print STRUN "use File::Path;\n";
  # get current working directory
  print STRUN "my \$dir = cwd();\n";
  # specify names of paths
  print STRUN "my \$results = \$dir . \"/results\";\n";
  print STRUN "my \$log = \$dir . \"/log\";\n";
  print STRUN "my \$harvester = \$dir . \"/harvester\";\n";
  # set starting k value, ending k value, and number of runs at each k
  print STRUN "my \$startk = $startk;\n";
  print STRUN "my \$endk = $endk;\n";
  print STRUN "my \$krun = $krun;\n";
  # remove directories if they exist
  print STRUN "rmtree( [\$results, \$log, \$harvester], 1, 1 );\n";
  # make directories
  print STRUN "mkpath( [\$results, \$log, \$harvester], 1, 0755 ) or die \"Couldn't make directories in \$dir: \$!\";\n";;
  print STRUN "chdir \$log;\n";
  print STRUN "mkpath( [";
  for( my $i = $startk; $i < $endk; $i++ ){
    print STRUN "'k$i', ";
  }
  print STRUN "'k$endk'], 1, 0755 ) or die \"Couldn't make directories in \$log: \$!\";\n";
  print STRUN "chdir \$dir;\n";
  print STRUN "for( my \$i = \$startk; \$i < \$endk+1; \$i++ ){\n";
  print STRUN "\tfor (my \$j = 1; \$j < \$krun+1; \$j++ ){\n";
  print STRUN "\t\tsystem( \"structure -K \$i -m mainparams -o results/$data.k\$i.run\$j 2>&1 | tee log/k\$i/$data.k\$i.run\$j.log\" );\n";
  print STRUN "\t}\n";
  print STRUN "}\n";
  print STRUN "system( \"zip -rj $data.results.zip results/\" );\n";
  print STRUN "exit;\n";    
  
  close STRUN;
}

#####################################################################################################

sub strcommand{
  
  my( $com, $data, $endk, $krun, $startk  ) = @_;
  
  # Generate array of random seeds
  # This is necessary because parallelizing structure in the manner used in this script can result in several independent runs starting simultaneously.
  # Structure generates random numbers based on the computer clock (if RANDOMIZE=1) so this means multiple independent runs could be identical if they start at the same time.
  # Perl is being used to generate an array of random number seeds that will be entered on the command line to remedy this problem.
  my @seeds;
  # Initialize the seed counter which will be increased in the loop that prints the structure commands
  my $seedcount = 0;
  # Determine number of seeds to generate.
  my $nruns = ((($endk - $startk)+1)*$krun);
  # create the array of random seeds
  for (my $i = 0; $i < $nruns; $i++){
    $seeds[$i] = int(rand(2000000000));
  }
  
  # Write the structurecommands.txt file
  open ( STCOM, '>', $com ) or die "Can't open $com: $!\n";
  for (my $i = $startk; $i < $endk+1; $i++ ){
    for (my $j = 1; $j < $krun+1; $j++ ){
      print STCOM "structure -K $i -m mainparams -o results/", $data,"_k", $i, "_run$j -D $seeds[$seedcount] 2>&1 | tee log/k$i/", "$data", "_k", $i, "_run$j.log\n";
      $seedcount++;
    }
  }
  
  close STCOM;
}

#####################################################################################################

sub evannoparse{
  
  # take name of evanno.txt file from input
  my( $efile ) = @_;
  
  # declare arrays to hold k-values and delta-k values
  my @kvalue;
  my @deltak;
  
  # open the evanno.txt file and get lines containing deltak values
  open( EFILE, $efile ) or die "Can't open $efile: $!\n";
  while( my $line = <EFILE> ){
    if( $line =~ /^\d/ ){
      chomp( $line );
      my @temp = split( /\t/, $line );
      push( @kvalue, $temp[0] );
      push( @deltak, $temp[6] );
    }
  }
  close EFILE;
  # remove the first and last elements of each array because they will be NA values for deltak.
  pop @kvalue;
  pop @deltak;
  shift @kvalue;
  shift @deltak;
  
  # get maximum deltak value
  my $index = &maxarray( \@deltak );
  my $clusters = $kvalue[$index];
  
  return $clusters;
}

#####################################################################################################

# This subroutine finds the maximum value in an array and returns its index

sub maxarray{
  my ($array) = @_;
  my $index = 0;
  if (not @$array) {
    die ("Empty array\n");
  }
  for (my $i = 0, my  $max = 0; $i < scalar(@$array); $i++) {
    if (@$array[$i] > $max) {
      $max = @$array[$i];
      $index = $i;
    }
  }
  
  return $index;
}
#####################################################################################################
# This subroutine calculates the mean of an array

sub mean{
    my($datain) = @_;
    if (not @$datain) {
	die("Empty array\n");
    }
    my $total = 0;
    for ( my $i = 0; $i < scalar (@$datain); $i++) {
	$total += @$datain[$i];
	# print "my new total is ", $total, "\n", 
    }
    my $sublength = scalar (@$datain);
    my $average = $total / $sublength;
    return $average;
}

#####################################################################################################

sub clumppwrite{
  
  # line to take input data
  my( $paramfile, $datatype, $indfile, $popfile, $outfile, $miscfile, $i, $inds, $krun, $method, $weight, $similarity, $greed ) = @_;
  
  # open file
  open( PARAMFILE, '>', $paramfile );
  
  print PARAMFILE "DATATYPE ", $datatype, "\n";
  print PARAMFILE "INDFILE ", $indfile, "\n";
  print PARAMFILE "POPFILE ", $popfile, "\n";
  print PARAMFILE "OUTFILE ", $outfile, "\n";
  print PARAMFILE "MISCFILE ", $miscfile, "\n";
  print PARAMFILE "K ", $i, "\n";
  print PARAMFILE "C ", $inds, "\n";
  print PARAMFILE "R ", $krun, "\n";
  print PARAMFILE "M ", $method, "\n";
  print PARAMFILE "W ", $weight, "\n";
  print PARAMFILE "S ", $similarity, "\n";
  print PARAMFILE "GREEDY_OPTION ", $greed, "\n";
  print PARAMFILE "REPEATS 100000\n";
  print PARAMFILE "PERMUTATIONFILE permutationfile\n";
  print PARAMFILE "PRINT_PERMUTED_DATA 1\n";
  print PARAMFILE "PERMUTED_DATAFILE perm_datafile\n";
  print PARAMFILE "PRINT_EVERY_PERM 0\n";
  print PARAMFILE "EVERY_PERMFILE every_permfile\n";
  print PARAMFILE "PRINT_RANDOM_INPUTORDER 0\n";
  print PARAMFILE "RANDOM_INPUTORDERFILE random_inputorderfile\n";
  print PARAMFILE "OVERRIDE_WARNINGS 0\n";
  print PARAMFILE "ORDER_BY_RUN 1\n";
  
  close PARAMFILE;
}

#####################################################################################################

sub distructwrite{
  
  my( $drawparams, $indivq, $popq, $outfile, $i, $numind, $pops, $bottomlabels, $blab, $toplabels, $tlab ) = @_;
  
  # open file
  open( DRAWPARAMS, '>', $drawparams );
  
  print DRAWPARAMS "#define INFILE_POPQ ", $popq, "\n";
  print DRAWPARAMS "#define INFILE_INDIVQ ", $indivq, "\n";
  print DRAWPARAMS "#define INFILE_LABEL_BELOW ", $bottomlabels, "\n";
  print DRAWPARAMS "#define INFILE_LABEL_ATOP ", $toplabels, "\n";
  print DRAWPARAMS "#define INFILE_CLUST_PERM /home/mussmann/local/src/distruct1.1/ColorBrewer/BrBG_", $i, "_div\n";
  print DRAWPARAMS "#define OUTFILE ", $outfile, "\n";
  print DRAWPARAMS "#define K ", $i, "\n";
  print DRAWPARAMS "#define NUMPOPS ", $pops, "\n";
  print DRAWPARAMS "#define NUMINDS ", $numind, "\n";
  print DRAWPARAMS "#define PRINT_INDIVS 1\n";
  print DRAWPARAMS "#define PRINT_LABEL_ATOP 1\n";
  print DRAWPARAMS "#define PRINT_LABEL_BELOW 0\n";
  print DRAWPARAMS "#define PRINT_SEP 1\n";
  print DRAWPARAMS "#define FONTHEIGHT 6\n";
  print DRAWPARAMS "#define DIST_ABOVE -110\n";
  print DRAWPARAMS "#define DIST_BELOW -50\n";
  print DRAWPARAMS "#define BOXHEIGHT 100\n";
  print DRAWPARAMS "#define INDIVWIDTH 2\n";
  print DRAWPARAMS "#define ORIENTATION 1\n";
  print DRAWPARAMS "#define XORIGIN 200\n";
  print DRAWPARAMS "#define YORIGIN 10\n";
  print DRAWPARAMS "#define XSCALE 1\n";
  print DRAWPARAMS "#define YSCALE 1\n";
  print DRAWPARAMS "#define ANGLE_LABEL_ATOP 270\n";
  print DRAWPARAMS "#define ANGLE_LABEL_BELOW 270\n";
  print DRAWPARAMS "#define LINEWIDTH_RIM 3\n";
  print DRAWPARAMS "#define LINEWIDTH_SEP 1\n";
  print DRAWPARAMS "#define LINEWIDTH_IND 1\n";
  print DRAWPARAMS "#define GRAYSCALE 0\n";
  print DRAWPARAMS "#define ECHO_DATA 1\n";
  print DRAWPARAMS "#define REPRINT_DATA 1\n";
  print DRAWPARAMS "#define PRINT_INFILE_NAME 0\n";
  print DRAWPARAMS "#define PRINT_COLOR_BREWER 1\n";
  # write commands to file
  
  # close file
  close DRAWPARAMS;
}
#####################################################################################################

sub qsort{
  
  my( $indfile, $popfile, $clusters ) = @_;
  
  # Declaring my hash of hashes of arrays
  # First hash key = population; second hash key = individual identifier; array = proportion of ancestry attributed to each population for each individual
  my %indhohoa; # hohoa to hold data from indivq file
  my %pophohoa; # hohoa to hold data from popq file
  my $ind = 1; # Individual counter to specify the individual identifier key (2nd hash) in the hohoa
  
  open ( INPUT, $indfile ) || die "Can't open $indfile: $!\n";
  while (my $line = <INPUT> ){
    # Regex to capture the string of probabilities
    if( $line =~ /.+:(.+)/ ){
      chomp ($line);
      # Assign the capture ($1) to a temporary variable name
      my $temp = $1;
      # Regex to remove leading and trailing whitespace
      $temp =~ s/^\s+|\s+$//g;
      # Split the string into the array @probs
      my @probs = split (/\s+/, $temp);
      # Use the subroutine maxarray to return the index of the greatest value in @probs.  The (array index + 1) is the population with the highest assignment probability ($pop)
      my $pop = ( &maxarray(\@probs) + 1);
      # add individual to the array
      unshift( @probs, $ind );
      # put arrays into %indhohoa to sort ind file
      $indhohoa{$pop}{$ind} = \@probs;
      # also put values into %pophohoa to recalculate values in pop file
      for (my $i = 1; $i < @probs; $i++){
	push ( @{$pophohoa{$pop}{$i}}, $probs[$i] ); 
      }
      # Increase the individual counter by 1
      $ind++;
    }
  }
  
  # Close the input file
  close INPUT;
  
  # Open the output indfile
  open (INDFILE, '>', $indfile) || die "Can't open $indfile: $!\n";
  
  my $counter = 1;
  # print populations in order, simple numeric sort of the population numbers
  foreach my $population ( sort {$a <=> $b } keys %indhohoa ){
    foreach my $individual ( sort { $indhohoa{$population}{$b}[$population] <=> $indhohoa{$population}{$a}[$population] } keys %{$indhohoa{$population}} ){
      print INDFILE $counter, "\t", $individual, "\t(0)\t", $population, "\t:\t";
      for (my $i = 1; $i < @{$indhohoa{$population}{$individual}}; $i++){
	print INDFILE $indhohoa{$population}{$individual}[$i], "\t";
      }
      print INDFILE "\n";
      $counter++;
    }
  }
  
  close INDFILE;
  
  print "Individual Q values have been written to $indfile \n";
  
  open (POPFILE, '>', $popfile) || die "Can't open $popfile: $!\n\n";
  foreach my $popassign (sort {$a <=> $b } keys %pophohoa ){
    print POPFILE $popassign, ":\t";
    my $length = 0;
    foreach my $popprob ( sort {$a <=> $b} keys %{$pophohoa{$popassign}} ){
      my $popavg = &mean( \@{$pophohoa{$popassign}{$popprob}} );
      my $rounded = sprintf( "%.4f", $popavg );
      print POPFILE $rounded, "\t";
      $length = scalar( @{$pophohoa{$popassign}{$popprob}} );
    }
    print POPFILE $length, "\n";
  }
  
  # Close the output indfile
  
  close POPFILE;
  
  print "Population Q values have been written to $popfile \n\n";
  
}

#####################################################################################################
# subroutine to produce a labels file for structure output

sub labels{

	my( $strfilepath, $harvesterpath ) = @_;
	if( -f $strfilepath ){
		my $labelfile = join( '/', $harvesterpath, "toplabels" );

		open(STR, $strfilepath) or die "Can't open $strfilepath: $!\n\n";
		my %hash;

		while( my $line = <STR> ){
			chomp( $line );
			my @temp = split( /\s/, $line );
			$hash{$temp[1]} = $temp[0];
		}

		close STR;
	
		open( OUT, '>', $labelfile ) or die "Can't open $labelfile: $!\n\n";
		foreach my $pop( sort {$a<=>$b} keys %hash ){
			my $name;
			if( $hash{$pop} =~ // ){
				$name = $1;
			}
			print OUT $pop, "\t", $name, "\n";
		}
  		close OUT;
	}else{
		print "Couldn't open $strfilepath: does not exist.\n\n";
	}


}

#####################################################################################################

# subroutine to parse the command line options

sub parsecom{ 
  
  my $params = shift( @_ );
  my %opts = %$params;
  
  # set default values for command line arguments
  
  if( $opts{w} and $opts{j} ){
    die "Program cannot be run using both options -w (overwrite) and -j (copy old runs).  Pick one or the other.\n\n";
  }
  
  my $data = $opts{f} or die "\nNo input file specified\n\n";
  my $startk = $opts{k} || "1";
  my $endk = $opts{K} || "10";
  my $burn = $opts{b} || "100000";
  my $gen = $opts{g} || "1000000";
  my $krun = $opts{r} || "10";
  my $numind = $opts{i} || die "\nSpecify number of individuals in input file\n\n";
  my $numloci = $opts{l} || die "\nSpecify number of loci in input file\n\n";
  my $miss = $opts{m} || "-9";
  my $excols = $opts{x} || "0";
  my $ploidy = $opts{Y} || "2";
  my $alpha = $opts{A} || "1.0";
  my $lambda = $opts{B} || "1.0";
  my $pops = $opts{X} || die "\nSpecify number of populations in input file\n\n";
  
  if( $endk < $startk){
    die "\nEnding K value (option -K) must be larger than starting K value (option -k).\n\n";
  }
  
  my $admix;
  if( $opts{a} ){
    $admix = 0;
  }else{
    $admix = 1;
  }
  
  my $row;
  if( $opts{t} ){
    $row = 0;
  }else{
    $row = 1;
  }
  
  my $uid;
  if( $opts{u} ){
    $uid = 0;
  }else{
    $uid = 1;
  }
  
  my $pop;
  if( $opts{P} ){
    $pop = 0;
  }else{
    $pop = 1;
  }
  
  my $flag;
  if( $opts{F} ){
    $flag = 0;
  }else{
    $flag = 1;
  }
  
  my $loc;
  if( $opts{L} ){
    $loc = 1;
  }else{
    $loc = 0;
  }
  
  my $pheno;
  if( $opts{Y} ){
    $pheno = 1;
  }else{
    $pheno = 0;
  }
  
  my $rec;
  if( $opts{R} ){
    $rec = 1;
  }else{
    $rec = 0;
  }
  
  my $map;
  if( $opts{M} ){
    $map = 1;
  }else{
    $map = 0;
  }
  
  my $phaseinfo;
  if( $opts{d} ){
    $phaseinfo = 1;
  }else{
    $phaseinfo = 0;
  }
  
  my $phased;
  if( $opts{e} ){
    $phased = 1;
  }else{
    $phased = 0;
  }
  
  my $mname;
  if( $opts{n} ){
    $mname = 0;
  }else{
    $mname = 1;
  }
  
  my $markov;
  if( $opts{v} ){
    $markov = 1;
  }else{
    $markov = 0;
  }
  
  my $link;
  if( $opts{I} ){
    $link = 1;
  }else{
    $link = 0;
  }
  
  my $usepop;
  if( $opts{U} ){
    $usepop = 1;
    if( $pop == 0 ){
      die "Must use -P flag if -U flag is used\n\n";
    }
  }else{
    $usepop = 0;
  }
  
  my $locprior;
  if( $opts{o} ){
    $locprior = 1;
  }else{
    $locprior = 0;
  }
  
  my $fst;
  if( $opts{s} ){
    $fst = 1;
  }else{
    $fst = 0;
  }
  
  my $infalpha;
  if( $opts{H} ){
    $infalpha = 0;
  }else{
    $infalpha = 1;
  }
  
  my $popalphas;
  if( $opts{z} ){
    $popalphas = 1;
  }else{
    $popalphas = 0;
  }
  
  my $inflam;
  if( $opts{N} ){
    $inflam = 1;
  }else{
    $inflam = 0;
  }
  
  my $popspeclam;
  if( $opts{O} ){
    $popspeclam = 1;
  }else{
    $popspeclam = 0;
  }
  
  # option to run Structure
  my $structure;
  if( $opts{S} ){
    $structure = 1;
  }else{
    $structure = 0;
  }
  
  # option to run the Evanno method
  my $evanno;
  if( $opts{E} ){
    $evanno = 1;
  }else{
    $evanno = 0;
  }
  
  # Structure will run using automatically generated random number seeds.  This is only changed by this script if Structure is run in parallel.
  my $randomize = 1;
  
  # clumpp options
  # M 
  my $method = $opts{c} || "3";
  # GREEDY_OPTION
  my $greed = $opts{G} || "2";
  
  # option to run CLUMPP
  my $clumpp;
  if( $opts{C} ){
    $clumpp = 1;
  }else{
    $clumpp = 0;
  }
  
  my $clumppall;
  if( $opts{q} ){
    $clumppall = 1;
  }else{
    $clumppall = 0;
  }
  
  # distruct options
  #
  
  # option to run DISTRUCT
  my $distruct;
  if( $opts{D} ){
    $distruct = 1;
  }else{
    $distruct = 0;
  }
  
  # file name for bottom labels for distruct
  my $bottomlabels = $opts{Z} || "bottomlabels";
  my $blab;
  if( $opts{Z} ){
    $blab = 1;
  }else{
    $blab = 0;
  }
  
  # file name for top labels for distruct
  my $toplabels = $opts{W} || "toplabels";
  my $tlab;
  if( $opts{W} ){
    $tlab = 1;
  }else{
    $tlab = 0;
  }
  
  my $sortq;
  if( $opts{Q} ){
	$sortq = 1;
  }else{
	$sortq = 0;
  }
  
  return( $data, $startk, $endk, $burn, $gen, $krun, $numind, $numloci, $miss, $excols, $ploidy, $alpha, $lambda, $admix, $row, $uid, $pop, $flag, $loc, $pheno, $rec, $map, $phaseinfo, $phased, $mname, $markov, $link, $usepop, $locprior, $fst, $infalpha, $popalphas, $inflam, $popspeclam, $randomize, $evanno, $structure, $pops, $method, $clumpp, $greed, $distruct, $bottomlabels, $blab, $toplabels, $tlab, $clumppall, $sortq );
}

#####################################################################################################
