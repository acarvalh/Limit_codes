#!/usr/bin/perl

# ----------------------------------------------------------------------------
#      MAIN PROGRAM
# ----------------------------------------------------------------------------

use Env;

#PG lettura dei parametri da cfg file
#PG --------------------------------
print "reading ".$ARGV[0]."\n" ;

open (USERCONFIG,$ARGV[0]) ;

while (<USERCONFIG>)
{
    chomp; 
    s/#.*//;                # no comments
    s/^\s+//;               # no leading white
    s/\s+$//;               # no trailing white
#    next unless length;     # anything left?
    my ($var, $value) = split(/\s*=\s*/, $_, 2);
    $User_Preferences{$var} = $value;
}

$BASEDir          = $User_Preferences{"BASEDir"};
$EXEName          = $User_Preferences{"EXEName"} ;
$JOBTemplate   = $User_Preferences{"JOBTemplate"} ;
$Mass             = $User_Preferences{"Mass"} ;
$NumOfSteps        = $User_Preferences{"NumOfSteps"} ;
$Var              = $User_Preferences{"Var"} ;
$Min              = $User_Preferences{"Min"} ;
$Max              = $User_Preferences{"Max"} ;
$Step             = $User_Preferences{"Step"} ;
$EXEName          = $User_Preferences{"EXEName"} ;
$QUEUE            = $User_Preferences{"QUEUE"} ;


print "BASEDir = "          .$BASEDir."\n" ;
print "EXEName = "          .$EXEName."\n" ;
print "JOBTemplate = "   .$JOBTemplate."\n" ;
print "Mass = "   .$Mass."\n" ;
print "NumOfSteps = "   .$NumOfSteps."\n" ;
print "Var = "   .$Var."\n" ;
print "Min = "   .$Min."\n" ;
print "Max = "   .$Max."\n" ;
print "Step = "   .$Step."\n" ;
print "EXEName = "   .$EXEName."\n" ;
print "QUEUE = "   .$QUEUE."\n" ;



$sampleJobListFile = "./lancia.sh";
open(SAMPLEJOBLISTFILE, ">", $sampleJobListFile);


#####################################################
# PG prepare the array containing the root files list
#####################################################


   system("cd ".$BASEDir."\n");
    
    chomp($_);
    
    $sample1 = "JOBS";
   
    print("Sample: ".$sample1."\n") ;  

    system ("rm -r ".$sample1." \n") ;
    system ("mkdir ".$sample1." \n") ;
    
  
    ################
    # loop over jobs 
    ################
    
    $jobID = -1;
    
    for($jobIt = 0; $jobIt < $NumOfSteps; ++$jobIt)
    { 
      for($jobIt2 = 0; $jobIt2 < $NumOfSteps; ++$jobIt2)
      { 
        $jobID = $jobID +1;

	$currDir = `pwd` ;
	chomp ($currDir) ;
    
	$jobDir = $currDir."/".$sample1."/JOB_".$jobID ;
	system ("mkdir ".$jobDir." \n") ;
    
	$tempBjob = $jobDir."/bjob_".$jobID.".sh" ;
	$command = "touch ".$tempBjob ;
	system ($command) ;
	$command = "chmod 777 ".$tempBjob ;
	system ($command) ;
        
        $cutMin = $Min - $jobIt*$Step;
        $cutMax = $Max + $jobIt2*$Step;

        $tempo1 = "./tempo1" ;
        $MACRO = $BASEDir."../R2GGBBFitter_mtot_range.cc";
        system ("cat ".$JOBTemplate." | sed -e s%MACRO%".$MACRO.
"%g > ".$tempo1) ;

	$tempo2 = "./tempo2" ;
        system ("cat ".$tempo1." | sed -e s%VAR%".$Var.                                                                 "%g > ".$tempo2) ;

        $tempo3 = "./tempo3" ;
        system ("cat ".$tempo2." | sed -e s%MIN%".$cutMin.                                                                "%g > ".$tempo3) ;

        $tempo4 = "./tempo4" ;
        system ("cat ".$tempo3." | sed -e s%MAX%".$cutMax.                                                                "%g > ".$tempo4) ;

        $tempo5 = "./tempo5" ;
        system ("cat ".$tempo4." | sed -e s%BASEDIR%".$BASEDir.                                                           "%g > ".$tempo5) ;
    
        $tempo6 = "./tempo6" ;
        system ("cat ".$tempo5." | sed -e s%MASS%".$Mass.                                                                 "%g > ".$tempo6) ;
    
	$JOBCfgFile = $jobDir."/".$EXEName ;
	system ("mv ".$tempo6." ".$JOBCfgFile) ;
        
        system ("rm ./tempo*") ;

        
    ######################
    # make job files
    ######################
    
        open (SAMPLEJOBFILE, ">", $tempBjob) or die "Can't open file ".$tempBjob;

        $command = "#!/bin/tcsh" ;
        print SAMPLEJOBFILE $command."\n";

        $command = "cd ".$BASEDir ;
        print SAMPLEJOBFILE $command."\n";

        $command = "setenv SCRAM_ARCH slc6_amd64_gcc462" ;
        print SAMPLEJOBFILE $command."\n";
    
        $command = "eval `scramv1 ru -csh`" ;
        print SAMPLEJOBFILE $command."\n";
    
        $command = "cd -" ;
        print SAMPLEJOBFILE $command."\n";

        $command = "root -l -q -b ".$JOBCfgFile ;
        print SAMPLEJOBFILE $command."\n";

        $command = "cp hgg.inputbkg_8TeV.root ".$jobDir."/";
        print SAMPLEJOBFILE $command."\n";

        $command = "cp hgg.mH".$Mass.".0_8TeV.inputsig.root ".$jobDir."/";
        print SAMPLEJOBFILE $command."\n";

        $command = "cp hgg.mH".$Mass.".0_8TeVrep.txt ".$jobDir."/";
        print SAMPLEJOBFILE $command."\n";

        $command = "combine -M Asymptotic ".$jobDir."/hgg.mH".$Mass.".0_8TeVrep.txt >> limit_".$Mass."_".$cutMin."_".$cutMax.".txt";
        print SAMPLEJOBFILE $command."\n";

        $command = "cp limit_".$Mass."_".$cutMin."_".$cutMax.".txt ".$BASEDir."/../";
        print SAMPLEJOBFILE $command."\n";

    ############
    # submit job
    ############

        $command = "bsub -cwd ".$jobDir." -q ".$QUEUE." ".$tempBjob."\n" ;
        print SAMPLEJOBLISTFILE $command."\n";
    
      }
    }
    
print "NumOfSteps = "   .$NumOfSteps."\n" ;
