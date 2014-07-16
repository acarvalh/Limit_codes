{
  gROOT->LoadMacro("MACRO");
  TString cuts_tmp = " && VAR > MIN && VAR < MAX";
  TString baseDir = "BASEDIR../";
  runfits(MASS,1,false,cuts_tmp,baseDir);
}
