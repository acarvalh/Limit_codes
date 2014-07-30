{
  gROOT->LoadMacro("MACRO");
  TString cuts_tmp0 = " && VAR > MIN && VAR < MAX";
  TString cuts_tmp1 = " && VAR > MIN && VAR < MAX";
  TString baseDir = "BASEDIR../";
  runfits(MASS,1,false,cuts_tmp0,cuts_tmp1,baseDir);
}
