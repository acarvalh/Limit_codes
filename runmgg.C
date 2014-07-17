{
  gROOT->LoadMacro("R2GGBBFitter_mgg_addhiggs.cc");
  TString cuts_tmp = "";
  TString baseDir = "";
  runfits(125.6,270,1,false,cuts_tmp,baseDir);
}
