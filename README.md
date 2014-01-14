Limit_codes
===========

On CMSSW_6_1_1 you install the higgs combination tools

https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit

......................................................................
cd CMSSW_6_1_1/myfolder/Limit_codes
cmsenv
mkdir datacards/
mkdir workspaces/
mkdir Minitrees/ (put limit trees inside)
(edit the path for limit trees on R2GBB... )
. ./fits_4body.sh (to 4 body fit*)
. ./fits.sh (to mgg fit*)

follow the limit trees convetion mXX.root
where XX is the mass in GeV -- eg:400
......................................................................




