using namespace RooFit;
using namespace RooStats ;

const Int_t NCAT = 2;

// declare the functions
void AddSigData(RooWorkspace*, Float_t);
void AddBkgData(RooWorkspace*);
void SigModelFit(RooWorkspace*, Float_t);
void MakePlots(RooWorkspace*, Float_t, RooFitResult* );
void MakeSigWS(RooWorkspace* w, const char* filename);
void MakeBkgWS(RooWorkspace* w, const char* filename);
void MakeDataCardonecat(RooWorkspace* w, const char* filename, const char* filename1);
void MakeDataCardREP(RooWorkspace* w, const char* filename, const char* filename1);
void MakeDataCardLnU(RooWorkspace* w, const char* filename, const char* filename1);
void SetParamNames(RooWorkspace*);
void SetConstantParams(const RooArgSet* params);

RooFitResult* fitresult[NCAT]; // container for the fit results
using namespace RooFit;
using namespace RooStats ;
const Int_t NCAT = 2;
bool addHiggs=true;
void AddSigData(RooWorkspace*, Float_t);
void SigModelFit(RooWorkspace*, Float_t);
void MakeSigWS(RooWorkspace* w, const char* filename);
void SetParamNames(RooWorkspace*);
void SetConstantParams(const RooArgSet* params);

RooFitResult* fitresult[NCAT]; // container for the fit results
RooFitResult* BkgModelFitBernstein(RooWorkspace*, Bool_t);

RooArgSet* defineVariables()
{
  RooRealVar* mgg = new RooRealVar("mgg","M(#gamma#gamma)",100,180,"GeV");
  //RooRealVar* mtot = new RooRealVar("mtot","M(#gamma#gammajj)",200,1600,"GeV");
  //RooRealVar* mjj = new RooRealVar("mjj","M(jj)",100,1600,"GeV");
  RooRealVar* evWeight = new RooRealVar("evWeight","HqT x PUwei",0,100,"");
  RooCategory* cut_based_ct = new RooCategory("cut_based_ct","event category 4") ;
  //
  cut_based_ct->defineType("cat4_0",0);
  cut_based_ct->defineType("cat4_1",1);
  //
  RooArgSet* ntplVars = new RooArgSet(*mgg, * cut_based_ct, *evWeight);
  ntplVars->add(*mgg);
  //ntplVars->add(*mtot);
  //ntplVars->add(*mjj);
  ntplVars->add(*cut_based_ct);
  return ntplVars;
}

void runfits(const Float_t mass=120, Int_t mode=1, Bool_t dobands = false)
{
  TString fileBaseName(TString::Format("hgg.hig.mH%.1f_8TeV.ggh", mass));
  TString card_name("models_test.rs"); // put the model parameters here!
  HLFactory hlf("HLFactory", card_name, false);

  RooFitResult* fitresults;
  bool cutbased=true;
  // the minitree to be addeed
  //
  TString ssignal = "/afs/cern.ch/work/o/obondu/public/forRadion/limitTrees/v28/v28_fitToMgg_noKinFit/ggh_m125_powheg_8TeV_m260.root";
  //
  RooWorkspace* w = hlf.GetWs();
  AddSigData(w, mass,ssignal);
  cout<<"SIGNAL ADDED"<<endl;
  SigModelFit(w, mass); // constructing signal pdf
  MakeSigWS(w, fileBaseName);
  cout<<" did ggH WS's"<<endl;
  cout<< "here"<<endl;
  return;
} // close runfits
////////////////////////////////////////////////////////////////////
// we add the data to the workspace in categories
void AddSigData(RooWorkspace* w, Float_t mass, TString signalfile) {
  cout << "================= Add Signal==============================" << endl;
  const Int_t ncat = NCAT;
  Float_t MASS(mass);
  // Luminosity:
  Float_t Lum = 19785.0; // pb-1
  RooRealVar lumi("lumi","lumi",Lum);
  w->import(lumi);
  RooArgSet* ntplVars = defineVariables();
  TFile sigFile(signalfile);
  TTree* sigTree = (TTree*) sigFile.Get("TCVARS");
  // common preselection cut
  TString mainCut("1");
  RooDataSet sigScaled(
"sigScaled",
"dataset",
sigTree,
*ntplVars,
mainCut,
"evWeight");
  cout << "======================================================================" << endl;
  RooDataSet* sigToFit[ncat];
  TString cut0 = " && 1>0";// "&& mtot > 955 && mtot < 1150 "; //
  TString cut1 = " && 1>0";//"&& mtot > 955 && mtot < 1150 "; // "&& 1>0";//
  //
  TString cutj0 = " && 1>0";//"&& mjj_wokinfit > 90 && mjj_wokinfit < 170 "; //"&& 1>0";//
  TString cutj1 = " && 1>0";//"&& mjj_wokinfit > 100 && mjj_wokinfit < 160 "; // "&& 1>0";//
  //
  // we take only mtot to fit to the workspace, we include the cuts
  sigToFit[0] = (RooDataSet*) sigScaled.reduce(
*w->var("mgg"),
mainCut+TString::Format(" && cut_based_ct==%d ",0)+cut0+cutj0);
  w->import(*sigToFit[0],Rename(TString::Format("Sig_cat%d",0)));
    //
  sigToFit[1] = (RooDataSet*) sigScaled.reduce(
*w->var("mgg"),
mainCut+TString::Format(" && cut_based_ct==%d ",1)+cut1+cutj1);
  w->import(*sigToFit[1],Rename(TString::Format("Sig_cat%d",1)));
  // Create full signal data set without categorization
  RooDataSet* sigToFitAll = (RooDataSet*) sigScaled.reduce(*w->var("mgg"),mainCut);
  cout << "======================================================================" << endl;
  w->import(*sigToFitAll,Rename("Sig"));
  // here we print the number of entries on the different categories
  cout << "========= the number of entries on the different categories ==========" << endl;
  cout << "---- one channel: " << sigScaled.sumEntries() << endl;
  for (int c = 0; c < ncat; ++c) {
    Float_t nExpEvt = sigToFitAll[c]->sumEntries();
    cout << TString::Format("nEvt exp. cat%d : ",c) << nExpEvt
<< TString::Format(" eff x Acc cat%d : ",c)
<< "%"
<< endl;
  } // close ncat
  cout << "======================================================================" << endl;
  sigScaled.Print("v");
  return;
} // end add signal function
/////////////////////////////////////////////////////////////////////
// we make the fit model
void SigModelFit(RooWorkspace* w, Float_t mass) {
  const Int_t ncat = NCAT;
  Float_t MASS(mass);
  //******************************************//
  // Fit signal with model pdfs
  //******************************************//
  // four categories to fit
  RooDataSet* sigToFit[ncat];
  RooAbsPdf* mggSig[ncat];
  // fit range
/*
const int minsigfit =mass - 120, maxsigfit=mass +120;
RooDataSet* sigToFit[ncat];
RooAbsPdf* mtotSig[ncat];
// fit range
Float_t minMassFit2(minfit2),minMassFit1(minsigfit),maxMassFit(maxsigfit);
*/
  Float_t minSigFit(115),maxSigFit(135);
  for (int c = 0; c < ncat; ++c) {
    // import sig and data from workspace
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    mggSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggSig_cat%d",c));
      ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->setVal(MASS);
    //RooRealVar* peak = w->var(TString::Format("mgg_sig_m0_cat%d",c));
    //peak->setVal(MASS);
    cout << "OK up to now..." <<MASS<< endl;
    // Fit model as M(x|y) to D(x,y)

    mggSig[c]->fitTo(*sigToFit[c],Range(minSigFit,maxSigFit),SumW2Error(kTRUE));
    cout << "old = " << ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal() << endl;

    double mPeak = ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal()+0.6; // shift the peak
    ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->setVal(mPeak); // shift the peak

    cout << "mPeak = " << mPeak << endl;
    cout << "new mPeak position = " << ((RooRealVar*) w->var(TString::Format("mgg_sig_m0_cat%d",c)))->getVal() << endl;

    // IMPORTANT: fix all pdf parameters to constant, why?
    w->defineSet(TString::Format("SigPdfParam_cat%d",c),
        RooArgSet(
*w->var(TString::Format("mgg_sig_m0_cat%d",c)),
*w->var(TString::Format("mgg_sig_sigma_cat%d",c)),
*w->var(TString::Format("mgg_sig_alpha_cat%d",c)),
*w->var(TString::Format("mgg_sig_n_cat%d",c)),
*w->var(TString::Format("mgg_sig_gsigma_cat%d",c)),
*w->var(TString::Format("mgg_sig_frac_cat%d",c))) );
    SetConstantParams(w->set(TString::Format("SigPdfParam_cat%d",c)));
  } // close for ncat
} // close signal model fit
///////////////////////////////////////////////////////////////
void MakeSigWS(RooWorkspace* w, const char* fileBaseName) {
  TString wsDir = "workspaces/";
  const Int_t ncat = NCAT;
  //**********************************************************************//
  // Write pdfs and datasets into the workspace before to save
  // for statistical tests.
  //**********************************************************************//
  RooAbsPdf* mggSigPdf[ncat];
  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  for (int c = 0; c < ncat; ++c) {
    mggSigPdf[c] = (RooAbsPdf*) w->pdf(TString::Format("mggSig_cat%d",c));
    wAll->import(*w->pdf(TString::Format("mggSig_cat%d",c)));
  }
  // (2) Systematics on energy scale and resolution
  // 1,1,1 statistical to be treated on the datacard
  wAll->factory("CMS_hgg_sig_m0_absShift[1,1,1]");
  wAll->factory("prod::CMS_hgg_sig_m0_cat0(mgg_sig_m0_cat0, CMS_hgg_sig_m0_absShift)");
  wAll->factory("prod::CMS_hgg_sig_m0_cat1(mgg_sig_m0_cat1, CMS_hgg_sig_m0_absShift)");
  // (3) Systematics on resolution
  wAll->factory("CMS_hgg_sig_sigmaScale[1,1,1]");
  wAll->factory("prod::CMS_hgg_sig_sigma_cat0(mgg_sig_sigma_cat0, CMS_hgg_sig_sigmaScale)");
  //wAll->factory("prod::CMS_hgg_sig_sigma_cat0(mgg_sig_sigma_cat0, CMS_hgg_sig_sigmaScale)");

  wAll->factory("prod::CMS_hgg_sig_sigma_cat1(mgg_sig_sigma_cat1, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_gsigma_cat0(mgg_sig_gsigma_cat0, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_gsigma_cat1(mgg_sig_gsigma_cat1, CMS_hgg_sig_sigmaScale)");
  // save the other parameters
/* for (int c = 0; c < ncat; ++c) {
wAll->factory(
TString::Format("CMS_hgg_sig_alpha_cat%d[%g,0.5,5]",
c, wAll->var(TString::Format("mgg_sig_alpha_cat%d",c))->getVal()));
wAll->factory(
TString::Format("CMS_hgg_sig_n_cat%d[%g,0.5,20]",
c, wAll->var(TString::Format("mgg_sig_n_cat%d",c))->getVal()));
wAll->factory(
TString::Format("CMS_hgg_sig_frac_cat%d[%g,0.0,1.0]",
c, wAll->var(TString::Format("mgg_sig_frac_cat%d",c))->getVal()));
}
*/
  // (4) do reparametrization of signal
  for (int c = 0; c < ncat; ++c) wAll->factory(
TString::Format("EDIT::CMS_hgg_sig_cat%d(mggSig_cat%d,",c,c) +
TString::Format(" mgg_sig_m0_cat%d=CMS_hgg_sig_m0_cat%d, ", c,c) +
TString::Format(" mgg_sig_sigma_cat%d=CMS_hgg_sig_sigma_cat%d, ", c,c) +
// TString::Format(" mgg_sig_alpha_cat%d=CMS_hgg_sig_alpha_cat%d, ", c,c) +
// TString::Format(" mgg_sig_n_cat%d=CMS_hgg_sig_n_cat%d, ", c,c) +
// TString::Format(" mgg_sig_frac_cat%d=CMS_hgg_sig_frac_cat%d, ", c,c) +
TString::Format(" mgg_sig_gsigma_cat%d=CMS_hgg_sig_gsigma_cat%d)", c,c)
  );
  TString filename(wsDir+TString(fileBaseName)+".inputsig.root");
  wAll->writeToFile(filename);
  cout << "Write signal workspace in: " << filename << " file" << endl;
  return;
} // close make signal WP
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
void SetConstantParams(const RooArgSet* params) {
  // set constant parameters for signal fit, ... NO IDEA !!!!
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }
} // close set const parameters
////////////////////////////////////////////////////////////////////////





