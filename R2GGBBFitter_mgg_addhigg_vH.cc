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
void MakePlots(RooWorkspace* w, Float_t Mass);

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
  style();
  TString fileBaseName(TString::Format("hgg.hig.mH%.1f_8TeV.vh", mass));
  TString card_name("models_test.rs"); // put the model parameters here!
  HLFactory hlf("HLFactory", card_name, false);

  RooFitResult* fitresults;
  bool cutbased=true;
  // the minitree to be addeed
  //
  TString ssignal = "/afs/cern.ch/work/a/acarvalh/CMSSW_6_1_1/src/code/Limit_codes/MiniTrees/v28/wzh_m125_8TeV_m260.root";
  //
  RooWorkspace* w = hlf.GetWs();
  AddSigData(w, mass,ssignal);
  cout<<"SIGNAL ADDED "<<ssignal<<endl;
  SigModelFit(w, mass); // constructing signal pdf
  MakeSigWS(w, fileBaseName);
  //MakePlots(w, mass);
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
    Float_t nExpEvt = sigToFit[c]->sumEntries();
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
  wAll->factory("prod::CMS_hgg_sig_m0_vh_cat0(mgg_sig_m0_cat0, CMS_hgg_sig_m0_absShift)");
  wAll->factory("prod::CMS_hgg_sig_m0_vh_cat1(mgg_sig_m0_cat1, CMS_hgg_sig_m0_absShift)");
  // (3) Systematics on resolution
  wAll->factory("CMS_hgg_sig_sigmaScale[1,1,1]");
  wAll->factory("prod::CMS_hgg_sig_sigma_vh_cat0(mgg_sig_sigma_cat0, CMS_hgg_sig_sigmaScale)");
  //wAll->factory("prod::CMS_hgg_sig_sigma_cat0(mgg_sig_sigma_cat0, CMS_hgg_sig_sigmaScale)");

  wAll->factory("prod::CMS_hgg_sig_sigma_vh_cat1(mgg_sig_sigma_cat1, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_gsigma_vh_cat0(mgg_sig_gsigma_cat0, CMS_hgg_sig_sigmaScale)");
  wAll->factory("prod::CMS_hgg_sig_gsigma_vh_cat1(mgg_sig_gsigma_cat1, CMS_hgg_sig_sigmaScale)");
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
TString::Format("EDIT::CMS_hgg_sig_vh_cat%d(mggSig_cat%d,",c,c) +
TString::Format(" mgg_sig_m0_cat%d=CMS_hgg_sig_m0_vh_cat%d, ", c,c) +
TString::Format(" mgg_sig_sigma_cat%d=CMS_hgg_sig_sigma_vh_cat%d, ", c,c) +
// TString::Format(" mgg_sig_alpha_cat%d=CMS_hgg_sig_alpha_cat%d, ", c,c) +
// TString::Format(" mgg_sig_n_cat%d=CMS_hgg_sig_n_cat%d, ", c,c) +
// TString::Format(" mgg_sig_frac_cat%d=CMS_hgg_sig_frac_cat%d, ", c,c) +
TString::Format(" mgg_sig_gsigma_cat%d=CMS_hgg_sig_gsigma_vh_cat%d)", c,c)
  );
  TString filename(wsDir+TString(fileBaseName)+".inputsig.root");
  wAll->writeToFile(filename);
  cout << "Write signal workspace in: " << filename << " file" << endl;
  return;
} // close make signal WP
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void MakePlots(RooWorkspace* w, Float_t Mass) {
  const Int_t ncat = NCAT;
  std::vector<TString> catdesc;
  catdesc.push_back(" 2 btag");
  catdesc.push_back(" 1 btag");
  catdesc.push_back("cat 2");
  catdesc.push_back("cat 3");
  // retrieve data sets from the workspace
  // RooDataSet* dataAll = (RooDataSet*) w->data("Data");
  //RooDataSet* signalAll = (RooDataSet*) w->data("Sig");
  //RooDataSet* higgsAll = (RooDataSet*) w->data("Hig");
  // blinded dataset
  // RooDataSet* data[ncat];
  RooDataSet* sigToFit[ncat];
  RooAbsPdf* mggGaussSig[ncat];
  RooAbsPdf* mggCBSig[ncat];
  RooAbsPdf* mggSig[ncat];
  //
  RooAbsPdf* mggBkg[ncat];
  for (int c = 0; c < ncat; ++c) {
  // data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("Sig_cat%d",c));
    mggGaussSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggGaussSig_cat%d",c));
    mggCBSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggCBSig_cat%d",c));
    mggSig[c] = (RooAbsPdf*) w->pdf(TString::Format("mggSig_cat%d",c));
    mggBkg[c] = (RooAbsPdf*) w->pdf(TString::Format("mggBkg_cat%d",c));
  } // close categories
  RooRealVar* mgg = w->var("mgg");
  mgg->setUnit("GeV");
  RooAbsPdf* mggGaussSigAll = w->pdf("mggGaussSig");
  RooAbsPdf* mggCBSigAll = w->pdf("mggCBSig");
  RooAbsPdf* mggSigAll = w->pdf("mggSig");
  //RooAbsPdf* mggBkgAll = w->pdf("mggBkg_cat1");
  //
  //****************************//
  // Plot mgg Fit results
  //****************************//
  // Set P.D.F. parameter names
  // WARNING: Do not use it if Workspaces are created
  // SetParamNames(w);
  Float_t minSigFit(120),maxSigFit(130);
  Float_t MASS(Mass);
  Int_t nBinsMass(20); // just need to plot
  //RooPlot* plotmggAll = mgg->frame(Range(minSigFit,maxSigFit),Bins(nBinsMass));
  //signalAll->plotOn(plotmggAll);
  gStyle->SetOptTitle(0);
  TCanvas* c1 = new TCanvas("c1","mgg",0,0,500,500);
  c1->cd(1);
  //********************************************//
  // Plot Signal Categories
  //****************************//
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  RooPlot* plotmgg[ncat];
  for (int c = 0; c < ncat; ++c) {
    plotmgg[c] = mgg->frame(Range(minSigFit,maxSigFit),Bins(nBinsMass));
    sigToFit[c]->plotOn(plotmgg[c],LineColor(kWhite),MarkerColor(kWhite));
    mggSig[c] ->plotOn(plotmgg[c]);
    double chi2n = plotmgg[c]->chiSquare(0) ;
    cout << "------------------------- Experimentakl chi2 = " << chi2n << endl;
    mggSig[c] ->plotOn(
plotmgg[c],
Components(TString::Format("GaussSig_cat%d",c)),
LineStyle(kDashed),LineColor(kGreen));
    mggSig[c] ->plotOn(
plotmgg[c],
Components(TString::Format("CBSig_cat%d",c)),
LineStyle(kDashed),LineColor(kRed));
    mggSig[c] ->paramOn(plotmgg[c]);
    sigToFit[c] ->plotOn(plotmgg[c]);
// TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
    TH1F *hist = new TH1F("hist", "hist", 400, minSigFit, maxSigFit);
    plotmgg[c]->SetTitle("CMS preliminary 19.7/fb ");
    plotmgg[c]->SetMinimum(0.0);
    plotmgg[c]->SetMaximum(1.40*plotmgg[c]->GetMaximum());
    plotmgg[c]->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV)");
    TCanvas* ctmp = new TCanvas("ctmp","Background Categories",0,0,500,500);
    plotmgg[c]->Draw();
    plotmgg[c]->Draw("SAME");
    TLegend *legmc = new TLegend(0.62,0.75,0.99,0.99);
    legmc->AddEntry(plotmgg[c]->getObject(5),"Simulation","LPE");
    legmc->AddEntry(plotmgg[c]->getObject(1),"Parametric Model","L");
    legmc->AddEntry(plotmgg[c]->getObject(2),"Gaussian Outliers","L");
    legmc->AddEntry(plotmgg[c]->getObject(3),"Crystal Ball component","L");
    legmc->SetHeader(" ");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();
    // float effS = effSigma(hist);
    TLatex *lat = new TLatex(
minSigFit+0.5,0.85*plotmgg[c]->GetMaximum(),
" Resonance - 300 GeV");
    lat->Draw();
    TLatex *lat2 = new TLatex(
minSigFit+1.5,0.75*plotmgg[c]->GetMaximum(),catdesc.at(c));
    lat2->Draw();
    ///////
    char myChi2buffer[50];
    sprintf(myChi2buffer,"#chi^{2}/ndof = %f",chi2n);
    TLatex* latex = new TLatex(0.52, 0.7, myChi2buffer);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    //latex -> Draw("same");
    ctmp->SaveAs(TString::Format("sigmodel_vh_cat%d.pdf",c));
    ctmp->SaveAs(TString::Format("sigmodel_vh_cat%d.png",c));
    //ctmp->SaveAs(TString::Format("sigmodel_cat%d.C",c));
  } // close categories
    return;
} // close makeplots signal
////////////////////////////////////////////////////////////////////////
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
void style(){
  TStyle *defaultStyle = new TStyle("defaultStyle","Default Style");
  defaultStyle->SetOptStat(0000);
  defaultStyle->SetOptFit(000);
  defaultStyle->SetPalette(1);
  /////// pad ////////////
  defaultStyle->SetPadBorderMode(1);
  defaultStyle->SetPadBorderSize(1);
  defaultStyle->SetPadColor(0);
  defaultStyle->SetPadTopMargin(0.05);
  defaultStyle->SetPadBottomMargin(0.13);
  defaultStyle->SetPadLeftMargin(0.13);
  defaultStyle->SetPadRightMargin(0.02);
  /////// canvas /////////
  defaultStyle->SetCanvasBorderMode(0);
  defaultStyle->SetCanvasColor(0);
  defaultStyle->SetCanvasDefH(600);
  defaultStyle->SetCanvasDefW(600);
  /////// frame //////////
  defaultStyle->SetFrameBorderMode(0);
  defaultStyle->SetFrameBorderSize(1);
  defaultStyle->SetFrameFillColor(0);
  defaultStyle->SetFrameLineColor(1);
  /////// label //////////
  defaultStyle->SetLabelOffset(0.005,"XY");
  defaultStyle->SetLabelSize(0.05,"XY");
  defaultStyle->SetLabelFont(42,"XY");
  /////// title //////////
  defaultStyle->SetTitleOffset(1.1,"X");
  defaultStyle->SetTitleSize(0.01,"X");
  defaultStyle->SetTitleOffset(1.25,"Y");
  defaultStyle->SetTitleSize(0.05,"Y");
  defaultStyle->SetTitleFont(42, "XYZ");
  /////// various ////////
  defaultStyle->SetNdivisions(505,"Y");
  defaultStyle->SetLegendBorderSize(0); // For the axis titles:

    defaultStyle->SetTitleColor(1, "XYZ");
    defaultStyle->SetTitleFont(42, "XYZ");
    defaultStyle->SetTitleSize(0.06, "XYZ");
 
    // defaultStyle->SetTitleYSize(Float_t size = 0.02);
    defaultStyle->SetTitleXOffset(0.9);
    defaultStyle->SetTitleYOffset(1.05);
    // defaultStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

    // For the axis labels:
    defaultStyle->SetLabelColor(1, "XYZ");
    defaultStyle->SetLabelFont(42, "XYZ");
    defaultStyle->SetLabelOffset(0.007, "XYZ");
    defaultStyle->SetLabelSize(0.04, "XYZ");

    // For the axis:
    defaultStyle->SetAxisColor(1, "XYZ");
    defaultStyle->SetStripDecimals(kTRUE);
    defaultStyle->SetTickLength(0.03, "XYZ");
    defaultStyle->SetNdivisions(510, "XYZ");
    defaultStyle->SetPadTickX(1); // To get tick marks on the opposite side of the frame
    defaultStyle->SetPadTickY(1);
    defaultStyle->cd();
  return;
}





