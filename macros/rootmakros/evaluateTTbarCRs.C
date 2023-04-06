// Decorrelate DNN tagger, see https://confluence.desy.de/x/0up4DQ
// author: F.Labe
// date: 24.09.2021
// Run it with following command:
// root -l -b -q decorrelatedTagger.C

void evaluateTTbarCRs(){

  Double_t cut_value = 0;

  // some style options
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleYOffset(1.20);
  gStyle->SetTitleXOffset(1.0);
  gStyle->SetPalette(55);

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.17);

  gROOT->ForceStyle();
  Double_t w = 800;
  Double_t h = 600;

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  c1_hist->SetLogz();

  // get histogram
  TString filename="TTbar";
  std::vector<TString> otherFiles = {"ST", "WJets", "QCD", "VV"};
  std::vector<TString> regions = {"ttbarCR_v1", "ttbarCR_v2", "ttbarCR_v3"};
  TString histname="pt_ST_rebinned";
  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/UL18/hadded/";
  TString fileprefix = "uhh2.AnalysisModuleRunner.MC.";

  for (auto region : regions) {

    // TTbar
    TFile *input = TFile::Open(path+fileprefix+filename+".root");
    if(!input) cout << "Empty file" << endl;
    TH1D *hist = (TH1D*)input->Get(region+"/"+histname);
    if(!hist) cout << "Empty hist" << endl;
    double ttbar_int = hist->Integral();

    bool first = true;
    TH1D *otherhist;
    for (auto otherFile : otherFiles) {
      TFile *otherinput = TFile::Open(path+fileprefix+otherFile+".root");
      TH1D *temphist = (TH1D*)otherinput->Get(region+"/"+histname);
      if(first) otherhist = temphist;
      else otherhist->Add(temphist);
    }

    double other_int = otherhist->Integral();

    std::cout << "ttbar purity in " << region << ": " << 100*ttbar_int / (ttbar_int+ other_int) << "%" << std::endl;

  }



}
