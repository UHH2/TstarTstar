// Efficiency histograms for trigger related studies
// author: F.Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q cutflowPlots.C

void kinCutPtDep(TString suffix = ""){

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

  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.17);

  gROOT->ForceStyle();
  Double_t w = 800;
  Double_t h = 600;

  // Reuseable buffers
  TFile *input;
  TH1D *hist;

  // Drawing Definitions
  TCanvas *canvas = new TCanvas("chist", "c", w, h);  

  // Defining paths
  TString path_pre = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Selection/";
  TString path_post = "/hadded/";
  TString histname = "N_jets";

  int year = 2018;

  // Defining Steps
  // AfterBtag
  TString step_before = "AfterBtag";
  TString step_after = "AfterdR";

  // Defining Samples
  TString fileprefix = "uhh2.AnalysisModuleRunner.MC.";
  TString sampleBasename = "TstarTstar_M-";

  // Defining Drawing options
  std::vector<int> colors = {1, 1, 1, 1};

  // ########################
  // ## Finish Definitions ##
  // ########################

  // Saving initial values
  std::vector<double> signal_before;
  for(uint i = 7; i <= 20; i++){
    if((year == 2017 || year == 2018) && i == 11) {
      signal_before.push_back(1);
      continue;
    }
    input = TFile::Open(path_pre+std::to_string(year)+path_post+fileprefix+sampleBasename+std::to_string(i*100)+".root");
    hist = (TH1D*)input->Get(step_before+"/"+histname);
    signal_before.push_back(hist->Integral());
  }

  // Saving after values
  std::vector<double> signal_after;
  for(uint i = 7; i <= 20; i++){
    if((year == 2017 || year == 2018) && i == 11) {
      signal_after.push_back(1);
      continue;
    }
    input = TFile::Open(path_pre+std::to_string(year)+path_post+fileprefix+sampleBasename+std::to_string(i*100)+".root");
    hist = (TH1D*)input->Get(step_after+"/"+histname);
    signal_after.push_back(hist->Integral());
  }

  TH1* hist_total = new TH1D("hist", "hist", 14, 650, 2050);

  for (uint i = 0; i < 14; i++){
    double val = signal_after.at(i)/signal_before.at(i);
    hist_total->SetBinContent(i+1, val);
  }

  // ##########################
  // ## Finish Filling Hists ##
  // ##########################

  // removing missing bins
  if(year == 2017 || year == 2018) hist_total->SetBinContent(5, 0);

  std::cout << "drawing..." << endl;

  hist_total->SetLineWidth(3);
  hist_total->SetMaximum(1.1);
  hist_total->SetMinimum(0.6);
  hist_total->SetTitle(step_before+" to "+step_after+" cut efficiency");
  hist_total->GetYaxis()->SetTitle("Efficiency");
  hist_total->GetXaxis()->SetTitle("Tstar mass");
  hist_total->Draw();

  canvas->SaveAs("eff_"+std::to_string(year)+"_"+step_before+"_to_"+step_after+".pdf");

}
