// Just a basic macro to create simple plots from given hists.
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2019
// Run it with following command:
// root -l -b -q GeneralPlotterMultiples.C

void GeneralPlotterMultiples(){

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

  gROOT->ForceStyle();

  TString path_pre = "/nfs/dust/cms/user/flabe/TstarTstar/data/";
  std::vector<TString> files = {"/DNN/UL17/hadded/uhh2.AnalysisModuleRunner.DATA.DATA.root", "/DNN/UL17/hadded/uhh2.AnalysisModuleRunner.MC.TTbar.root", "/DNN/UL17/hadded/uhh2.AnalysisModuleRunner.MC.ST.root", "/DNN/UL17/hadded/uhh2.AnalysisModuleRunner.DATA.datadrivenBG.root", "datadriven_direct/uhh2.AnalysisModuleRunner.MC.Other.root"};
  std::vector<TString> steps = {"newTaggerCR", "newTaggerCR", "newTaggerCR", "newTaggerCR", "reco"};
  std::vector<TString> hists = {"pt_ST_rebinned", "pt_ST_rebinned", "pt_ST_rebinned", "pt_ST_rebinned", "pt_ST_rebinned"};
  std::vector<TString> labels = {"DATA", "TTbar", "ST", "datadriven", "datadriven_real"};

  bool doNormalize = true;

  Double_t w = 800;
  Double_t h = 600;

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  c1_hist->SetLogy();
  TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.2,1.0,1.0);
  TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.2);
  pad1->Draw();
  pad1->SetLogy();
  pad2->Draw();
  auto legend = new TLegend(0.1+0.6,0.8,0.3+0.6,0.9);


  pad1->cd();
  for(int i = 0; i < files.size(); i++){
    std::cout << "Plotting: " << files.at(i) << " // " << steps.at(i) << " // " << hists.at(i) << endl;
    TFile *input = TFile::Open(path_pre+files.at(i));
    TH1D *hist = (TH1D*)input->Get(steps.at(i)+"/"+hists.at(i)); //histogram
    if(!hist) cout<<"Hist is empty"<<endl;
    if(doNormalize) hist->Scale(1./hist->Integral());
    hist->SetMarkerStyle(0);
    hist->SetLineColor(i+1);
    hist->SetMarkerColor(i+1);
    hist->SetLineStyle(1);
    if(i == 0) hist->Draw("hist");
    else hist->Draw("hist same");
    legend->AddEntry(hist, labels.at(i));
  }

  legend->Draw("same");

  // ratios to first
  pad2->cd();
  TFile *input_base = TFile::Open(path_pre+files.at(0));
  TH1D *hist_base = (TH1D*)input_base->Get(steps.at(0)+"/"+hists.at(0)); //histogram
  if(doNormalize) hist_base->Scale(1./hist_base->Integral());
  for(int i = 1; i < files.size(); i++){
    std::cout << "Ratio for: " << files.at(i) << " // " << steps.at(i) << " // " << hists.at(i) << endl;
    TFile *input = TFile::Open(path_pre+files.at(i));
    TH1D *hist = (TH1D*)input->Get(steps.at(i)+"/"+hists.at(i)); //histogram
    if(!hist) cout<<"Hist is empty"<<endl;
    if(doNormalize) hist->Scale(1./hist->Integral());
    hist->Divide(hist_base);
    hist->SetMarkerStyle(0);
    hist->SetLineColor(i+1);
    hist->SetLineStyle(1);
    hist->GetYaxis()->SetRangeUser(0.1, 3);
    if(i == 0) hist->Draw("hist");
    else hist->Draw("hist same");
    //legend->AddEntry(hist, labels.at(i));
  }

  c1_hist->SaveAs("plots/plot.pdf");
}
