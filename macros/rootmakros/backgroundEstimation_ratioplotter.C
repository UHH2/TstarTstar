// Just a basic macro to create ratio plot of two hists.
// author: F. Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q ratioplotter.C

void backgroundEstimation_ratioplotter(){

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
  gStyle->SetPadBottomMargin(0.4);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.1);

  gROOT->ForceStyle();
  Double_t w = 800;
  Double_t h = 600;

  TString MCpath = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/UL17/hadded/";
  TString filename_base = "uhh2.AnalysisModuleRunner.MC.";
  std::vector<TString> samples_to_replace = {"WJets", "QCD", "VV"};
  TString MCfolder = "newTaggerCR";
  TString MChistname = "pt_ST_rebinned";

  TString datadrivenPath = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/UL17/hadded/";
  TString datadrivenFilename = "uhh2.AnalysisModuleRunner.DATA.datadrivenBG.root";
  TString datadrivenFolder = "newTaggerCR";
  TString datadrivenHistname = "pt_ST_rebinned";

  TString label = "S_{T}";

  TCanvas *canvas = new TCanvas("canvas", "c", w, h);

  TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.35,1.0,1.0);
  TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.35);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  pad1->SetLogy();

  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);

  // getting and adding previous stuff
  bool first = true;
  TH1D *hist_MC;
  for (auto sample : samples_to_replace) {
    TFile *input_MC = TFile::Open(MCpath+filename_base+sample+".root");
    if (first) hist_MC = (TH1D*)input_MC->Get(MCfolder+"/"+MChistname);
    else hist_MC->Add((TH1D*)input_MC->Get(MCfolder+"/"+MChistname));
    first = false;
  }
  if(!hist_MC) std::cout << "MC does not exist!" << std::endl;


  // getting datadriven stuff
  TFile *input_datadriven = TFile::Open(datadrivenPath+datadrivenFilename);
  TH1D *hist_datadriven = (TH1D*)input_datadriven->Get(datadrivenFolder+"/"+datadrivenHistname); //histogram
  if(!hist_datadriven) std::cout << "DD does not exist!" << std::endl;

  TH1D* hist_DD_clone = (TH1D*) hist_datadriven->Clone();

  hist_MC->SetTitle("");
  hist_MC->GetXaxis()->SetTickLength(0);
  hist_MC->GetXaxis()->SetLabelOffset(999);
  hist_MC->GetYaxis()->SetTitleOffset(0.75);
  hist_MC->GetYaxis()->SetLabelSize(0.075);
  hist_MC->GetYaxis()->SetTitleSize(0.075);
  hist_MC->GetYaxis()->SetTitle("events");
  hist_MC->SetMarkerStyle(20);
  hist_MC->SetMarkerColor(1);
  hist_MC->SetLineColor(1);
  hist_MC->Draw("hist");
  hist_datadriven->SetMarkerStyle(20);
  hist_datadriven->SetMarkerColor(2);
  hist_datadriven->SetLineColor(2);
  hist_datadriven->Draw("hist same");

  auto legend = new TLegend(0.55,0.7,0.78,0.88);
  gStyle->SetLegendTextSize(0.05);
  legend->AddEntry(hist_MC,"MC-based","l");
  legend->AddEntry(hist_datadriven,"datadriven","l");
  legend->SetBorderSize(0);
  legend->Draw();

  pad2->cd();


  // plot ratio
  hist_DD_clone->Divide(hist_MC);
  hist_DD_clone->SetTitle("");
  hist_DD_clone->GetXaxis()->SetTitle("S_{T} [GeV]");
  hist_DD_clone->GetYaxis()->SetTitle("ratio");
  hist_DD_clone->SetLineColor(2);

  // axis styline
  hist_DD_clone->GetXaxis()->SetLabelSize(0.1);
  hist_DD_clone->GetXaxis()->SetTitleSize(0.15);
  hist_DD_clone->GetYaxis()->SetLabelSize(0.1);
  hist_DD_clone->GetYaxis()->SetNdivisions(505);
  hist_DD_clone->GetYaxis()->SetTitleOffset(0.);
  hist_DD_clone->GetYaxis()->SetTitleSize(10);

  hist_DD_clone->Draw("");
  canvas->SaveAs("plots/BGest_eval.pdf");

}
