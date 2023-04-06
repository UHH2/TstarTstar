// Just a basic macro to create ratio plot of two hists.
// author: F. Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q ratioplotter.C

void ULcompare_ratioplotter(){

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

  TString sample = "DATA.DATA";

  TString beforePath = "/nfs/dust/cms/user/flabe/TstarTstar/data/Preselection/2018/hadded/";
  TString afterPath = "/nfs/dust/cms/user/flabe/TstarTstar/data/Preselection/UL18/hadded/";
  TString filename = "uhh2.AnalysisModuleRunner."+sample+".root";
  TString histname = "pt_ele";
  TString label = "p^{e}_{T} [GeV]";
  //label = histname;

  TString beforeFolder = "AfterLep";
  TString afterFolder = "AfterLep";

  TCanvas *canvas = new TCanvas("canvas", "c", w, h);

  TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.35,1.0,1.0);
  TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.35);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  pad1->SetLogy();

  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);

  TFile *input_EOY = TFile::Open(beforePath+filename);
  TFile *input_UL = TFile::Open(afterPath+filename);

  TH1D *hist_num = (TH1D*)input_UL->Get(afterFolder+"/"+histname); //histogram
  TH1D *hist_denom = (TH1D*)input_EOY->Get(beforeFolder+"/"+histname); //histogram
  if(!hist_num) std::cout << "Numerator does not exist!" << std::endl;
  if(!hist_denom) std::cout << "Denominator does not exist!" << std::endl;

  TH1D* hist_num_clone = (TH1D*) hist_num->Clone();

  hist_num->SetTitle("");
  hist_num->GetXaxis()->SetTickLength(0);
  hist_num->GetXaxis()->SetLabelOffset(999);
  hist_num->GetYaxis()->SetTitleOffset(0.75);
  hist_num->GetYaxis()->SetLabelSize(0.075);
  hist_num->GetYaxis()->SetTitleSize(0.075);
  hist_num->GetYaxis()->SetTitle("events");
  hist_num->SetMarkerStyle(20);
  hist_num->SetMarkerColor(1);
  hist_num->SetLineColor(1);
  hist_num->Draw("hist");
  hist_denom->SetMarkerStyle(20);
  hist_denom->SetMarkerColor(2);
  hist_denom->SetLineColor(2);
  hist_denom->Draw("hist same");

  auto legend = new TLegend(0.55,0.7,0.78,0.88);
  gStyle->SetLegendTextSize(0.05);
  legend->AddEntry(hist_num,"UL","l");
  legend->AddEntry(hist_denom,"EOY","l");
  legend->SetBorderSize(0);
  legend->Draw();

  pad2->cd();


  // plot ratio
  hist_num_clone->Divide(hist_denom);
  hist_num_clone->SetTitle("");
  hist_num_clone->GetXaxis()->SetTitle(label);
  hist_num_clone->GetYaxis()->SetTitle("ratio");
  hist_num_clone->SetLineColor(1);

  // axis styline
  hist_num_clone->GetXaxis()->SetLabelSize(0.1);
  hist_num_clone->GetXaxis()->SetTitleSize(0.15);
  hist_num_clone->GetYaxis()->SetLabelSize(0.1);
  hist_num_clone->GetYaxis()->SetNdivisions(505);
  hist_num_clone->GetYaxis()->SetTitleOffset(0.);
  hist_num_clone->GetYaxis()->SetTitleSize(10);


  hist_num_clone->GetYaxis()->SetRangeUser(0.5,1.5);

  hist_num_clone->Draw("");
  canvas->SaveAs("plots/ULcomparison_"+histname+"_"+beforeFolder+"_"+sample+".pdf");

}
