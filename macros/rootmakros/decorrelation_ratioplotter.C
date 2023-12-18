// Just a basic macro to create ratio plot of two hists.
// author: F. Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q ratioplotter.C

void decorrelation_ratioplotter(){

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

  TString sample = "MC.ST";

  bool drop_first_bin_for_norm = false;
  int bin_to_drop = 6; // todo check this

  TString beforePath = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/hadded/";
  TString afterPath = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/hadded/";
  TString filename = "uhh2.AnalysisModuleRunner."+sample+".root";
  TString beforeHistname = "pt_ST_nominal";
  TString afterHistname = "pt_ST_nominal";
  TString label = "S_{T} [GeV]";

  TString beforeFolder = "SignalRegion_total";
  TString afterFolder = "ValidationRegion_total";

  TCanvas *canvas = new TCanvas("canvas", "c", w, h);

  TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.35,1.0,1.0);
  TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.35);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  pad1->SetLogy();

  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);

  pad1->SetTickx();
  pad1->SetTicky();

  TFile *input_EOY = TFile::Open(beforePath+filename);
  TFile *input_UL = TFile::Open(afterPath+filename);

  TH1D *hist_num = (TH1D*)input_UL->Get(beforeFolder+"/"+beforeHistname); //histogram
  TH1D *hist_denom = (TH1D*)input_EOY->Get(afterFolder+"/"+afterHistname); //histogram
  if(!hist_num) std::cout << "Numerator does not exist!" << std::endl;
  if(!hist_denom) std::cout << "Denominator does not exist!" << std::endl;

  double scaling_factor_num = hist_num->Integral();
  double scaling_factor_denom = hist_denom->Integral();

  if(drop_first_bin_for_norm) {

    std::cout << hist_num->GetBinContent(bin_to_drop) << " " << hist_denom->GetBinContent(bin_to_drop) << std::endl;

    scaling_factor_num -= hist_num->GetBinContent(bin_to_drop);
    scaling_factor_denom -= hist_denom->GetBinContent(bin_to_drop);
  } 

  hist_num->Scale(1./scaling_factor_num);
  hist_denom->Scale(1./scaling_factor_denom);
  
  TH1D* hist_num_clone = (TH1D*) hist_num->Clone();

  hist_num->SetTitle("");
  hist_num->GetXaxis()->SetTickLength(0.02);
  hist_num->GetXaxis()->SetLabelOffset(999);
  hist_num->GetYaxis()->SetTitleOffset(0.75);
  hist_num->GetYaxis()->SetLabelSize(0.075);
  hist_num->GetYaxis()->SetTitleSize(0.075);
  hist_num->GetYaxis()->SetTitle("Events");
  hist_num->GetYaxis()->SetRangeUser(2e-7, 20);
  hist_num->GetYaxis()->SetTickLength(0.02);
  hist_num->SetMarkerStyle(20);
  hist_num->SetMarkerColor(1);
  hist_num->SetLineColor(1);
  hist_num->Draw("hist");
  hist_denom->SetMarkerStyle(20);
  hist_denom->SetMarkerColor(2);
  hist_denom->SetLineColor(2);
  hist_denom->Draw("hist same");

  auto legend = new TLegend(0.62,0.7,0.85,0.88);
  gStyle->SetLegendTextSize(0.05);
  legend->AddEntry(hist_num,"Signal region","l");
  legend->AddEntry(hist_denom,"Validation region","l");
  legend->SetBorderSize(0);
  legend->Draw();

  // draw CMS Work in Progress text
  TString cmstext = "CMS";
  TLatex *text2 = new TLatex(3.5, 24, cmstext);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.16);
  text2->SetTextFont(62);
  text2->SetTextSize(0.07);
  text2->SetY(0.3+0.55);
  text2->Draw();
  TString preltext = "Work in Progress";
  TLatex *text3 = new TLatex(3.5, 24, preltext);
  text3->SetNDC();
  text3->SetTextAlign(13);
  text3->SetX(0.16);
  text3->SetTextFont(52);
  text3->SetTextSize(0.055);
  text3->SetY(0.225+0.55);
  text3->Draw();

  TLatex *text4 = new TLatex(3.5, 24, "Simulation");
  text4->SetNDC();
  text4->SetTextAlign(13);
  text4->SetX(0.23);
  text4->SetTextFont(52);
  text4->SetTextSize(0.055);
  text4->SetY(0.225+0.55+0.0675);
  text4->Draw();


  pad2->cd();
  pad2->SetTickx();
  pad2->SetTicky();

  // plot ratio
  hist_num_clone->Divide(hist_denom);
  hist_num_clone->SetTitle("");
  hist_num_clone->GetXaxis()->SetTitle(label);
  hist_num_clone->GetYaxis()->SetTitle("Ratio");
  hist_num_clone->SetLineColor(2);

  // axis styline
  hist_num_clone->GetXaxis()->SetLabelSize(0.13);
  hist_num_clone->GetXaxis()->SetTitleSize(0.15);
  hist_num_clone->GetXaxis()->SetLabelOffset(0.01);

  hist_num_clone->GetYaxis()->SetLabelSize(0.13);
  hist_num_clone->GetYaxis()->SetLabelOffset(0.01);
  hist_num_clone->GetYaxis()->SetNdivisions(505);
  hist_num_clone->GetYaxis()->SetTitleOffset(0.375);
  hist_num_clone->GetYaxis()->SetTitleSize(0.15);

  hist_num_clone->GetYaxis()->SetRangeUser(0.5, 2);

  //gPad->SetLogy();

  hist_num_clone->Draw("");

  // storing this histogram to a file to use as uncertainty later
  TFile *output = TFile::Open("files/decorrelationComparison_" + sample + ".root", "recreate");
  hist_num_clone->SetName("decorrelation_uncertainty");
  hist_num_clone->Write();
  output->Close();

  if(drop_first_bin_for_norm) canvas->SaveAs("plots/decorrelationComparison_"+beforeHistname+"_"+beforeFolder+"_"+sample+"_firstDropped.pdf");
  else canvas->SaveAs("plots/decorrelationComparison_"+beforeHistname+"_"+beforeFolder+"_"+sample+".pdf");

}
