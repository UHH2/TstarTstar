// Just a basic macro to create ratio plot of two hists.
// author: F. Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q ratioplotter.C

void ratioplotter(){

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

  TString sample = "TTbar";

  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/UL16postVFP/hadded/";
  TString filename = "uhh2.AnalysisModuleRunner.MC."+sample+".root";
  TString histname = "pt_ST";
  TString label = "S_{T}";

  //TString folder_num = "AfterDNNcut_06";
  //TString folder_denom  = "notDNNcut_06";
  //TString folder_denom_2  = "notDNNcut_03";
  TString suffix = "";
  TString folder_num = "newTaggerSR";
  TString folder_denom  = "newTaggerCR";
  TString folder_denom_2  = "newTaggerCR";

  TCanvas *canvas = new TCanvas("canvas", "c", w, h);
  canvas->SetLogy();

  TFile *input = TFile::Open(path+filename);
  TH1D *hist_num = (TH1D*)input->Get(folder_num+suffix+"/"+histname); //histogram
  TH1D *hist_denom = (TH1D*)input->Get(folder_denom+suffix+"/"+histname); //histogram
  TH1D *hist_denom_2 = (TH1D*)input->Get(folder_denom_2+suffix+"/"+histname); //histogram
  if(!hist_num) std::cout << "Numerator does not exist!" << std::endl;
  if(!hist_denom) std::cout << "Denominator does not exist!" << std::endl;

  // normalize both hists
  hist_num->Scale(1/hist_num->Integral());
  hist_denom->Scale(1/hist_denom->Integral());
  hist_denom_2->Scale(1/hist_denom_2->Integral());

  hist_num->SetTitle("");
  hist_num->GetXaxis()->SetTitle(label);
  hist_num->GetYaxis()->SetTitle("events [a.u.]");
  hist_num->SetMarkerStyle(20);
  hist_num->SetMarkerColor(1);
  hist_num->SetLineColor(1);
  hist_num->Draw("hist");
  hist_denom->SetMarkerStyle(20);
  hist_denom->SetMarkerColor(2);
  hist_denom->SetLineColor(2);
  hist_denom->Draw("hist same");
  hist_denom_2->SetMarkerStyle(20);
  hist_denom_2->SetMarkerColor(4);
  hist_denom_2->SetLineColor(4);
  //hist_denom_2->Draw("hist same");

  auto legend = new TLegend(0.55,0.7,0.78,0.88);
  gStyle->SetLegendTextSize(0.05);
  legend->AddEntry(hist_num,"SR","l");
  legend->AddEntry(hist_denom,"CR","l");
  //legend->AddEntry(hist_denom_2,"Same CR","l");
  legend->Draw();

  canvas->SaveAs("plots/STcomparison_SRCR_"+sample+suffix+".pdf");

  canvas->SetLogy(false);

  // plot ratio
  TH1D *hist_num_2 = (TH1D*) hist_num->Clone();
  hist_num->Divide(hist_denom);
  hist_num_2->Divide(hist_denom_2);
  hist_num->GetYaxis()->SetTitle("ratio");
  hist_num->SetLineColor(2);
  hist_num_2->SetLineColor(4);
  hist_num->Draw("hist");
  //hist_num_2->Draw("hist same");
  canvas->SaveAs("plots/STratio_SRCR_"+sample+suffix+".pdf");


}
