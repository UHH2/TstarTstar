// Just a basic macro to create ratio plot of two hists.
// author: F. Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q ratioplotter.C

void ratioplotter_2(){

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

  TString path_1 = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/hadded/";
  TString path_2 = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/currentMAIN/";
  TString filename = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TString histname = "pt_ST";
  TString label = "S_{T}";
  TString outputfilename = "ST";

  TString folder_num = "AfterDNNcut_06";
  TString folder_denom  = "notDNNcut_03";

  TCanvas *canvas = new TCanvas("canvas", "c", w, h);
  canvas->SetLogy();

  // new
  TH1D *ratio_new;
  {
    TFile *input = TFile::Open(path_1+filename);
    TH1D *hist_num = (TH1D*)input->Get(folder_num+"/"+histname); //histogram
    TH1D *hist_denom = (TH1D*)input->Get(folder_denom+"/"+histname); //histogram
    if(!hist_num) std::cout << "Numerator does not exist!" << std::endl;
    if(!hist_denom) std::cout << "Denominator does not exist!" << std::endl;

    // normalize both hists
    hist_num->Scale(1/hist_num->Integral());
    hist_denom->Scale(1/hist_denom->Integral());

    hist_num->SetTitle("Flat reweighting");
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

    auto legend = new TLegend(0.55,0.7,0.78,0.88);
    gStyle->SetLegendTextSize(0.05);
    legend->AddEntry(hist_num,"Above 0.6","l");
    legend->AddEntry(hist_denom,"Below 0.3","l");
    legend->Draw();

    canvas->SaveAs("plots/comparison_new_"+outputfilename+".pdf");

    // plot ratio
    hist_num->Divide(hist_denom);
    hist_num->SetTitle("Flat reweighting");
    hist_num->GetYaxis()->SetTitle("ratio");
    hist_num->Draw("hist");
    canvas->SaveAs("plots/ratio_new_"+outputfilename+".pdf");

    ratio_new = hist_num;
  }

  // previous
  TH1D *ratio_previous;
  {
    TFile *input = TFile::Open(path_2+filename);
    TH1D *hist_num = (TH1D*)input->Get(folder_num+"/"+histname); //histogram
    TH1D *hist_denom = (TH1D*)input->Get(folder_denom+"/"+histname); //histogram
    if(!hist_num) std::cout << "Numerator does not exist!" << std::endl;
    if(!hist_denom) std::cout << "Denominator does not exist!" << std::endl;

    // normalize both hists
    hist_num->Scale(1/hist_num->Integral());
    hist_denom->Scale(1/hist_denom->Integral());

    hist_num->SetTitle("Previous reweighting");
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

    auto legend = new TLegend(0.55,0.7,0.78,0.88);
    gStyle->SetLegendTextSize(0.05);
    legend->AddEntry(hist_num,"Above 0.6","l");
    legend->AddEntry(hist_denom,"Below 0.3","l");
    legend->Draw();

    canvas->SaveAs("plots/comparison_old_"+outputfilename+".pdf");

    // plot ratio
    hist_num->Divide(hist_denom);
    hist_num->SetTitle("Previous reweighting");
    hist_num->GetYaxis()->SetTitle("ratio");
    hist_num->Draw("hist");
    canvas->SaveAs("plots/ratio_old_"+outputfilename+".pdf");

    ratio_previous = hist_num;
  }

  ratio_previous->Draw("hist");
  ratio_new->Draw("hist same");
  ratio_new->SetMarkerStyle(20);
  ratio_new->SetMarkerColor(2);
  ratio_new->SetLineColor(2);

  auto legend = new TLegend(0.55,0.7,0.78,0.88);
  gStyle->SetLegendTextSize(0.05);
  legend->AddEntry(ratio_previous,"Old model","l");
  legend->AddEntry(ratio_new,"Flat reweighting","l");
  legend->Draw();

  canvas->SaveAs("plots/both_ratios_"+outputfilename+".pdf");

  ratio_new->Divide(ratio_previous);
  ratio_new->GetYaxis()->SetTitle("ratio of ratios");
  ratio_new->Draw("hist");
  canvas->SaveAs("plots/ratio_of_ratios_"+outputfilename+".pdf");

}
