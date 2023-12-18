

void checkPUdist() {

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

  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/UL18/hadded/";

  TString sample_MC = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TString sample_data = "uhh2.AnalysisModuleRunner.DATA.DATA.root";

  TString folder = "AfterTrigger";
  TString histname = "N_pv";

  // actual histograms
  TString path_actual = "/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/common/UHH2-data/UL18/";
  TString actual_MC = "MyMCPileupHistogram_UL18.root";
  TString actual_MC_folder = "input_Event";
  TString actual_MC_histname = "N_TrueInteractions";

  TFile *file_MC = TFile::Open(path+sample_MC);
  if(!file_MC) std::cout << "MC file does not exist" << std::endl;
  TFile *file_data = TFile::Open(path+sample_data);
  if(!file_data) std::cout << "data file does not exist" << std::endl;
  TH1D *hist_MC = (TH1D*)file_MC->Get(folder+"/"+histname);
  if(!hist_MC) std::cout << "MC hist does not exist" << std::endl;
  TH1D *hist_data = (TH1D*)file_data->Get(folder+"/"+histname);
  if(!hist_data) std::cout << "data hist does not exist" << std::endl;

  TFile *file_actual_MC = TFile::Open(path_actual+actual_MC);
  if(!file_actual_MC) std::cout << "actual MC file does not exist" << std::endl;
  TH1D *hist_actual_MC = (TH1D*)file_actual_MC->Get(actual_MC_folder+"/"+actual_MC_histname);

  hist_MC->Scale(1/hist_MC->Integral());
  hist_data->Scale(1/hist_data->Integral());
  hist_actual_MC->Scale(1/hist_actual_MC->Integral());

  TCanvas *canvas = new TCanvas("chist", "c", 600, 500);

  hist_MC->SetLineColor(4);
  hist_MC->Draw("hist");
  hist_data->SetLineStyle(2);
  hist_data->Draw("hist same");
  hist_actual_MC->SetLineColor(3);
  hist_actual_MC->Draw("same");

  auto legend = new TLegend(0.55,0.7,0.78,0.88);
  gStyle->SetLegendTextSize(0.05);
  legend->AddEntry(hist_MC,"MC","l");
  legend->AddEntry(hist_data,"data","l");
  legend->AddEntry(hist_actual_MC,"actual MC","l");

  legend->SetBorderSize(0);
  legend->Draw();


  canvas->SaveAs("plots/checkPU_"+histname+"_"+folder+".pdf");

}
