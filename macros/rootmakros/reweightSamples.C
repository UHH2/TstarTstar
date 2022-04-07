// A short script to prepare the Root files for the DNN training
// author: F.Labe
// date: 30.06.2010
// Run it with following command:
// root -l -b -q reweightSamples.C

void reweightSamples(){

  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/ModelStudies/";
  TString spin12Name = "uhh2.AnalysisModuleRunner.MC.TstarTstarToTgluonTgluon_top3decays.root";
  TString spin32Name = "uhh2.AnalysisModuleRunner.MC.TstarTstarToTgluonTgluon_Spin32.root";

  TCanvas *canvas = new TCanvas("chist", "c", 800, 600);

  TString plotdir = "plots";
  TString datadir = ".";

  // ##########################################
  // ## Calculating ST ratio for reweighting ##
  // ##########################################

  std::cout << "Start getting pt histogram." << endl;

  // get 12
  TH1D *pt_12;
  {
    TFile *file = new TFile(path+spin12Name, "READ");
    pt_12 = (TH1D*)file->Get("ModelHists/Tstar_pt");
  }

  // get BG ST hist
  TH1D *pt_32;
  {
    TFile *file = new TFile(path+spin32Name, "READ");
    pt_32 = (TH1D*)file->Get("ModelHists/Tstar_pt");
  }

  // Saving the hists.
  TFile *file = new TFile("ModelStudies_pt_weights.root", "RECREATE");
  pt_12->SetName("spin12_Tstar_pt");
  pt_32->SetName("spin32_Tstar_pt");
  pt_12->Write();
  pt_32->Write();

  file->Close();
  delete file;

}