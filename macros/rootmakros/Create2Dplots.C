// Store 2D histograms filled in UHH2 model as pdf
// author: A.Karavdina
// date: 23.09.2019
// Run it with following command:
// root -l -b -q Create2Dplots.C

//void Create2Dplots(TString filename="uhh2.AnalysisModuleRunner.MC.MC_TstarTstarToTgammaTgamma_M-700_Run2016v3.root", TString subpath="SemiLepTTBarMatch",TString histname="pt_mu_pt_ak8jet1", TString label="TstarTstarToTgammaTgamma_M-700_Run2016v3"){
void Create2Dplots(TString filename="uhh2.AnalysisModuleRunner.MC.MC_TstarTstarToTgammaTgamma_M-700_Run2016v3.root", TString label="TstarTstarToTgammaTgamma_M-700_Run2016v3", TString subpath="After2D",TString histname="pt_mu_pt_ak8jet1"){
  // gStyle->SetOptStat(0);
  // gStyle->SetTitleSize(0.045,"x");  
  // gStyle->SetTitleSize(0.045,"y");
  // gStyle->SetTitleYOffset(0.9);

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
  Double_t w = 600;
  Double_t h = 600;
  
  // bool isWeights = true;
  //  bool isWeights = false;

 //Files after selection
 TString path = "/nfs/dust/cms/user/karavdia/TstarTstar/102X_v1/Preselection/RunII_2016_MuonHihjPtId_cutBasedPhotonIDlooseFall17/semileptonic/";
 // TString gl_label = "QCDHT_slimmedMETsModifiedMET";
 TFile *input = TFile::Open(path+filename);
 TH2D *hist = (TH2D*)input->Get(subpath+"/"+histname);

 TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
 // c1_hist->SetLogz(1);
 hist->GetZaxis()->SetRangeUser(0,5e-4);
 hist->Draw("colz");
 c1_hist->SaveAs(label+"_"+subpath+"_"+histname+".pdf");

}
