// Efficiency histograms for trigger related studies
// author: A.Karavdina
// date: 24.09.2019
// Run it with following command:
// root -l -b -q TriggerEffPlots.C

void photonEffPlots(TString filename="uhh2.AnalysisModuleRunner.MC.MC_TstarTstarToTgammaTgluon_M-1500_Run2016v3.root", TString label="TstarTstarToTgammaTgluon_M-1500_Run2016v3", TString subpath="SemiLepTTBarMatchGENRECO",TString histname="Pt_photon"){
  
  TString photonID[5] = {"cutBasedPhotonIDlooseFall17", "cutBasedPhotonIDmediumFall17", "cutBasedPhotonIDtightFall17", "mvaPhoIDwp90Fall17", "mvaPhoIDwp80Fall17"};
  
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

  //Files after selection
  //We expect histograms filled with and without trigger selection stored in the same file
  for (int i = 0; i < 5; i++){
    TString path = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Preselection/RunII_2016_MuonHihjPtId_"+photonID[i]+"_nonIsoandIsoHLT_addTTBarRECO/semileptonic/";
    TString path_denom = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Preselection/RunII_2016_MuonHihjPtId_noPhotonID_nonIsoandIsoHLT_addTTBarRECO/semileptonic/";
    TFile *input = TFile::Open(path+filename);
    TH1D *hist_eff = (TH1D*)input->Get(subpath+"/"+histname);//histogram after photon selection
    TFile *input_denom = TFile::Open(path_denom+filename);
    TH1D *hist_denom = (TH1D*)input_denom->Get(subpath+"/"+histname);//histogram before photon selection
    if(!hist_eff || !hist_denom) cout<<"Hists are empty"<<endl;;
    if(!hist_eff || !hist_denom) return;
    if(hist_denom->GetEntries()>0 && hist_eff->GetEntries()>0) hist_eff->Divide(hist_denom);
    TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
    hist_eff->GetXaxis()->SetTitle(hist_eff->GetTitle());
    hist_eff->GetYaxis()->SetTitle("Efficiency");
    hist_eff->GetYaxis()->SetRangeUser(0,1.5);
    hist_eff->SetTitle("");
    hist_eff->SetMarkerStyle(20);
    hist_eff->SetMarkerColor(1);
    hist_eff->Draw();
    c1_hist->SaveAs("PhoIDEff_"+label+"_"+photonID[i]+"_"+subpath+"_"+histname+".pdf");
  }
}
