// Efficiency histograms for trigger related studies
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2019
// Run it with following command:
// root -l -b -q photonPTPlots.C

void photonPTPlots(TString filename="uhh2.AnalysisModuleRunner.MC.MC_TstarTstarToTgammaTgluon_M-1500_Run2016v3.root", TString label="TstarTstarToTgammaTgluon_M-1500_Run2016v3", TString subpath="After_TstarTstar_Reco",TString histname="pt_photon"){
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
  Double_t w = 1000;
  Double_t h = 600;

  TString path = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Preselection/RunII_2016_MuonHihjPtId_mvaPhoIDwp90Fall17_nonIsoandIsoHLT_addTTBarRECO/semileptonic/";

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  TLegend *leg = new TLegend(0.5,0.7,0.83,0.90);

  TString hists[] = {"pt_photon", "pt_photon_1", "pt_photon_2"};

  int i = 1;
  for(auto histname: hists){
    TFile *input = TFile::Open(path+filename);
    TH1D *hist = (TH1D*)input->Get(subpath+"/"+histname);//histogram
    if(!hist) cout<<"Hist is empty"<<endl;;

    leg->AddEntry(hist, histname, "l");    
    hist->SetTitle("Photon PT");
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(i);
    hist->SetLineColor(i);
    hist->GetXaxis()->SetTitle("PT_{#gamma}");
    hist->Draw("hist same");
    i++;
  }
  
  leg->Draw();


  c1_hist->SaveAs("PhoPT_"+label+"_"+subpath+"_"+histname+".pdf");
}
