// Efficiency histograms for trigger related studies
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2019
// Run it with following command:
// root -l -b -q photonSelecPlots.C

void photonSelecPlots(TString filename="uhh2.AnalysisModuleRunner.MC.MC_TstarTstarToTgammaTgluon_M-1500_Run2016v3.root", TString label="TstarTstarToTgammaTgluon_M-1500_Run2016v3", TString subpath="After_TstarTstar_Reco",TString histname="N_photon"){
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

  //TODO: should use some kind of for loop if this is actually used in any way!
  //Files after selection
  
  TString crit[5] = {"mvaPhoIDwp90"};

  TString path_pre = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Preselection/RunII_2016_MuonHihjPtId_";
  TString path_post = "Fall17_nonIsoandIsoHLT_addTTBarRECO/semileptonic/";
  
  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  TLegend *leg = new TLegend(0.5,0.7,0.83,0.90);

  for (int i = 0; i<1; i++){
    TFile *input = TFile::Open(path_pre+crit[i]+path_post+filename);
    TH1D *hist = (TH1D*)input->Get(subpath+"/"+histname);//histogram
    if(!hist) cout<<"Hist is empty"<<endl;;

    leg->AddEntry(hist,crit[i],"pl");    
    hist->SetTitle("");
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(i+1);
    if (i == 0) {
      hist->GetXaxis()->SetTitle("N_{#gamma}");
      hist->Draw();
    }
    else {
      hist->GetXaxis()->SetTitle("N_{#gamma}");
      hist->Draw("same");
    }
  }
  
  //leg->Draw();


  c1_hist->SaveAs("N_Gamma_"+label+"_"+subpath+"_"+histname+".pdf");
}
