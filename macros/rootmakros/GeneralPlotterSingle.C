// Just a basic macro to create simple plots from given hists.
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2019
// Run it with following command:
// root -l -b -q GeneralPlotterSingle.C

void GeneralPlotterSingle(TString channel = "tgtg", TString masspoint="1600", TString subpath="NoCuts_GEN", TString histname="Pt_tstar_gen"){

  Double_t w = 800;
  Double_t h = 600;

  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/";
  TString filename = "uhh2.AnalysisModuleRunner.MC.MC_TstarTstarTo";
  if(channel == "tgtg"){filename += "TgluonTgluon";}
  else if(channel == "tgtgamma"){filename += "TgammaTgluon";}
  else {cout<<"Channel invalid!"<<endl;}
  filename += "_M-"+masspoint+"_Run2016v3.root";

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);

  TFile *input = TFile::Open(path+filename);
  TH1D *hist = (TH1D*)input->Get(subpath+"/"+histname); //histogram
  if(!hist) cout<<"Hist is empty"<<endl;;
  hist->GetXaxis()->SetTitle("M_{T*} [GeV]");
  hist->GetYaxis()->SetTitle("events");
  hist->SetTitle(channel + ": " + histname + " (" + subpath + ")");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(1);
  hist->SetLineColor(1);
  hist->Draw("hist");

  c1_hist->SaveAs(channel+"_M-"+masspoint+"_"+subpath+"_"+histname+".pdf");
}
