// Just a basic macro to create simple plots from given hists.
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2019
// Run it with following command:
// root -l -b -q FitDeltaM.C

void FitTopMass(TString subpath="RecoPlots_GEN", TString histname="Mtoplep"){
  
  Double_t w = 800;
  Double_t h = 600;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetHistFillStyle(1);
  gStyle->SetNdivisions(505, "X");

  TString path = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Analysis/";
  TString filename = "FullSignal.root";
  
  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);

  gPad->SetBottomMargin(0.16);
  gPad->SetLeftMargin(0.16);

  TFile *input = TFile::Open(path+filename);
  TH1D *hist = (TH1D*)input->Get(subpath+"/"+histname); //histogram
  if(!hist) cout<<"Hist is empty"<<endl;
  TF1 *fit = new TF1("fit", "gaus", 0, 500);
  hist->Rebin(2);
  hist->Fit(fit, "", "SAME", 125, 220);
  hist->GetXaxis()->SetTitle("m_{t} [GeV]");
  hist->GetYaxis()->SetTitle("Events");
  hist->SetTitle("");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(1);
  hist->SetFillColor(kGray);
  hist->SetLineColor(1);

  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetLabelFont(42);

  hist->GetXaxis()->SetNdivisions(505);
  hist->GetXaxis()->SetLabelSize(0.055); // 0.045                                      
  hist->GetXaxis()->SetLabelOffset(0.01);
  hist->GetXaxis()->SetTickLength(0.03);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.2);

  hist->GetYaxis()->SetTitleOffset(1.2); // 1.8               
  hist->GetYaxis()->SetTitleSize(0.06); // 0.05      
  hist->GetYaxis()->SetLabelSize(0.05); // 0.045                                                                                    
  hist->GetYaxis()->SetTickLength(0.02);
  hist->GetYaxis()->SetLabelOffset(0.011);

  TPaveText *pt = new TPaveText(0.65, 0.65, 0.975, 0.95, "blNDC");
  pt->SetBorderSize(1);
  pt->SetFillColor(0);
  pt->AddText("Fit parameters");
  pt->AddLine(0, 0.75, 1, 0.75);
  pt->AddText((TString)"#chi^{2} / ndf = 63 / 2");
  pt->AddText("mean = 173.8 #pm  0.3");
  pt->AddText("width = 25.6 #pm  0.3");
  for (auto str : {"ndf", "mean", "width"}){
    pt->GetLineWith(str)->SetTextAlign(12);
  }

  hist->Draw("hist");
  pt->Draw();
  fit->Draw("same");

  c1_hist->SaveAs("Fit_"+subpath+"_"+histname+".pdf");
}
