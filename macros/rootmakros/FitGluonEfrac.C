// Just a basic macro to create simple plots from given hists.
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2019
// Run it with following command:
// root -l -b -q FitDeltaM.C

void FitGluonEfrac(){
  
  Double_t w = 800;
  Double_t h = 600;

  gStyle->SetOptStat(0);
  gStyle->SetHistFillStyle(1);
  gStyle->SetNdivisions(505, "XY");
  gStyle->SetOptFit(0);

  TString path = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/GenInfo/";
  TString filename = "MC.FullSignal.root";
  TString subpath = "GenHists";
  
  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);

  gPad->SetBottomMargin(0.16);
  gPad->SetLeftMargin(0.16);

  std::cout << "Doing AK4" << endl;
  TFile *input = TFile::Open(path+filename);
  TH1D *hist = (TH1D*)input->Get(subpath+"/gluon_Efrac_AK4"); //histogram
  if(!hist) cout<<"Hist is empty"<<endl;
  TF1 *fit = new TF1("fit", "gaus", 0, 2);
  hist->Fit(fit, "", "SAME", 0, 2);
  hist->GetXaxis()->SetTitle("E_{jet}/E_{gluon}");
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
  pt->AddText((TString)"#chi^{2} / ndf = 2.73 #times 10^{5} / 17");
  pt->AddText("mean = 0.85 #pm  0.00");
  pt->AddText("width = 0.16 #pm  0.00");
  for (auto str : {"ndf", "mean", "width"}){
    pt->GetLineWith(str)->SetTextAlign(12);
  }

  hist->Draw("hist");
  pt->Draw();
  fit->Draw("same");

  c1_hist->SaveAs("Gluon_Efrac_AK4.pdf");

  std::cout << "Doing HOTVR" << endl;
  TH1D *hist2 = (TH1D*)input->Get(subpath+"/"+"gluon_Efrac_HOTVR"); //histogram
  if(!hist2) cout<<"Hist is empty"<<endl;
  TF1 *fit2 = new TF1("fit", "gaus", 0, 2);
  hist2->Fit(fit2, "", "SAME", 0, 2);
  hist2->GetXaxis()->SetTitle("E_{jet}/E_{gluon}");
  hist2->GetYaxis()->SetTitle("Events");
  hist2->SetTitle("");
  hist2->SetMarkerStyle(20);
  hist2->SetMarkerColor(1);
  hist2->SetFillColor(kGray);
  hist2->SetLineColor(1);

  hist2->GetXaxis()->SetTitleFont(42);
  hist2->GetXaxis()->SetLabelFont(42);
  hist2->GetYaxis()->SetTitleFont(42);
  hist2->GetYaxis()->SetLabelFont(42);

  hist2->GetXaxis()->SetLabelSize(0.055); // 0.045                                                                                
  hist2->GetXaxis()->SetLabelOffset(0.01);
  hist2->GetXaxis()->SetTickLength(0.03);
  hist2->GetXaxis()->SetTitleSize(0.05);                                                                                        
  hist2->GetXaxis()->SetTitleOffset(1.2);                                                                                       

  hist2->GetYaxis()->SetTitleOffset(1.2); // 1.8                                                                                  
  hist2->GetYaxis()->SetTitleSize(0.06); // 0.05                                                                                  
  hist2->GetYaxis()->SetLabelSize(0.05); // 0.045                                                                                 
  hist2->GetYaxis()->SetTickLength(0.02);
  hist2->GetYaxis()->SetLabelOffset(0.011);

  TPaveText *pt2 = new TPaveText(0.65, 0.65, 0.975, 0.95, "blNDC");
  pt2->SetBorderSize(1);
  pt2->SetFillColor(0);
  pt2->AddText("Fit parameters");
  pt2->AddLine(0, 0.75, 1, 0.75);
  pt2->AddText((TString)"#chi^{2} / ndf = 1.72 #times 10^{5} / 17");
  pt2->AddText("mean = 0.87 #pm  0.00");
  pt2->AddText("width = 0.14 #pm  0.00");
  for (auto str : {"ndf", "mean", "width"}){
    pt2->GetLineWith(str)->SetTextAlign(12);
  }


  hist2->Draw("hist");
  pt2->Draw();
  fit2->Draw("same");
  c1_hist->SaveAs("Gluon_Efrac_HOTVR.pdf");
}
