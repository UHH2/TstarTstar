// Draw hatched line
// author: F.Labe
// date: 02.02.2021
// Run it with following command:
// root -l -q hatchedLine.C

void hatchedLine(){

  gROOT->SetBatch(kTRUE); // to not open canvas and get XQuartz in the way

  Double_t w = 420;
  Double_t h = 400;

  // Drawing Definitions
  TCanvas *canvas = new TCanvas("chist", "c", w, h);


  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.16);  gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.065);
  gPad->SetTicky();

  canvas->SetLogy();

  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/hadded/";
  TString filename = "uhh2.AnalysisModuleRunner.MC.TstarTstar_M-1000.root";
  TString subpath = "DNN";
  TString histname = "DNN_output";

  TFile *input = TFile::Open(path+filename);
  TH1D *hist = (TH1D*)input->Get(subpath+"/"+histname); //histogram
  if(!hist) cout<<"Hist is empty"<<endl;;
  hist->GetXaxis()->SetTitle("S_{T} [GeV]");
  hist->GetYaxis()->SetTitle("Events");
  hist->SetTitle("");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(0);
  hist->SetLineColor(0);

  hist->SetStats(0);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetAxisColor(0);
  hist->GetYaxis()->SetLabelColor(0);
  hist->GetYaxis()->SetTitleColor(0);

  hist->GetXaxis()->SetLabelSize(0.055); // 0.045
  hist->GetXaxis()->SetLabelOffset(0.01);
  hist->GetXaxis()->SetTickLength(0.03);
  //hist->GetXaxis()->SetTitleSize(0.05);
  //hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetXaxis()->SetAxisColor(0);
  hist->GetXaxis()->SetTitleColor(0);
  hist->GetXaxis()->SetLabelColor(0);


  hist->GetYaxis()->SetTitleOffset(1.5); // 1.8
  hist->GetYaxis()->SetTitleSize(0.06); // 0.05
  hist->GetYaxis()->SetLabelSize(0.05); // 0.045
  hist->GetYaxis()->SetTickLength(0.02);
  hist->GetYaxis()->SetLabelOffset(0.011);

  hist->Draw("hist");

  // Draw line after certain bin to split presel, sel,
  Double_t x[2], y[2];
  x[0] = 0.7;
  x[1] = 0.7;
  y[0] = 0;
  y[1] = 1000;
  TGraph *gr = new TGraph(2,x,y);
  gr->SetLineColor(1);
  gr->SetLineWidth(-202);
  gr->SetFillStyle(3004);
  gr->SetFillColor(1);
  gr->Draw("same");

  canvas->SaveAs("plots/hatchedLine.pdf");


}
