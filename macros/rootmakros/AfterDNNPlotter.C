// plot stuff
// author: F.Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q GeneralPlotterSingle.C

void AfterDNNPlotter(){
  
  // config //
  TString subpath = "DNNeval";
  TString histname = "jetmass";

  Double_t w = 400;
  Double_t h = 400;

  // Drawing Definitions 
  TCanvas *canvas = new TCanvas("chist", "c", w, h);
  TLegend *leg = new TLegend(0.58,0.82,0.75,0.94);

  leg->SetTextFont(42);
  leg->SetTextSize(0.035);

  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.16);  gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.05);
  canvas->SetLogy();

  TString source = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/AfterDNN/FullSignal.root";
  
  TFile *input = TFile::Open(source);
  if(!input) cout << "Input file is empty" << endl;
  TH1D *hist1 = (TH1D*)input->Get(subpath+"_lowDNN/"+histname); 
  TH1D *hist2 = (TH1D*)input->Get(subpath+"_highDNN/"+histname); 
  if(!hist1) cout<<"Hist 1 is empty"<<endl;
  if(!hist2) cout<<"Hist 2 is empty"<<endl;
  
  hist1->Scale(1/hist1->Integral());
  hist2->Scale(1/hist2->Integral());

  hist1->SetMaximum(0.8);
  hist1->SetMinimum(1e-3);

  hist1->SetStats(0);
  hist1->GetXaxis()->SetTitleFont(42);
  hist1->GetXaxis()->SetLabelFont(42);
  hist1->GetYaxis()->SetTitleFont(42);
  hist1->GetYaxis()->SetLabelFont(42);

  hist1->GetXaxis()->SetLabelSize(0.055); // 0.045                                                                                
  hist1->GetXaxis()->SetLabelOffset(0.008);
  hist1->GetXaxis()->SetTickLength(0.03);
  hist1->GetXaxis()->SetTitleSize(0.06);   
  hist1->GetXaxis()->SetTitleOffset(1.2);                                                                                       

  hist1->GetYaxis()->SetTitleOffset(1.4); // 1.8                                                                                  
  hist1->GetYaxis()->SetTitleSize(0.06); // 0.05                                                                                  
  hist1->GetYaxis()->SetLabelSize(0.055); // 0.045                                                                                 
  hist1->GetYaxis()->SetTickLength(0.02);
  hist1->GetYaxis()->SetLabelOffset(0.011);

  hist1->GetYaxis()->SetTitle("Events / total events");
  hist1->GetXaxis()->SetTitle(hist1->GetTitle());
  //hist1->GetXaxis()->SetTitle("N_{AK4}"); // for custom
  hist1->GetXaxis()->SetNdivisions(505);
  hist1->SetTitle("");
  hist1->SetMarkerStyle(20);
  hist1->SetMarkerColor(1);
  hist1->SetLineColor(1);
  hist1->SetLineWidth(2);
  hist1->Draw("hist");

  hist2->SetMarkerStyle(20);
  hist2->SetMarkerColor(1);
  hist2->SetLineColor(2);
  hist2->SetLineWidth(2);
  hist2->Draw("hist same");

  leg->AddEntry(hist1, "Below DNN threshold", "l");
  leg->AddEntry(hist2, "Above DNN threshold", "l");
  leg->SetBorderSize(0);
  leg->Draw();

  canvas->SaveAs("AfterDNN_"+subpath+"_"+histname+".pdf");
}
