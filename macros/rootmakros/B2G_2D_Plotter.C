// Plot 2D histograms from UHH2
// author: F. Labe
// date: 23.09.2019
// Run it with following command:
// root -l -b -q general2Dplotter.C

void B2G_2D_Plotter(){

  gROOT->SetBatch(kTRUE); // to not open canvas and get XQuartz in the way

  // defining histogram to read
  TString histname = "2D_HT_mjj";


  // input paths, filenames...
  TString path = "/nfs/dust/cms/user/flabe/B2G_Trigger_contact/data/";
  TString fileprefix = "uhh2.AnalysisModuleRunner.";
  TString filename = "DATA.SingleMuon2018_RunD_TrimMassTriggerTest_2018";
  TString signalSampleFilename = "MC.NMSSM_UL18";
  TString before_path = "TrimMass_before";
  TString after_path = "TrimMass_afterAll";
  TString after_path_noTrimMass30 = "TrimMass_TrimMass30removed";
  TString title = "Total";

  TFile *input = TFile::Open(path+fileprefix+filename+".root");
  if(!input) cout << "Empty file" << endl;
  TH2D *hist = (TH2D*)input->Get(before_path+"/"+histname);
  if(!hist) cout << "Empty hist" << endl;

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
  Double_t h = 400;
  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  c1_hist->SetLogz();

  if(histname.Contains("_HT_")) hist->GetXaxis()->SetTitle("H_{T}");
  if(histname.Contains("_mj_")) hist->GetXaxis()->SetTitle("m_{j}");
  hist->GetXaxis()->SetNdivisions(505);
  //hist->GetYaxis()->SetTitle("DNN output");
  hist->SetTitle("");

  hist->Draw("colz");

  // draw Lumi text
  TString infotext = TString::Format("Run 2018D (%d TeV)", 13);
  TLatex *text = new TLatex(3.5, 24, infotext);
  text->SetNDC();
  text->SetTextAlign(33);
  text->SetX(0.88);
  text->SetTextFont(42);
  text->SetY(1);
  text->SetTextSize(0.045);
  text->Draw();

  // draw CMS Work in Progress text
  TString cmstext = "CMS";
  TLatex *text2 = new TLatex(3.5, 24, cmstext);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.185);
  text2->SetTextFont(62);
  text2->SetTextSize(0.05);
  text2->SetY(1);
  text2->Draw();
  TString preltext = "Work in Progress";
  TLatex *text3 = new TLatex(3.5, 24, preltext);
  text3->SetNDC();
  text3->SetTextAlign(13);
  text3->SetX(0.25);
  text3->SetTextFont(52);
  text3->SetTextSize(0.035);
  text3->SetY(0.99);
  text3->Draw();

  // Draw line after certain bin to split presel, sel,
  /**
  Double_t x[2], y[2];
  x[0] = 500;
  x[1] = 500;
  y[0] = 0;
  y[1] = 1;
  TGraph *gr = new TGraph(2,x,y);
  gr->SetLineColor(1);
  gr->SetLineWidth(2);
  gr->Draw("same");
  **/

  c1_hist->SaveAs("plots/B2G_"+histname+".pdf");

  // now the efficiency

  TH2D *hist_eff = (TH2D*)input->Get(after_path+"/"+histname);
  if(!hist_eff) cout << "Empty hist" << endl;

  hist_eff->Divide(hist);

  if(histname.Contains("_HT_")) hist_eff->GetXaxis()->SetTitle("H_{T}");
  if(histname.Contains("_mj_")) hist_eff->GetXaxis()->SetTitle("m_{j}");
  hist_eff->GetXaxis()->SetNdivisions(505);
  //hist->GetYaxis()->SetTitle("DNN output");
  hist_eff->SetTitle("");

  hist_eff->Draw("colz");

  // draw texts
  text->Draw();
  text2->Draw();
  text3->Draw();

  // Draw line after certain bin to split presel, sel,
  /**
  Double_t x[2], y[2];
  x[0] = 500;
  x[1] = 500;
  y[0] = 0;
  y[1] = 1;
  TGraph *gr = new TGraph(2,x,y);
  gr->SetLineColor(1);
  gr->SetLineWidth(2);
  gr->Draw("same");
  **/

  c1_hist->SaveAs("plots/B2G_eff_all_"+histname+".pdf");

  // now the other efficiency

  TH2D *hist_eff_noTrimMass30 = (TH2D*)input->Get(after_path_noTrimMass30+"/"+histname);
  if(!hist_eff_noTrimMass30) cout << "Empty hist" << endl;

  hist_eff_noTrimMass30->Divide(hist);

  if(histname.Contains("_HT_")) hist_eff_noTrimMass30->GetXaxis()->SetTitle("H_{T}");
  if(histname.Contains("_mj_")) hist_eff_noTrimMass30->GetXaxis()->SetTitle("m_{j}");
  hist_eff_noTrimMass30->GetXaxis()->SetNdivisions(505);
  //hist->GetYaxis()->SetTitle("DNN output");
  hist_eff_noTrimMass30->SetTitle("");

  hist_eff_noTrimMass30->Draw("colz");

  // draw texts
  text->Draw();
  text2->Draw();
  text3->Draw();

  // Draw line after certain bin to split presel, sel,
  /**
  Double_t x[2], y[2];
  x[0] = 500;
  x[1] = 500;
  y[0] = 0;
  y[1] = 1;
  TGraph *gr = new TGraph(2,x,y);
  gr->SetLineColor(1);
  gr->SetLineWidth(2);
  gr->Draw("same");
  **/

  c1_hist->SaveAs("plots/B2G_eff_noTrimMass30_"+histname+".pdf");

  // now the difference

  hist_eff->Add(hist_eff_noTrimMass30,-1);
  hist_eff->Draw("colz");

  // draw texts
  text->Draw();
  text2->Draw();
  text3->Draw();

  c1_hist->SaveAs("plots/B2G_gain_TrimMass30_"+histname+".pdf");


}
