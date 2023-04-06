// Decorrelate DNN tagger, see https://confluence.desy.de/x/0up4DQ
// author: F.Labe
// date: 24.09.2021
// Run it with following command:
// root -l -b -q decorrelatedTagger.C

void evaluateTagger(){

  Double_t cut_value = 0;

  // some style options
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
  Double_t w = 800;
  Double_t h = 600;

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  c1_hist->SetLogz();

  // get histogram
  TString filename="TstarTstar_M-1800";
  TString subpath="DNN_newTagger";
  TString histname="2D_newTagger_ST";
  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/UL16postVFP/hadded/";
  TString fileprefix = "uhh2.AnalysisModuleRunner.MC.";
  TFile *input = TFile::Open(path+fileprefix+filename+".root");
  if(!input) cout << "Empty file" << endl;
  TH2D *hist = (TH2D*)input->Get(subpath+"/"+histname);
  if(!hist) cout << "Empty hist" << endl;

  // flip
  TH2D *hist2 = new TH2D("oldrebin",hist->GetTitle(), 40, 0, 4000, 50, -1, 1);
  TAxis *xaxis = hist->GetXaxis(); TAxis *yaxis = hist->GetYaxis();
  for (int j=1; j<=yaxis->GetNbins();j++) {
    for (int i=1; i<=xaxis->GetNbins();i++) {
      hist2->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(yaxis->GetNbins()-j),hist->GetBinContent(i,j));
    }
  }
  hist = hist2;

  // calculate DNN cut value for each bin in x
  int binsX = 100;
  TH1D *efficiency = new TH1D("efficiency", "efficiency", 40, 0, 4000);

  for(uint x = 0; x <= binsX; x++){
    TH1D* projection =  hist->ProjectionY("proj", x, x);
    Double_t total = projection->Integral();
    int bin = projection->FindBin(cut_value);
    Double_t partial = projection->Integral(0, bin);
    if(partial > 0 && total > 0) efficiency->Fill(efficiency->GetBinCenter(x), partial/total);
  }

  efficiency->SetTitle("");
  efficiency->GetXaxis()->SetTitle("S_{T} [GeV]");
  efficiency->GetYaxis()->SetTitle("efficiency");
  efficiency->GetXaxis()->SetNdivisions(505);
  efficiency->GetXaxis()->SetRangeUser(500, 4000);
  efficiency->GetYaxis()->SetRangeUser(0, 1);
  efficiency->Draw("hist");
  c1_hist->Update();

  Double_t mean = efficiency->Integral()/40.;
  std::cout << "Mean: " << mean << std::endl;

  TLine *meanLine = new TLine(0, mean, 4000, mean);
  //meanLine->Draw("same");

  // draw Lumi text
  TString infotext = TString::Format("%3.1f fb^{-1} (%d TeV)", 137., 13);
  TLatex *text = new TLatex(3.5, 24, infotext);
  text->SetNDC();
  text->SetTextAlign(33);
  text->SetX(0.82);
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
  text3->SetX(0.27);
  text3->SetTextFont(52);
  text3->SetTextSize(0.035);
  text3->SetY(0.986);
  text3->Draw();

  c1_hist->SaveAs("plots/effi_"+filename+"_"+histname+".pdf");


}
