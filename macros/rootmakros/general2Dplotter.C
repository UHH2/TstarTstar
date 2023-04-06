// Plot 2D histograms from UHH2
// author: F. Labe
// date: 23.09.2019
// Run it with following command:
// root -l -b -q general2Dplotter.C

void general2Dplotter(TString filename="TstarTstar", TString subpath="DNN_newTagger", TString histname="2D_newTagger_ST"){
  // gStyle->SetOptStat(0);
  // gStyle->SetTitleSize(0.045,"x");
  // gStyle->SetTitleSize(0.045,"y");
  // gStyle->SetTitleYOffset(0.9);

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
 TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/hadded/";
 //TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/TriggerSF/";
 TString fileprefix = "uhh2.AnalysisModuleRunner.MC.";
 //TString fileprefix = "";

 TFile *input = TFile::Open(path+fileprefix+filename+".root");
 if(!input) cout << "Empty file" << endl;
 TH2D *hist = (TH2D*)input->Get(subpath+"/"+histname);
 if(!hist) cout << "Empty hist" << endl;

 TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
 c1_hist->SetLogz();
 //hist->GetZaxis()->SetRangeUser(0,1);

 // scale some stuff
 /**
 std::cout << "Attention! Reweighting by row is done!!" << std::endl;
 int binsX = 50;
 int binsY = 50;

 for(uint x = 1; x <= binsX; x++){
   double sum = 0;
   for(uint y = 1; y <= binsY; y++){
     sum += hist->GetBinContent(x, y);
   }
   cout << x << "   " << sum << endl;
   if(sum > 0){
     for(uint y = 1; y <= binsY; y++){
       hist->SetBinContent(x, y, hist->GetBinContent(x, y)/sum);
     }
   }
 }
 **/

 //hist->GetXaxis()->SetTitle("S_{T} [GeV]");
 hist->GetXaxis()->SetNdivisions(505);
 //hist->GetYaxis()->SetTitle("DNN output");
 //hist->GetZaxis()->SetRangeUser(10, 4e2);
 hist->SetTitle("");

 hist->Draw("colz");

 // draw Lumi text
 TString infotext = TString::Format("%3.1f fb^{-1} (%d TeV)", 137., 13);
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
 text3->SetX(0.29);
 text3->SetTextFont(52);
 text3->SetTextSize(0.035);
 text3->SetY(0.986);
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

 c1_hist->SaveAs("plots/"+filename+"_"+subpath+"_"+histname+".pdf");

}
