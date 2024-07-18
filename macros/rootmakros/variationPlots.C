

void drawCMSandLumi(TString year) {

  double lumi = 0;
  if(year == "UL16") lumi = 36.31;
  else if(year == "UL17") lumi = 41.48;
  else if(year == "UL18") lumi = 59.83;
  else lumi = 137.62 ;

  // draw Lumi text
  TString infotext = TString::Format("%3.1f fb^{-1} (%d TeV)", lumi, 13);
  if(lumi > 100.) infotext = TString::Format("%3.0f fb^{-1} (%d TeV)", lumi, 13);
  TLatex *text = new TLatex(3.5, 24, infotext);
  text->SetNDC();
  text->SetTextAlign(33);
  text->SetX(0.93);
  text->SetTextFont(42);
  text->SetY(0.975);
  text->SetTextSize(0.045);
  text->Draw();

  // draw CMS Work in Progress text
  TString cmstext = "CMS";
  TLatex *text2 = new TLatex(3.5, 24, cmstext);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.21);
  text2->SetTextFont(62);
  text2->SetTextSize(0.05);
  text2->SetY(0.9);
  text2->Draw();
  TString preltext = "Work in Progress";
  TLatex *text3 = new TLatex(3.5, 24, preltext);
  text3->SetNDC();
  text3->SetTextAlign(13);
  text3->SetX(0.21);
  text3->SetTextFont(52);
  text3->SetTextSize(0.035);
  text3->SetY(0.85);
  text3->Draw();

}

void variationPlots(){

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

  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.07);

  gROOT->ForceStyle();

  TString sample = "TTbar";
  TString year = "UL18";
  TString channel = "mu";
  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/" + year + "/";
  TString macro_path = "/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/";
  TString filename = "uhh2.AnalysisModuleRunner.MC."+sample+".root";
  TString region = "SignalRegion";
  TString foldername = region + "_" + channel;
  TString histname = "pt_ST";
  bool doJECJER = true;
  bool doPDFandScale = false;

  std::vector<TString> variations_elec = {"sfelec_id", "sfelec_reco", "sfelec_trigger"};
  std::vector<TString> variations_mu = {"sfmu_id", "sfmu_iso", "sfmu_trigger"};
  std::vector<TString> variations_btagging = {"btagging_total", "btagging_hf", "btagging_hfstats1", "btagging_hfstats2", "btagging_lf", "btagging_lfstats1", "btagging_lfstats2", "btagging_cferr1", "btagging_cferr2"};
  std::vector<TString> variations_other = {"pu", "prefiring", "decorrelation"};
  std::vector<int> colors = {1, 2,3,4,6,7,8,9, 11, 12};

  // open main file
  TFile *main_file = TFile::Open(path+"/hadded/"+filename);
  if(!main_file) std::cout << "Main file does not exist" << std::endl;
  TH1D *main_hist = (TH1D*)main_file->Get(foldername+"/"+histname+"_nominal");
  if(!main_hist) std::cout << "Main hist does not exist" << std::endl;

  TCanvas *can = new TCanvas("chist", "c", 600, 600);
  TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  auto line = TLine(600.1,1,6000,1);

  if(doJECJER) {
    // JER
    TFile *file_JER_up = TFile::Open(path+"JER_up/hadded/"+filename);
    if(!file_JER_up) std::cout << "Main file does not exist" << std::endl;
    TH1D *hist_JER_up = (TH1D*)file_JER_up->Get(foldername+"/"+histname+"_nominal");
    if(!hist_JER_up) std::cout << "Main hist does not exist" << std::endl;
    TFile *file_JER_down = TFile::Open(path+"JER_down/hadded/"+filename);
    if(!file_JER_down) std::cout << "Main file does not exist" << std::endl;
    TH1D *hist_JER_down = (TH1D*)file_JER_down->Get(foldername+"/"+histname+"_nominal");
    if(!hist_JER_down) std::cout << "Main hist does not exist" << std::endl;

    hist_JER_up->Divide(main_hist);
    hist_JER_down->Divide(main_hist);

    hist_JER_up->SetLineColor(colors.at(0));
    hist_JER_down->SetLineColor(colors.at(0));
    hist_JER_up->SetLineWidth(3);
    hist_JER_down->SetLineWidth(3);
    hist_JER_up->SetLineStyle(2);
    hist_JER_down->SetLineStyle(3);

    // styling first
    hist_JER_up->GetYaxis()->SetRangeUser(0.7, 1.3);
    hist_JER_up->GetXaxis()->SetRangeUser(600, 6000);
    hist_JER_up->GetXaxis()->SetTitle( hist_JER_up->GetTitle() );
    hist_JER_up->GetYaxis()->SetTitle( "variation / nominal" );
    hist_JER_up->SetTitle("");

    hist_JER_up->Draw("hist");
    hist_JER_down->Draw("hist same");

    leg->AddEntry(hist_JER_up,"JER","l");

    // JEC
    TFile *file_JEC_up = TFile::Open(path+"JEC_up/hadded/"+filename);
    if(!file_JEC_up) std::cout << "Main file does not exist" << std::endl;
    TH1D *hist_JEC_up = (TH1D*)file_JEC_up->Get(foldername+"/"+histname+"_nominal");
    if(!hist_JEC_up) std::cout << "Main hist does not exist" << std::endl;
    TFile *file_JEC_down = TFile::Open(path+"JEC_down/hadded/"+filename);
    if(!file_JEC_down) std::cout << "Main file does not exist" << std::endl;
    TH1D *hist_JEC_down = (TH1D*)file_JEC_down->Get(foldername+"/"+histname+"_nominal");
    if(!hist_JEC_down) std::cout << "Main hist does not exist" << std::endl;

    hist_JEC_up->Divide(main_hist);
    hist_JEC_down->Divide(main_hist);

    hist_JEC_up->SetLineColor(colors.at(1));
    hist_JEC_down->SetLineColor(colors.at(1));
    hist_JEC_up->SetLineWidth(3);
    hist_JEC_down->SetLineWidth(3);
    hist_JEC_up->SetLineStyle(2);
    hist_JEC_down->SetLineStyle(3);

    hist_JEC_up->Draw("hist same");
    hist_JEC_down->Draw("hist same");

    leg->AddEntry(hist_JEC_up,"JEC","l");

    // line
    line = TLine(600.1,1,6000,1);
    line.Draw("same");
    leg->Draw();

    can->SaveAs("plots/variations_" + sample + "_JECJER_" + year + "_" + channel + ".pdf");
  }

  if(doPDFandScale) {

    leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);

    // PDF
    TFile *file_PDF = TFile::Open(macro_path+"/PDF/"+region+"_PDF_"+year+"_" + channel + "_"+sample+".root");
    if(!file_PDF) std::cout << "Main file does not exist" << std::endl;
    TH1D *hist_PDF_up = (TH1D*)file_PDF->Get(sample+"_PDF_up");
    if(!hist_PDF_up) std::cout << "PDF hist up does not exist" << std::endl;
    TH1D *hist_PDF_down = (TH1D*)file_PDF->Get(sample+"_PDF_down");
    if(!hist_PDF_down) std::cout << "PDF hist down does not exist" << std::endl;

    hist_PDF_up->Divide(main_hist);
    hist_PDF_down->Divide(main_hist);

    hist_PDF_up->SetLineColor(colors.at(0));
    hist_PDF_down->SetLineColor(colors.at(0));
    hist_PDF_up->SetLineWidth(3);
    hist_PDF_down->SetLineWidth(3);
    hist_PDF_up->SetLineStyle(2);
    hist_PDF_down->SetLineStyle(3);

    // styling first
    hist_PDF_up->GetYaxis()->SetRangeUser(0.5, 2.1);
    hist_PDF_up->GetXaxis()->SetRangeUser(600, 6000);
    hist_PDF_up->GetXaxis()->SetTitle( hist_PDF_up->GetTitle() );
    hist_PDF_up->GetYaxis()->SetTitle( "variation / nominal" );
    hist_PDF_up->SetTitle("");

    hist_PDF_up->Draw("hist");
    hist_PDF_down->Draw("hist same");

    leg->AddEntry(hist_PDF_up, "PDF", "l");

    // scale
    TFile *file_scale = TFile::Open(macro_path+"/scale/"+region+"_scale_"+year+"_" + channel + "_"+sample+".root");
    if(!file_scale) std::cout << "Main file does not exist" << std::endl;
    TH1D *hist_scale_up = (TH1D*)file_scale->Get(sample+"_scale_up");
    if(!hist_scale_up) std::cout << "scale hist up does not exist" << std::endl;
    TH1D *hist_scale_down = (TH1D*)file_scale->Get(sample+"_scale_down");
    if(!hist_scale_down) std::cout << "scale hist down does not exist" << std::endl;

    hist_scale_up->Divide(main_hist);
    hist_scale_down->Divide(main_hist);

    hist_scale_up->SetLineColor(colors.at(1));
    hist_scale_down->SetLineColor(colors.at(1));
    hist_scale_up->SetLineWidth(3);
    hist_scale_down->SetLineWidth(3);
    hist_scale_up->SetLineStyle(2);
    hist_scale_down->SetLineStyle(3);

    hist_scale_up->Draw("hist same");
    hist_scale_down->Draw("hist same");

    leg->AddEntry(hist_scale_up,"MC scale","l");

    // line
    line = TLine(600.1,1,6000,1);
    line.Draw("same");
    leg->Draw();

    drawCMSandLumi(year);
    can->SaveAs("plots/variations_" + sample + "_PDFandScale_" + year + "_" + channel + ".pdf");
  }

  std::vector<TString> variations;

  // elec
  variations = variations_elec;
  leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  for (int i = 0; i < variations.size(); i++) {
    TString variation = variations.at(i);
    TH1D *hist_up = (TH1D*)main_file->Get(foldername+"/"+histname+"_"+variation+"Up");
    TH1D *hist_down = (TH1D*)main_file->Get(foldername+"/"+histname+"_"+variation+"Down");

    hist_up->Divide(main_hist);
    hist_down->Divide(main_hist);

    hist_up->SetLineColor(colors.at(i));
    hist_down->SetLineColor(colors.at(i));
    hist_up->SetLineWidth(3);
    hist_down->SetLineWidth(3);
    hist_up->SetLineStyle(2);
    hist_down->SetLineStyle(3);

    // styling first
    hist_up->GetYaxis()->SetRangeUser(0.9, 1.1);
    hist_up->GetXaxis()->SetRangeUser(600, 6000);
    hist_up->GetXaxis()->SetTitle( hist_up->GetTitle() );
    hist_up->GetYaxis()->SetTitle( "variation / nominal" );
    hist_up->SetTitle("");

    if(i == 0) hist_up->Draw("hist");
    else hist_up->Draw("hist same");
    hist_down->Draw("hist same");

    leg->AddEntry(hist_up,variation,"l");

  }
  drawCMSandLumi(year);
  line.Draw("same");
  leg->Draw();
  can->SaveAs("plots/variations_" + sample + "_ele_" + year + "_" + channel + ".pdf");

  // mu
  leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  variations = variations_mu;
  for (int i = 0; i < variations.size(); i++) {
    TString variation = variations.at(i);
    TH1D *hist_up = (TH1D*)main_file->Get(foldername+"/"+histname+"_"+variation+"Up");
    TH1D *hist_down = (TH1D*)main_file->Get(foldername+"/"+histname+"_"+variation+"Down");

    hist_up->Divide(main_hist);
    hist_down->Divide(main_hist);

    hist_up->SetLineColor(colors.at(i));
    hist_down->SetLineColor(colors.at(i));
    hist_up->SetLineWidth(3);
    hist_down->SetLineWidth(3);
    hist_up->SetLineStyle(2);
    hist_down->SetLineStyle(3);

    hist_up->GetYaxis()->SetRangeUser(0.95, 1.05);
    hist_up->GetXaxis()->SetRangeUser(600, 6000);
    hist_up->GetXaxis()->SetTitle( hist_up->GetTitle() );
    hist_up->GetYaxis()->SetTitle( "variation / nominal" );
    hist_up->SetTitle("");

    if(i == 0) hist_up->Draw("hist");
    else hist_up->Draw("hist same");
    hist_down->Draw("hist same");

    leg->AddEntry(hist_up,variation,"l");

  }
  drawCMSandLumi(year);
  line.Draw("same");
  leg->Draw();
  can->SaveAs("plots/variations_" + sample + "_mu_" + year + "_" + channel + ".pdf");

  // btagging
  leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  variations = variations_btagging;
  for (int i = 0; i < variations.size(); i++) {
    TString variation = variations.at(i);
    TH1D *hist_up = (TH1D*)main_file->Get(foldername+"/"+histname+"_"+variation+"Up");
    TH1D *hist_down = (TH1D*)main_file->Get(foldername+"/"+histname+"_"+variation+"Down");

    hist_up->Divide(main_hist);
    hist_down->Divide(main_hist);

    hist_up->SetLineColor(colors.at(i));
    hist_down->SetLineColor(colors.at(i));
    hist_up->SetLineWidth(3);
    hist_down->SetLineWidth(3);
    hist_up->SetLineStyle(2);
    hist_down->SetLineStyle(3);

    hist_up->GetYaxis()->SetRangeUser(0., 2);
    hist_up->GetXaxis()->SetRangeUser(600, 6000);
    hist_up->GetXaxis()->SetTitle( hist_up->GetTitle() );
    hist_up->GetYaxis()->SetTitle( "variation / nominal" );
    hist_up->SetTitle("");

    if(i == 0) hist_up->Draw("hist");
    else hist_up->Draw("hist same");
    hist_down->Draw("hist same");

    leg->AddEntry(hist_up,variation,"l");

  }
  drawCMSandLumi(year);
  line.Draw("same");
  leg->Draw();
  can->SaveAs("plots/variations_" + sample + "_btagging_" + year + "_" + channel + ".pdf");

  // other
  leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  variations = variations_other;
  for (int i = 0; i < variations.size(); i++) {
    TString variation = variations.at(i);
    TH1D *hist_up = (TH1D*)main_file->Get(foldername+"/"+histname+"_"+variation+"Up");
    TH1D *hist_down = (TH1D*)main_file->Get(foldername+"/"+histname+"_"+variation+"Down");

    hist_up->Divide(main_hist);
    hist_down->Divide(main_hist);

    hist_up->SetLineColor(colors.at(i));
    hist_down->SetLineColor(colors.at(i));
    hist_up->SetLineWidth(3);
    hist_down->SetLineWidth(3);
    hist_up->SetLineStyle(2);
    hist_down->SetLineStyle(3);

    hist_up->GetYaxis()->SetRangeUser(0.6, 1.4);
    hist_up->GetXaxis()->SetRangeUser(600, 6000);
    hist_up->GetXaxis()->SetTitle( hist_up->GetTitle() );
    hist_up->GetYaxis()->SetTitle( "variation / nominal" );
    hist_up->SetTitle("");

    if(i == 0) hist_up->Draw("hist");
    else hist_up->Draw("hist same");
    hist_down->Draw("hist same");

    leg->AddEntry(hist_up,variation,"l");

  }
  drawCMSandLumi(year);
  line.Draw("same");
  leg->Draw();
  can->SaveAs("plots/variations_" + sample + "_other_" + year + "_" + channel + ".pdf");


}
