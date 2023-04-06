// this makro takes the "before" hist (from before_path) and the "after Trigger" hists (afterPaths) and calculates efficiencies by deviding.
// Everything else is just input reading and a lot of styling options

void B2G_TrimMass_Efficiency(){

  gROOT->SetBatch(kTRUE); // to not open canvas and get XQuartz in the way

  // defining histogram to read
  TString histname = "lep_pt_aboveHT_600";
  bool doSignalSample = false;
  bool doComparison = false; // only one can be true, doComparison is stronger

  // input paths, filenames...
  TString path = "/nfs/dust/cms/user/flabe/B2G_Trigger_contact/data/";
  TString fileprefix = "uhh2.AnalysisModuleRunner.";
  TString filename = "DATA.SingleMuon2018_RunD_TrimMassTriggerTest_2018";
  TString signalSampleFilename = "MC.NMSSM_UL18";
  TString before_path = "LeptonCross_before";
  std::vector<TString> afterPaths = {"LeptonCross_afterAll", "LeptonCross_afterWithoutCross"};
  std::vector<TString> titles = {"All", "Without cross-trigger"};
  std::vector<int> colors = {kGreen+2,kRed};

  // set x label
  TString xAxisLabel;
  xAxisLabel = "p^{e}_{T} [GeV]";

  // open file
  TFile *file = TFile::Open(path+fileprefix+filename+".root");

  // doing the quick alpha check
  /**
  TH1D *hist_alpha_ref = (TH1D*)file->Get("h_alpha_ref/oneBin");
  TH1D *hist_alpha_mu = (TH1D*)file->Get("h_alpha_mu/oneBin");
  TH1D *hist_alpha_jets = (TH1D*)file->Get("h_alpha_jets/oneBin");
  TH1D *hist_alpha_both= (TH1D*)file->Get("h_alpha_both/oneBin");

  double eff_mu = hist_alpha_mu->Integral() / hist_alpha_ref->Integral();
  double eff_jets = hist_alpha_jets->Integral() / hist_alpha_ref->Integral();
  double eff_both = hist_alpha_both->Integral() / hist_alpha_ref->Integral();
  std::cout << "Efficiencies: " << std::endl;
  std::cout << "Single muon trigger (reference): " << eff_mu << std::endl;
  std::cout << "Combined jets triggers: " << eff_jets << std::endl;
  std::cout << "Both together: " << eff_both << std::endl;

  double alpha = (eff_mu * eff_jets)/eff_both;
  std::cout << "The alpha parameter is " << alpha << "." << std::endl;
  **/

  // read before hist histograms
  TH1D *hist_before = (TH1D*)file->Get(before_path+"/"+histname);

  // plotting parameters
  Double_t w = 600;
  Double_t h = 400;
  TCanvas *canvas = new TCanvas("chist", "c", w, h);

  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.16);  gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.09);
  gPad->SetTicky();

  TLegend *leg = new TLegend(0.6, 0.75, 0.9, 0.95);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);

  // calculate and plot efficiency hists
  std::vector<TGraphAsymmErrors> effs = {TGraphAsymmErrors(), TGraphAsymmErrors(), TGraphAsymmErrors(), TGraphAsymmErrors(), TGraphAsymmErrors()};
  int i = 0;
  for (const auto & subpath : afterPaths) {
    std::cout << "Processing " << subpath << "." << std::endl;
    TH1D *hist = (TH1D*)file->Get(subpath+"/"+histname);

    // this line is the actual efficiency calculation!
    effs.at(i).Divide(hist, hist_before, "cl=0.68 b(1,1) mode");

    effs.at(i).GetXaxis()->SetTitle(xAxisLabel);
    effs.at(i).GetYaxis()->SetTitle("Efficiency");
    effs.at(i).SetTitle("");
    effs.at(i).SetMarkerStyle(1);
    effs.at(i).SetLineWidth(2);
    effs.at(i).SetLineColor(colors.at(i));
    effs.at(i).SetMarkerColor(colors.at(i));

    // Drawing
    if(i == 0) {
      // cosmetics
      effs.at(i).GetXaxis()->SetTitleFont(42);
      effs.at(i).GetXaxis()->SetLabelFont(42);
      effs.at(i).GetYaxis()->SetTitleFont(42);
      effs.at(i).GetYaxis()->SetLabelFont(42);

      effs.at(i).GetXaxis()->SetLabelSize(0.055); // changed from 0.045
      effs.at(i).GetXaxis()->SetLabelOffset(0.008);
      effs.at(i).GetXaxis()->SetTickLength(0.03);
      effs.at(i).GetXaxis()->SetTitleSize(0.06); // changed from 0.05
      effs.at(i).GetXaxis()->SetTitleOffset(1.2);

      effs.at(i).GetYaxis()->SetTitleOffset(0.9); // changed from 1.8
      effs.at(i).GetYaxis()->SetTitleSize(0.06); // changed from 0.05
      effs.at(i).GetYaxis()->SetLabelSize(0.055); // changed from 0.045
      effs.at(i).GetYaxis()->SetTickLength(0.02);
      effs.at(i).GetYaxis()->SetLabelOffset(0.011);

      effs.at(i).GetXaxis()->SetNdivisions(505);
      effs.at(i).GetYaxis()->SetRangeUser(0,1.5);
      //effs.at(i).GetXaxis()->SetRangeUser(0,2000);

      effs.at(i).Draw("ap");
    }
    else {
      effs.at(i).Draw("p same");
    }
    leg->AddEntry(&effs.at(i), titles.at(i), "l");

    std::cout << "Done." << std::endl;
    i++;
  }

  leg->SetBorderSize(0);
  leg->Draw("same");

  // draw Lumi text
  TString infotext = TString::Format("Run 2018D (%d TeV)", 13);
  TLatex *text = new TLatex(3.5, 24, infotext);
  text->SetNDC();
  text->SetTextAlign(33);
  text->SetX(0.91);
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

  TString outputname = "plots/B2G_LeptonHT_Eff_"+histname+".pdf";
  canvas->SaveAs(outputname);

}
