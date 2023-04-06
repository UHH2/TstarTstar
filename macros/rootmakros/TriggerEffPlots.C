// Efficiency histograms for trigger related studies
// author: A.Karavdina, changes F.Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q TriggerEffPlots.C

void TriggerEffPlots(){

  std::vector<TString> filenames = {"uhh2.AnalysisModuleRunner.MC.TTbar.root", "uhh2.AnalysisModuleRunner.DATA.DATA.root"};
  std::vector<TString> labels = {"TTbar", "DATA"};

  Double_t w = 400;
  Double_t h = 400;

  TCanvas *m_can = new TCanvas("canvas", "canvas", w, h);
  TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
  leg->SetBorderSize(0);

  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.16);  gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.05);

  // ###################### end cosmetics #######################

  TString channel = "mu";
  TString vsVar = "pt";
  TString subpath_pre="AfterCorrections";
  TString subpath_post="AfterTriggerSF";
  TString suffix = "afterSF";

  //Files after selection
  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/UL18/hadded/";

  std::vector<TGraphAsymmErrors> graphs = {TGraphAsymmErrors(), TGraphAsymmErrors()};

  for (int ifile = 0; ifile < filenames.size(); ifile++) {
    TString filename = filenames.at(ifile);
    TString label = labels.at(ifile);

    TFile *input = TFile::Open(path+filename);
    TH1D *hist_trigger = (TH1D*)input->Get(subpath_post+"_"+channel+"/"+vsVar+"_"+channel);
    TH1D *hist_denom = (TH1D*)input->Get(subpath_pre+"_"+channel+"/"+vsVar+"_"+channel);

    hist_trigger->Rebin(1);
    hist_denom->Rebin(1);

    // checking if everything is proper
    for (int i = 0; i < hist_trigger->GetXaxis()->GetNbins()+2; i++) {
      if(hist_trigger->GetBinContent(i) > hist_denom->GetBinContent(i)) {
        std::cout << "Bin " << i << " Before: " << hist_denom->GetBinContent(i) << " - After: " << hist_trigger->GetBinContent(i) << std::endl;
        hist_trigger->SetBinContent(i, hist_denom->GetBinContent(i));
      }
    }

    if(!hist_trigger || !hist_denom) cout<<"Hists are empty"<<endl;;
    if(!hist_trigger || !hist_denom) return;

    //TEfficiency eff;
    if(hist_denom->GetEntries()>0 && hist_trigger->GetEntries()>0) graphs.at(ifile).Divide(hist_trigger, hist_denom, "cl=0.68 b(1,1) mode");

  }

  auto eff_MC = graphs.at(0);
  auto eff_DATA = graphs.at(1);

  // cosmetics
  eff_MC.GetXaxis()->SetTitle(channel + " " + vsVar);
  eff_MC.GetYaxis()->SetTitle("Efficiency");
  eff_MC.SetTitle("");

  eff_MC.GetXaxis()->SetTitleFont(42);
  eff_MC.GetXaxis()->SetLabelFont(42);
  eff_MC.GetYaxis()->SetTitleFont(42);
  eff_MC.GetYaxis()->SetLabelFont(42);

  eff_MC.GetXaxis()->SetLabelSize(0.055); // changed from 0.045
  eff_MC.GetXaxis()->SetLabelOffset(0.008);
  eff_MC.GetXaxis()->SetTickLength(0.03);
  eff_MC.GetXaxis()->SetTitleSize(0.06); // changed from 0.05
  eff_MC.GetXaxis()->SetTitleOffset(1.2);

  eff_MC.GetYaxis()->SetTitleOffset(1.4); // changed from 1.8
  eff_MC.GetYaxis()->SetTitleSize(0.06); // changed from 0.05
  eff_MC.GetYaxis()->SetLabelSize(0.055); // changed from 0.045
  eff_MC.GetYaxis()->SetTickLength(0.02);
  eff_MC.GetYaxis()->SetLabelOffset(0.011);

  eff_MC.GetYaxis()->SetRangeUser(0,1.5);
  //eff_MC.GetXaxis()->SetRangeUser(0,1000);

  eff_MC.SetMarkerStyle(1);
  eff_MC.SetLineWidth(2);
  eff_MC.SetMarkerColor(1);
  eff_MC.Draw("AP");
  leg->AddEntry(&eff_MC, labels.at(0), "l");

  eff_DATA.SetMarkerStyle(1);
  eff_DATA.SetLineWidth(2);
  eff_DATA.SetMarkerColor(2);
  eff_DATA.SetLineColor(2);
  eff_DATA.Draw("P SAME");
  leg->AddEntry(&eff_DATA, labels.at(1), "l");

  TString infotext = TString::Format("%3.1f fb^{-1} (%d TeV)", 137., 13);
  TLatex *text = new TLatex(3.5, 24, infotext);
  text->SetNDC();
  text->SetTextAlign(33);
  text->SetX(0.95);
  text->SetTextFont(42);
  text->SetY(1);
  text->SetTextSize(0.045);
  //text->Draw();

  leg->Draw();

  if(suffix == "") m_can->SaveAs("plots/TrgEff_"+vsVar+"_"+channel+".pdf");
  else m_can->SaveAs("plots/TrgEff_"+vsVar+"_"+channel+"_"+suffix+".pdf");
}
