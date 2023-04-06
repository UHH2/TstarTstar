
void btagSFcheck() {

  TString year = "UL16preVFP";
  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/"+year+"/hadded/";
  TString filename_base = "uhh2.AnalysisModuleRunner.MC.";

  std::vector<TString> samples = {"TTbar", "WJets", "ST", "QCD", "VV", "DYJets", "TstarTstar"};

  TString histname2D = "pt_HT_N_jet_rebinned";
  std::vector<TString> hists_to_crosscheck = {"pt_HT", "N_jets", "DeepJetscore"};

  TString folder_before = "BeforeBCorrections";
  TString folder_after = "AfterBCorrections";
  TString folder_crosscheck = "AfterBYieldCorrections";

  bool writeSFsToFile = false;

  // main loop, we are doing this for each sample
  std::vector<TH2D*> histograms_to_store;
  for (const auto sample : samples) {
    std::cout << "Processing " << sample << "..." << std::endl;

    TFile *input = TFile::Open(path+filename_base+sample+".root");
    if(!input) cout << "Empty file" << endl;
    TH2D *hist_before = (TH2D*)input->Get(folder_before+"/"+histname2D);
    if(!hist_before) cout << "Empty hist before" << endl;
    TH2D *hist_after = (TH2D*)input->Get(folder_after+"/"+histname2D);
    if(!hist_after) cout << "Empty hist after" << endl;

    // clone it, divide, then close the file
    TH2D *hist_ratio = (TH2D*)hist_before->Clone();
    hist_ratio->Divide(hist_after);

    // drawing the 2D histogram
    TCanvas *can_2D = new TCanvas("can_2D", "c", 500, 500);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    hist_ratio->SetTitle("");
    hist_ratio->SetName(sample);
    hist_ratio->GetYaxis()->SetRangeUser(4,20);
    hist_ratio->GetXaxis()->SetRangeUser(0,4000);
    hist_ratio->Draw("colz");

    // draw sample text
    TLatex *text = new TLatex(3.5, 24, year + " " + sample);
    text->SetNDC();
    text->SetTextAlign(33);
    text->SetX(0.88);
    text->SetTextFont(42);
    text->SetY(0.95);
    text->SetTextSize(0.045);
    text->Draw();

    can_2D->SaveAs("plots/btagYieldSFs_"+sample+"_"+year+".pdf");

    // storing to vector for saving later
    histograms_to_store.push_back(hist_ratio);

    // finally, we are doing some cross-checks
    for (auto histname_to_crosscheck : hists_to_crosscheck) {
      TH1D *hist_crosscheck_before = (TH1D*)input->Get(folder_before+"/"+histname_to_crosscheck);
      TH1D *hist_crosscheck_after = (TH1D*)input->Get(folder_after+"/"+histname_to_crosscheck);
      TH1D *hist_crosscheck_crosscheck = (TH1D*)input->Get(folder_crosscheck+"/"+histname_to_crosscheck);

      TCanvas *can_crosscheck = new TCanvas("can_crosscheck", "c", 500, 500);
      TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.35,1.0,1.0);
      TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.35);
      pad1->Draw();
      pad2->Draw();
      pad1->cd();
      pad1->SetLogy();

      pad1->SetBottomMargin(0);
      pad2->SetTopMargin(0);

      hist_crosscheck_before->SetLineColor(1);
      hist_crosscheck_after->SetLineColor(2);
      hist_crosscheck_crosscheck->SetLineColor(3);

      hist_crosscheck_before->GetXaxis()->SetTitle( hist_crosscheck_before->GetTitle() );
      hist_crosscheck_before->GetXaxis()->SetTitle("events");
      hist_crosscheck_before->SetTitle("");

      hist_crosscheck_before->Draw("hist");
      hist_crosscheck_after->Draw("hist same");
      hist_crosscheck_crosscheck->Draw("hist same");

      auto legend = new TLegend(0.55,0.7,0.78,0.88);
      gStyle->SetLegendTextSize(0.05);
      legend->AddEntry(hist_crosscheck_before,"before corrections","l");
      legend->AddEntry(hist_crosscheck_after,"after b-tagging SF","l");
      legend->AddEntry(hist_crosscheck_crosscheck,"after yield correction","l");

      legend->SetBorderSize(0);
      legend->Draw();

      // draw sample text
      TLatex *text = new TLatex(3.5, 24, year + "_" + sample);
      text->SetNDC();
      text->SetTextAlign(33);
      text->SetX(0.88);
      text->SetTextFont(42);
      text->SetY(0.95);
      text->SetTextSize(0.045);
      text->Draw();

      pad2->cd();

      TH1D* hist_crosscheck_after_ratio = (TH1D*) hist_crosscheck_after->Clone();
      TH1D* hist_crosscheck_crosscheck_ratio = (TH1D*) hist_crosscheck_crosscheck->Clone();

      hist_crosscheck_after_ratio->Divide(hist_crosscheck_before);
      hist_crosscheck_crosscheck_ratio->Divide(hist_crosscheck_before);

      hist_crosscheck_after_ratio->SetTitle("");
      hist_crosscheck_crosscheck_ratio->SetTitle("");

      hist_crosscheck_after_ratio->Draw();
      hist_crosscheck_crosscheck_ratio->Draw("same");

      auto line = TLine(0,1,hist_crosscheck_after_ratio->GetXaxis()->GetXmax(),1);
      line.SetLineStyle(2);
      line.Draw("same");

      can_crosscheck->SaveAs("plots/btagYieldSFs_crosscheck_"+sample+"_"+year+"_"+histname_to_crosscheck+".pdf");

    }

  }

  if(writeSFsToFile) {
    // outputting the 2D histograms for scaling
    TFile *output = TFile::Open("files/btagYieldSFs_"+year+".root", "RECREATE");
    for (auto histogram : histograms_to_store) {
      histogram->Write();
    }
  }


}
