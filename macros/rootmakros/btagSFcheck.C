
void btagSFcheck() {

  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.12);

  TString year = "UL18";
  TString channel = ""; // pre-pend a _ for mu and ele!
  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/"+year+"/hadded/";
  TString filename_base = "uhh2.AnalysisModuleRunner.MC.";

  //std::vector<TString> samples = {"TTbar", "WJets", "ST", "QCD", "VV", "DYJets", "TstarTstar", "TstarTstar_Spin32"};
  //std::vector<TString> samples = {"TstarTstar", "TstarTstar_Spin32"};
  std::vector<TString> samples = {"TTbar"};

  TString histname2D = "pt_HT_N_jet_rebinned";
  std::vector<TString> hists_to_crosscheck = {"pt_HT", "N_jets", "DeepJetscore"};
  std::vector<TString> labels_to_replace = {"H_{T} [GeV]", "N(small-radius jets)", "DeepJet score"};

  if(channel != "") {
    hists_to_crosscheck.push_back("pt"+channel);
    labels_to_replace.push_back("pt"+channel);
  }

  TString folder_before = "HOTVRcut"+channel;
  TString folder_after = "bcorrections"+channel;
  TString folder_crosscheck = "byield"+channel;

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
    TCanvas *can_2D = new TCanvas("can_2D", "c", 600, 600);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    hist_ratio->SetTitle("");
    hist_ratio->SetName(sample);
    hist_ratio->GetYaxis()->SetRangeUser(4,20);
    hist_ratio->GetXaxis()->SetRangeUser(0,4000);

    // label replace
    if(hist_ratio->GetXaxis()->GetTitle() == TString("H_{T}")) hist_ratio->GetXaxis()->SetTitle("H_{T} [GeV]");
    if(hist_ratio->GetYaxis()->GetTitle() == TString(" N(jet)")) hist_ratio->GetYaxis()->SetTitle("N(small-radius jets)");

    // styling
    hist_ratio->GetXaxis()->SetTitleSize(0.05);
    hist_ratio->GetYaxis()->SetTitleSize(0.05);

    hist_ratio->Draw("colz");

    // draw sample text
    TString printsample = "todo";
    if (sample == "TTbar") printsample = "t#bar{t}";
    TString printchannel = "";
    if (channel == "_ele") printchannel = "e";
    if (channel == "_mu") printchannel = "#mu";
    TLatex *text = new TLatex(3.5, 24, year + " " + printsample + " " + printchannel);
    text->SetNDC();
    text->SetTextAlign(13);
    text->SetX(0.765);
    text->SetTextFont(42);
    text->SetTextSize(0.045);
    text->SetY(0.956);
    text->Draw();

    /**
    TString cmstext = "CMS";
    TLatex *text2 = new TLatex(3.5, 24, cmstext);
    text2->SetNDC();
    text2->SetTextAlign(13);
    text2->SetX(0.15);
    text2->SetTextFont(62);
    text2->SetTextSize(0.065);
    text2->SetY(0.96);
    text2->Draw();
    **/

    TString preltext = "Private Work (CMS Simulation)";
    TLatex *text3 = new TLatex(3.5, 24, preltext);
    text3->SetNDC();
    text3->SetTextAlign(13);
    text3->SetX(0.15);
    text3->SetTextFont(52);
    text3->SetTextSize(0.045);
    text3->SetY(0.946);
    text3->Draw();

    can_2D->SaveAs("plots/btagyield/btagYieldSFs_"+channel+"_"+sample+"_"+year+".pdf");

    // storing to vector for saving later
    histograms_to_store.push_back(hist_ratio);

    // finally, we are doing some cross-checks
    int index = 0;
    for (auto histname_to_crosscheck : hists_to_crosscheck) {
      TH1D *hist_crosscheck_before = (TH1D*)input->Get(folder_before+"/"+histname_to_crosscheck);
      TH1D *hist_crosscheck_after = (TH1D*)input->Get(folder_after+"/"+histname_to_crosscheck);
      TH1D *hist_crosscheck_crosscheck = (TH1D*)input->Get(folder_crosscheck+"/"+histname_to_crosscheck);

      TCanvas *can_crosscheck = new TCanvas("can_crosscheck", "c", 600, 600);
      TPad *pad1 = new TPad("pad1", "The pad 80% of the height", 0.0, 0.35, 1.0, 1.0);
      TPad *pad2 = new TPad("pad2", "The pad 20% of the height", 0.0, 0.0, 1.0, 0.35);
      pad1->Draw();
      pad2->Draw();
      pad1->cd();
      pad1->SetTickx();
      pad1->SetTicky();
      pad1->SetLogy();

      pad1->SetBottomMargin(0);
      pad2->SetTopMargin(0);

      pad1->SetTopMargin(0.125);
      pad2->SetBottomMargin(0.3);

      hist_crosscheck_before->SetLineColor(1);
      hist_crosscheck_after->SetLineColor(2);
      hist_crosscheck_crosscheck->SetLineColor(4);

      hist_crosscheck_before->GetXaxis()->SetTitle( hist_crosscheck_before->GetTitle() );
      hist_crosscheck_before->GetYaxis()->SetTitle("Events");
      hist_crosscheck_before->GetYaxis()->SetTitleSize(0.08);
      hist_crosscheck_before->GetYaxis()->SetLabelSize(0.06);
      hist_crosscheck_before->GetYaxis()->SetTitleOffset(0.9);
      hist_crosscheck_before->GetXaxis()->SetTitleSize(0);
      hist_crosscheck_before->GetXaxis()->SetLabelSize(0);
      hist_crosscheck_before->SetTitle("");

      hist_crosscheck_before->Draw("hist");
      hist_crosscheck_after->Draw("hist same");
      hist_crosscheck_crosscheck->Draw("hist same");

      auto legend = new TLegend(0.53,0.65,0.8,0.835);
      gStyle->SetLegendTextSize(0.05);
      legend->AddEntry(hist_crosscheck_before,"before corrections","l");
      legend->AddEntry(hist_crosscheck_after,"after b-tagging SF","l");
      legend->AddEntry(hist_crosscheck_crosscheck,"after yield correction","l");

      legend->SetBorderSize(0);
      legend->Draw();

      // draw sample text
      TString printsample = "todo";
      if (sample == "TTbar") printsample = "t#bar{t}";
      TLatex *text = new TLatex(3.5, 24, year + " " + printsample);
      text->SetNDC();
      text->SetTextAlign(33);
      text->SetX(0.88);
      text->SetTextFont(42);
      text->SetY(0.95);
      text->SetTextSize(0.055);
      text->Draw();

      /**
      TString cmstext = "CMS";
      TLatex *text2 = new TLatex(3.5, 24, cmstext);
      text2->SetNDC();
      text2->SetTextAlign(13);
      text2->SetX(0.15);
      text2->SetTextFont(62);
      text2->SetTextSize(0.08);
      text2->SetY(0.955);
      text2->Draw();
      **/

      TString preltext = "Private Work (CMS simulation)";
      TLatex *text3 = new TLatex(3.5, 24, preltext);
      text3->SetNDC();
      text3->SetTextAlign(13);
      text3->SetX(0.15);
      text3->SetTextFont(52);
      text3->SetTextSize(0.06);
      text3->SetY(0.94);
      text3->Draw();

      pad2->cd();
      pad2->SetTickx();
      pad2->SetTicky();

      TH1D* hist_crosscheck_after_ratio = (TH1D*) hist_crosscheck_after->Clone();
      TH1D* hist_crosscheck_crosscheck_ratio = (TH1D*) hist_crosscheck_crosscheck->Clone();

      hist_crosscheck_after_ratio->Divide(hist_crosscheck_before);
      hist_crosscheck_crosscheck_ratio->Divide(hist_crosscheck_before);

      hist_crosscheck_after_ratio->GetXaxis()->SetTitle(hist_crosscheck_after_ratio->GetTitle());
      hist_crosscheck_after_ratio->GetXaxis()->SetTitleSize(0.13);
      hist_crosscheck_after_ratio->GetXaxis()->SetLabelSize(0.1);
      hist_crosscheck_after_ratio->GetXaxis()->SetLabelOffset(0.02);
      hist_crosscheck_after_ratio->GetXaxis()->SetTitleOffset(1.1);

      hist_crosscheck_after_ratio->GetYaxis()->SetTitle("After / before");
      hist_crosscheck_after_ratio->GetYaxis()->SetTitleSize(0.15);
      hist_crosscheck_after_ratio->GetYaxis()->SetLabelSize(0.1125);
      hist_crosscheck_after_ratio->GetYaxis()->SetNdivisions(505, kTRUE);
      hist_crosscheck_after_ratio->GetYaxis()->SetTitleOffset(0.48);

      // replace x axis titles with something more nice
      hist_crosscheck_after_ratio->GetXaxis()->SetTitle(labels_to_replace.at(index));

      hist_crosscheck_after_ratio->SetTitle("");
      hist_crosscheck_crosscheck_ratio->SetTitle("");

      if (histname_to_crosscheck == "DeepJetscore") {
        hist_crosscheck_after_ratio->GetYaxis()->SetRangeUser(0.75, 1.5);
      }

      hist_crosscheck_after_ratio->Draw();
      hist_crosscheck_crosscheck_ratio->Draw("same");

      auto line = TLine(0,1,hist_crosscheck_after_ratio->GetXaxis()->GetXmax(),1);
      line.SetLineStyle(2);
      line.Draw("same");

      can_crosscheck->SaveAs("plots/btagyield/btagYieldSFs_crosscheck_"+sample+"_"+year+"_"+channel+"_"+histname_to_crosscheck+".pdf");

      index++;

    }

  }

  if(writeSFsToFile) {
    // outputting the 2D histograms for scaling
    TFile *output = TFile::Open("files/btagyield/btagYieldSFs_"+year+"_"+channel+".root", "RECREATE");
    for (auto histogram : histograms_to_store) {
      histogram->Write();
    }
  }


}
