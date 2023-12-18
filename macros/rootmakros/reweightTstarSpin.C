

void reweightTstarSpin() {

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

  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.17);

  gROOT->ForceStyle();

  std::vector<TString> masspoints = {"700", "800", "900", "1000", "1200", "1300", "1400", "1500", "1600", "1700", "1800", "1900", "2000"};

  TString histname = "AfterCommon_gen/Pt_tstar_gen";

  std::vector<TH1D*> ratios;

  TCanvas *canvas = new TCanvas("canvas", "c", 600, 600);


  for (auto masspoint : masspoints) {
    TString path_to_UL = "/nfs/dust/cms/user/flabe/TstarTstar/data/Preselection/UL16postVFP/hadded/uhh2.AnalysisModuleRunner.MC.TstarTstar_M-"+masspoint+".root";
    TString path_to_EOY = "/nfs/dust/cms/user/flabe/TstarTstar/data/Preselection/2016/hadded/uhh2.AnalysisModuleRunner.MC.TstarTstar_M-"+masspoint+".root";

    // getting the histograms
    TFile *input_UL = TFile::Open(path_to_UL);
    TFile *input_EOY = TFile::Open(path_to_EOY);
    TH1D *hist_UL = (TH1D*)input_UL->Get(histname);
    TH1D *hist_EOY = (TH1D*)input_EOY->Get(histname);

    // normalizing
    hist_UL->Scale(1/hist_UL->Integral());
    hist_EOY->Scale(1/hist_EOY->Integral());

    // rebinning
    const int nbins = 31;
    double bins[nbins] = {0, 60, 120, 180, 240, 300, 360, 420, 480, 540, 600, 660, 720, 780, 840, 900, 960, 1020, 1080, 1140, 1200, 1260, 1320, 1380, 1440, 1500, 1620, 1740, 1860, 2100, 3000};

    hist_UL = (TH1D*)hist_UL->Rebin(nbins-1, "rebinned UL", bins);
    hist_EOY = (TH1D*)hist_EOY->Rebin(nbins-1, "rebinned EOY", bins);

    // getting the ratio
    TH1D *hist_ratio = (TH1D*)hist_EOY->Clone();
    hist_ratio->Divide(hist_UL);
    hist_ratio->SetName("ratio_"+masspoint);
    hist_ratio->GetXaxis()->SetTitle("p^{T*}_{T} [GeV]");
    hist_ratio->SetTitle("");

    // plotting the stuff
    canvas->SetLogy();

    hist_UL->SetTitle("");
    hist_UL->GetXaxis()->SetTitle("p_{T} T*");
    hist_UL->GetYaxis()->SetTitle("events [a.u.]");
    hist_UL->SetMarkerStyle(20);
    hist_UL->SetMarkerColor(1);
    hist_UL->SetLineColor(1);
    hist_UL->Draw("hist");
    hist_EOY->SetMarkerStyle(20);
    hist_EOY->SetMarkerColor(2);
    hist_EOY->SetLineColor(2);
    hist_EOY->Draw("hist same");

    auto legend = new TLegend(0.55,0.7,0.78,0.88);
    gStyle->SetLegendTextSize(0.05);
    legend->AddEntry(hist_UL,"spin 1/2","l");
    legend->AddEntry(hist_EOY,"spin 3/2","l");
    legend->Draw();

    canvas->SaveAs("plots/spincomparison_tstarpt_"+masspoint+".pdf");

    ratios.push_back(hist_ratio);

  }

  //canvas->SetLogy(false);

  TFile *output = TFile::Open("files/spinFactors.root", "RECREATE");

  // plot ratio
  int i = 1;
  for (auto hist_ratio : ratios) {

    // prevent negatives (no idea where those are coming from!)
    for(int j=1; j < hist_ratio->GetXaxis()->GetNbins()+1; j++){
      auto content = hist_ratio->GetBinContent(j);
      if(content < 0) hist_ratio->SetBinContent(j, hist_ratio->GetBinContent(j-1));
    }

    hist_ratio->GetYaxis()->SetTitle("ratio");
    hist_ratio->SetLineColor(i);
    if(i == 1) hist_ratio->Draw("hist");
    else hist_ratio->Draw("hist same");
    i++;

    hist_ratio->Write();

  }

  canvas->SaveAs("plots/spincomparison_tstarpt_ratio.pdf");

}
