

void reweightTstarSpin2D() {

  std::vector<TString> masspoints = {"700", "800", "900", "1000", "1200", "1300", "1400", "1500", "1600", "1700", "1800", "1900", "2000"};

  TString histname_hist = "Pt_M_tstar_2D_gen";
  TString histname = "AfterCommon_gen/"+histname_hist;

  std::vector<TH2D*> ratios;

  for (auto masspoint : masspoints) {
    TString path_to_UL = "/nfs/dust/cms/user/flabe/TstarTstar/data/Preselection/UL16postVFP/hadded/uhh2.AnalysisModuleRunner.MC.TstarTstar_M-"+masspoint+".root";
    TString path_to_EOY = "/nfs/dust/cms/user/flabe/TstarTstar/data/Preselection/2016/hadded/uhh2.AnalysisModuleRunner.MC.TstarTstar_M-"+masspoint+".root";

    // getting the histograms
    TFile *input_UL = TFile::Open(path_to_UL);
    TFile *input_EOY = TFile::Open(path_to_EOY);
    TH2D *hist_UL = (TH2D*)input_UL->Get(histname);
    TH2D *hist_EOY = (TH2D*)input_EOY->Get(histname);

    // normalizing
    hist_UL->Scale(1/hist_UL->Integral());
    hist_EOY->Scale(1/hist_EOY->Integral());

    // getting the ratio
    TH2D *hist_ratio = (TH2D*)hist_EOY->Clone();
    hist_ratio->Divide(hist_UL);
    hist_ratio->SetName("ratio_"+masspoint);

    ratios.push_back(hist_ratio);

  }

  //canvas->SetLogy(false);

  TFile *output = TFile::Open("files/spinFactors_"+histname_hist+".root", "RECREATE");

  // plot ratio
  int i = 1;
  for (auto hist_ratio : ratios) {

    // prevent negatives (no idea where those are coming from!)
    for(int j=1; j < hist_ratio->GetXaxis()->GetNbins()+1; j++){
      for(int k=1; k < hist_ratio->GetYaxis()->GetNbins()+1; k++){
        auto content = hist_ratio->GetBinContent(j, k);
        if(content < 0) hist_ratio->SetBinContent(j, k, hist_ratio->GetBinContent(j-1, k-1));
      }
    }

    hist_ratio->Write();

  }


}
