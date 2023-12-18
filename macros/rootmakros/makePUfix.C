
void makePUfix() {

  TString year = "UL18";
  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/"+year+"/hadded/";
  TString filename_base = "uhh2.AnalysisModuleRunner.MC.";

  std::vector<TString> samples = {"TTbar", "WJets", "ST", "QCD", "VV"};

  TString filename_data = "uhh2.AnalysisModuleRunner.DATA.DATA.root";

  TString folder = "AfterTrigger";
  TString histname = "N_pv";

  bool writeSFsToFile = true;

  // main loop, we are doing this for each sample
  std::vector<TH2D*> histograms_to_store;
  for (const auto sample : samples) {
    std::cout << "Processing " << sample << "..." << std::endl;

    TFile *input_MC = TFile::Open(path+filename_base+sample+".root");
    if(!input_MC) cout << "Empty file MC" << endl;
    TH1D *hist_MC = (TH1D*)input_MC->Get(folder+"/"+histname);
    if(!hist_MC) cout << "Empty hist MC" << endl;

    TFile *input_data = TFile::Open(path+filename_data);
    if(!input_data) cout << "Empty file data" << endl;
    TH1D *hist_data = (TH1D*)input_data->Get(folder+"/"+histname);
    if(!hist_data) cout << "Empty hist data" << endl;

    hist_MC->Scale(1/hist_MC->Integral());
    hist_data->Scale(1/hist_data->Integral());

    // clone it, divide, then close the file
    TH2D *hist_ratio = (TH2D*)hist_data->Clone();
    hist_ratio->Divide(hist_MC);

    hist_ratio->SetName(sample);

    histograms_to_store.push_back(hist_ratio);

  }


  if(writeSFsToFile) {
    // outputting the 2D histograms for scaling
    TFile *output = TFile::Open("PUfakeSF_"+year+".root", "RECREATE");
    for (auto histogram : histograms_to_store) {
      histogram->Write();
    }

  }

}
