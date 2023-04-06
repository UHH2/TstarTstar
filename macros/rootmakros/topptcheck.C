

void topptcheck(){

  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/UL18/hadded";
  std::vector<TString> samples = {"TTbar", "WJets", "VV", "QCD", "ST"};

  TString beforeFolder = "AfterCorrections";
  TString afterFolder = "Aftertopptreweighting";

  TString histname = "N_jets";

  // ##########



  TH1D *hist_before;
  TH1D *hist_after;

  bool first = true;
  for (auto sample : samples) {
    TString filename = "uhh2.AnalysisModuleRunner.MC."+sample+".root";
    TFile *file = TFile::Open(path+"/"+filename);
    if(!file) std::cout << "File does not exist!" << std::endl;

    if(first) {
      hist_before = (TH1D*)file->Get(beforeFolder+"/"+histname); //histogram
      hist_after = (TH1D*)file->Get(afterFolder+"/"+histname); //histogram
    } else {
      hist_before->Add( (TH1D*)file->Get(beforeFolder+"/"+histname) );
      hist_after->Add( (TH1D*)file->Get(afterFolder+"/"+histname) );
    }

    first = false;
  }

  double before = hist_before->Integral();
  double after = hist_after->Integral();

  // ##########

  std::cout << "Integral before: " << before << std::endl;
  std::cout << "Integral after: " << after << std::endl;
  std::cout << "Thus, the scaling factor should be: " << before/after << std::endl;


}
