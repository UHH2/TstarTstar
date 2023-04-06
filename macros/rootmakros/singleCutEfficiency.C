

void singleCutEfficiency(){

  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/Analysis/hadded/";
  TString fileprefix = "uhh2.AnalysisModuleRunner.MC.";
  TString sample = "TTbar";

  TString step_1 = "main";
  TString step_2 = "reco";
  TString histname_1 = "pt_ST";
  TString histname_2 = "M_Tstar_gHOTVR";

  TFile *input;
  TH1D *hist1, *hist2;

  input = TFile::Open(path+fileprefix+sample+".root");
  hist1 = (TH1D*)input->Get(step_1+"/"+histname_1);
  hist2 = (TH1D*)input->Get(step_2+"/"+histname_2);

  double int1 = hist1->Integral();
  double int2 = hist2->Integral();

  std::cout << "Sample: " << sample << std::endl;
  std::cout << "main: " << int1 << std::endl;
  std::cout << "reco: " << int2 << std::endl;
  std::cout << "Ratio: " << int2/int1 << std::endl;
}
