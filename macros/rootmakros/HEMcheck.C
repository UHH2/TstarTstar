

void HEMcheck(){

  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/UL18/hadded";
  TString sample = "DATA.DATA";

  TString beforeFolder = "AfterTrigger";
  TString afterFolder = "AfterHEMcleaning";

  TString histname = "N_jets";

  // ##########

  TString filename = "uhh2.AnalysisModuleRunner."+sample+".root";

  TFile *file = TFile::Open(path+"/"+filename);
  if(!file) std::cout << "File does not exist!" << std::endl;

  TH1D *hist_before = (TH1D*)file->Get(beforeFolder+"/"+histname); //histogram
  TH1D *hist_after = (TH1D*)file->Get(afterFolder+"/"+histname); //histogram
  if(!hist_before) std::cout << "Before hist does not exist!" << std::endl;
  if(!hist_after) std::cout << "Denominator does not exist!" << std::endl;

  double before = hist_before->Integral();
  double after = hist_after->Integral();

  // ##########

  std::cout << "Evaluating HEM fix for " << filename << std::endl;
  std::cout << "Integral before: " << before << std::endl;
  std::cout << "Integral after: " << after << std::endl;
  std::cout << "Thus, the scaling factor should be: " << after/before << std::endl;



}
