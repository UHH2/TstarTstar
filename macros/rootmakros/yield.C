// Efficiency histograms for trigger related studies
// author: F.Labe
// date: 24.09.2019
// Run it with following command:
// root -l -q yield.C

void yield(){

  TString filename = "TTbar";
  TString year = "2018";

  // Reuseable buffers
  TFile *input;
  TH1D *hist;

  // Defining paths
  //TString pathPresel = "/nfs/dust/cms/user/flabe/TstarTstar/data/Preselection/"+year+"/hadded/";
  TString pathPresel = "/nfs/dust/cms/user/flabe/TstarTstar/data/Analysis/"+year+"/hadded/";
  //TString pathReco = "/nfs/dust/cms/user/flabe/TstarTstar/data/Analysis/"+year+"/hadded/";
  TString pathReco = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/"+year+"/hadded/";
  TString fileprefix = "uhh2.AnalysisModuleRunner.MC.";
  TString histname = "N_jets";

  // Defining Steps
  TString pre = "main";
  TString post = "AfterDNNcut_08";

  // Saving initial values
  input = TFile::Open(pathPresel+fileprefix+filename+".root");
  hist = (TH1D*)input->Get(pre+"/"+histname);
  double before = hist->Integral();

  if(filename == "WJets" && year == "2017"){
    std::cout << "Adding second file as this is wjets" << std::endl;
    input = TFile::Open(pathPresel+fileprefix+filename+"_1.root");
    hist = (TH1D*)input->Get(pre+"/"+histname);
    before += hist->Integral();
  }

  // Saving initial values
  input = TFile::Open(pathReco+fileprefix+filename+".root");
  hist = (TH1D*)input->Get(post+"/"+histname);
  double after = hist->Integral();

  double efficiency = after/before;

  std::cout << filename << ": " << efficiency << std::endl;

}
