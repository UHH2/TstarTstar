// Efficiency histograms for trigger related studies
// author: F.Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q cutflowPlots.C

void totalEff(){

  // Reuseable buffers
  TFile *input;
  TH1D *hist;

  TString pathReco = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/hadded/";
  TString pathDNN = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/hadded/";
  TString fileprefix = "uhh2.AnalysisModuleRunner.MC.";
  TString histname = "N_jets";

  TString before =  "AfterST";
  TString after = "AfterDNNcut_06_UGLYFIX";

  // Defining Samples
  std::vector<TString> signalSamples = {"TstarTstar_M-700", "TstarTstar_M-800", "TstarTstar_M-900", "TstarTstar_M-1000", "TstarTstar_M-1200", "TstarTstar_M-1300", "TstarTstar_M-1400", "TstarTstar_M-1500", "TstarTstar_M-1600", "TstarTstar_M-1700", "TstarTstar_M-1800", "TstarTstar_M-1900", "TstarTstar_M-2000"};
  std::vector<TString> BGSamples = {"TTbar", "ST", "WJets", "QCD", "VV"};

  std::vector<TString> signal_labels = {"TstarTstar_M-700", "TstarTstar_M-800", "TstarTstar_M-900", "TstarTstar_M-1000", "TstarTstar_M-1200", "TstarTstar_M-1300", "TstarTstar_M-1400", "TstarTstar_M-1500", "TstarTstar_M-1600", "TstarTstar_M-1700", "TstarTstar_M-1800", "TstarTstar_M-1900", "TstarTstar_M-2000"};
  std::vector<TString> BG_labels = {"TTbar", "ST", "WJets", "QCD", "VV"};

  // ########################
  // ## Finish Definitions ##
  // ########################

  // Saving initial values
  std::vector<double> initial_signal;
  std::vector<double> initial_BG;
  double initial_BG_sum = 0;
  for(const auto & sample : signalSamples){
    input = TFile::Open(pathReco+fileprefix+sample+".root");
    if (!input) cout << "First file not found" << endl;
    hist = (TH1D*)input->Get(before+"/"+histname);
    if (!hist) cout << "First hist not found" << endl;
    initial_signal.push_back(hist->Integral());
  }
  for(const auto & sample : BGSamples){
    input = TFile::Open(pathReco+fileprefix+sample+".root");
    if (!input) cout << "Second file not found" << endl;
    hist = (TH1D*)input->Get(before+"/"+histname);
    if (!hist) cout << "second hist not found" << endl;
    initial_BG.push_back(hist->Integral());
    initial_BG_sum += hist->Integral();
  }

  // Saving result values
  std::vector<double> result_signal;
  std::vector<double> result_BG;
  double result_BG_sum = 0;
  for(const auto & sample : signalSamples){
    input = TFile::Open(pathDNN+fileprefix+sample+".root");
    if (!input) cout << "third file not found" << endl;
    hist = (TH1D*)input->Get(after+"/"+histname);
    if (!hist) cout << "third hist not found" << endl;
    result_signal.push_back(hist->Integral());
  }
  for(const auto & sample : BGSamples){
    input = TFile::Open(pathDNN+fileprefix+sample+".root");
    if (!input) cout << "FOurth file not found" << endl;
    hist = (TH1D*)input->Get(after+"/"+histname);
    if (!hist) cout << "fourth hist not found" << endl;
    result_BG.push_back(hist->Integral());
    result_BG_sum += hist->Integral();
  }


  // Output DNN efficiency
  for(uint i = 0; i < signalSamples.size(); i++){
    double val = result_signal.at(i)/initial_signal.at(i);
    std::cout << "DNN efficiency for " << signal_labels.at(i) << ": " << val << std::endl;
  }

  for(uint i = 0; i < BGSamples.size(); i++){
    double val = result_BG.at(i)/initial_BG.at(i);
    std::cout << "DNN efficiency for " << BG_labels.at(i) << ": " << val << std::endl;
  }

  {
    double val = result_BG_sum/initial_BG_sum;
    std::cout << "DNN efficiency for BG sum: " << val << std::endl;

  }

}
