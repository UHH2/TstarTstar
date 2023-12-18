
void combineErrorOnce(TString file, TString histpath, TString errorbandfile, TString errorbandpath, TString outfilename) {

  // get histogram
  TFile *hist_file = TFile::Open(file);
  if(!hist_file) std::cout << "File " << file << " does not exist" << std::endl;
  TH1D *histogram = (TH1D*)hist_file->Get(histpath);
  if(!histogram) std::cout << "Hist at " << histpath << " does not exist" << std::endl;

  // load error band
  TFile *errorband_file = TFile::Open(errorbandfile);
  if(!errorband_file) std::cout << "File " << errorbandfile << "does not exist" << std::endl;
  TH1D *errorband = (TH1D*)errorband_file->Get(errorbandpath);
  if(!errorband) std::cout << "Hist at " << errorbandpath << " does not exist" << std::endl;

  // now lets loop over all bins of the histogram and add the error
  for (int i = 0; i <= histogram->GetXaxis()->GetNbins(); i++) {

    // get existing bin error
    auto staterror = histogram->GetBinError(i);

    // get x position value at bin
    auto x_pos = histogram->GetXaxis()->GetBinCenter(i);

    // get fit value at bin
    auto hist_val = histogram->GetBinContent(i);
    if(hist_val == 0) continue;

    // get fit value at bin
    auto fit_val = errorband->GetBinContent(errorband->GetXaxis()->FindBin(x_pos));

    // get error from error band
    auto banderror = errorband->GetBinError(errorband->GetXaxis()->FindBin(x_pos));

    // do error propagation of multiplication of two values
    auto relativeErrorA = staterror/hist_val;
    auto relativeErrorB = banderror/fit_val;

    auto totalrelativeerror = TMath::Sqrt( relativeErrorA*relativeErrorA + relativeErrorB*relativeErrorB );

    // setting bin error
    histogram->SetBinError(i, totalrelativeerror * hist_val);
  }

  // save resulting histogram
  TFile *output = TFile::Open(outfilename, "RECREATE");
  histogram->Write();

}

void backgroundEstimation_combineStatError(){

  // do the following twice for both functions
  // 1) load the histogram from the output file
  // 2) load the corresponding fit functions error band
  // 3) add error from band to stat error from hist in quadrature
  // 4) store result somewhere for combine to use

  TString year = "UL16postVFP";
  TString region = "newTaggerCR";

  // first, the nominal case
  TString nominalFile = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN_datadriven/"+year+"/hadded/uhh2.AnalysisModuleRunner.DATA.datadrivenBG.root";
  TString histpath = region+"/pt_ST_rebinned";
  TString errorbandFile = "/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/alphaFunction.root";
  TString errorbandpathnominal = "hint2";
  TString outfilenominal = "files/bgest_"+year+"_nominal.root";

  combineErrorOnce(nominalFile, histpath, errorbandFile, errorbandpathnominal, outfilenominal);

  TString variationFile = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN_datadriven_variation/"+year+"/hadded/uhh2.AnalysisModuleRunner.DATA.datadrivenBG.root";
  TString errorbandpathvariation = "hint";
  TString outfilevariation = "files/bgest_"+year+"_variation.root";
  combineErrorOnce(variationFile, histpath, errorbandFile, errorbandpathvariation, outfilevariation);

}
