#pragma once
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/AnalysisModule.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <TF1.h>
#include <TGraph.h>
#include <TFile.h>

using namespace std;
using namespace uhh2;

Particle GetLepton(uhh2::Event & event);

// stolen form Alex, not sure yet what trhis does!

class CorrectionFactor: public uhh2::AnalysisModule{
public:

  explicit CorrectionFactor(uhh2::Context &,  const std::string &, const std::string &, bool allHad_, const TString &year_);
  virtual bool process(uhh2::Event & ) override;

private:
  void get_function(const TString &);
  void get_additionalSYS();
  double get_factor(double, double);
  uhh2::Event::Handle<std::vector<TopJet>>h_oldjets;
  uhh2::Event::Handle<std::vector<TopJet>>h_newjets;
  uhh2::Event::Handle<vector<float>>h_cor_factor_had;
  uhh2::Event::Handle<vector<float>>h_cor_factor_lep;
  //double arr[6][12]; // values from table
  //double par[12][3]; // values with function parameters in 12 eta bins
  double eta_bins[13] = {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930, 2.322, 2.65, 3.139, 5.191};
  std::vector<TF1*> CentralCorrectionFunctions;
  std::vector<TGraph*> UpDownCorrectionGraphs;
  TGraph *AdditionalSys;
  bool CorUp, CorDown;
  bool allHad;
  string year;

};
