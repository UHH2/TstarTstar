#pragma once

// this is a test 2

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>

#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/common/include/MCWeight.h"

enum ParticleID { H=25, W=24, Z=23, t=6, b=5, c=4, s=3, d=2, u=1, g=21};
inline bool FindInString(const std::string& search, const std::string& str) {return str.find(search)!=std::string::npos ;}

// Generic Class for Applying SFs
class ScaleFactorsFromHistos : public uhh2::AnalysisModule {

public:
  void LoadHisto(TFile* file, std::string name, std::string hname);
  double Evaluator(std::string hname, double var);

protected:
  std::unordered_map<std::string, std::unique_ptr<TH1F> > histos;

};

// Apply Theory weights //TODO make multiple inheritance
class NLOCorrections : public ScaleFactorsFromHistos {

public:
  explicit NLOCorrections(uhh2::Context& ctx);
  virtual bool process(uhh2::Event&) override;
  double GetPartonObjectPt(uhh2::Event& event, ParticleID objID);

private:
  bool is_Wjets, is_Zjets, is_DY, is_Znn, is2016;

  // Theory Corrections
  const bool do_EWK = true;
  const bool do_QCD_EWK = false;
  const bool do_QCD_NLO  = true;
  const bool do_QCD_NNLO = false;

};

// HEM issue addressation for MC
class EtaPhiEventCleanerMC: public uhh2::AnalysisModule {
public:
  EtaPhiEventCleanerMC(uhh2::Context& ctx, float weight_factor, float min_eta, float max_eta, float min_phi, float max_phi, std::string jetCollection = "topjets", bool doJets = true, bool doElectrons = true, bool doMuons = true);
  virtual bool process(uhh2::Event& event) override;

private:
  float weight_factor;
  float min_eta;
  float max_eta;
  float min_phi;
  float max_phi;
  std::string jetCollection;
  bool doJets;
  bool doElectrons;
  bool doMuons;
  uhh2::Event::Handle<std::vector<Jet> > h_jets;
  uhh2::Event::Handle<std::vector<TopJet> > h_topjets;
};

// Class which inherits from EtaPhiEventCleanerMC, specific to the HEM15/16 issues.
// This applies to 2018 only.
// The values are specified in https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
class HEMCleanerMCScale: public EtaPhiEventCleanerMC {
public:
  static constexpr float min_eta_HEM = -2.964;
  static constexpr float max_eta_HEM = -1.305;
  static constexpr float min_phi_HEM = -1.6;
  static constexpr float max_phi_HEM = -0.87;
  static constexpr float weight_factor = 0.352;

  HEMCleanerMCScale(uhh2::Context& ctx, std::string jetCollection, bool doJets=true, bool doElectrons=true, bool doMuons=true);

private:
  std::string jetCollection;
};

class TstarTstarSpinScale: public uhh2::AnalysisModule {
public:

  explicit TstarTstarSpinScale(uhh2::Context& ctx, TString path);
  virtual bool process(uhh2::Event&) override;

private:
  bool doSpinReweighting = false;
  TH1D *scaleFactorHist;

};
