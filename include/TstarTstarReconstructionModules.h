#pragma once

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/NtupleObjects.h>
#include <UHH2/core/include/LorentzVector.h>
#include <UHH2/common/include/TTbarGen.h>

#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/ReconstructionHypothesis.h"
#include "UHH2/common/include/ReconstructionHypothesisDiscriminators.h"
#include "UHH2/TstarTstar/include/ReconstructionTstarHypothesis.h"

class TopTagMassWindow {
public:

  explicit TopTagMassWindow(double mlower=105., double mupper=220.);

  bool operator()(const TopJet & topjet, const uhh2::Event & event) const;

 private:
  double m_mlower;
  double m_mupper;
};

class TopTagSubbtag {
public:

  explicit TopTagSubbtag(double btag=0.5);

  bool operator()(const TopJet & topjet, const uhh2::Event & event) const;

 private:
  double m_btag;
};


class ttbarChi2Discriminator : uhh2::AnalysisModule{

public:
  explicit ttbarChi2Discriminator(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

private:
  uhh2::Event::Handle< std::vector<ReconstructionHypothesis> > h_ttbar_hyps_;
  uhh2::Event::Handle<ReconstructionHypothesis> h_recohyp_;
  uhh2::Event::Handle<bool> h_is_ttbar_reconstructed_;
  uhh2::Event::Handle<int> h_ttag_jet_pos;
  float mtoplep_, mtoplep_ttag_;
  float sigmatoplep_, sigmatoplep_ttag_;
  float mtophad_, mtophad_ttag_;
  float sigmatophad_, sigmatophad_ttag_;

};


class TstarTstar_tgtg_TopTag_Reconstruction : uhh2::AnalysisModule{

 public:
  explicit TstarTstar_tgtg_TopTag_Reconstruction(uhh2::Context&, const NeutrinoReconstructionMethod&, TopJetId id=CMSTopTag());
  virtual bool process(uhh2::Event&) override;

 private:
  NeutrinoReconstructionMethod m_neutrinofunction;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
  uhh2::Event::Handle<std::vector<ReconstructionTstarHypothesis>> h_tstartstar_hyp_vector;
  uhh2::Event::Handle<int> h_flag_toptagevent;

  TopJetId topjetID_;
};

class TstarTstar_tgtg_AK4_Reconstruction : uhh2::AnalysisModule{

 public:
  explicit TstarTstar_tgtg_AK4_Reconstruction(uhh2::Context&, const NeutrinoReconstructionMethod&, TopJetId id=CMSTopTag());
  virtual bool process(uhh2::Event&) override;

 private:
  NeutrinoReconstructionMethod m_neutrinofunction;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
  uhh2::Event::Handle<std::vector<ReconstructionTstarHypothesis>> h_tstartstar_hyp_vector;
  uhh2::Event::Handle<int> h_flag_toptagevent;

  TopJetId topjetID_;
};


class TstarTstar_Discrimination : uhh2::AnalysisModule{

 public:
  explicit TstarTstar_Discrimination(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

 private:
  uhh2::Event::Handle<std::vector<ReconstructionTstarHypothesis>> h_tstartstar_hyp_vector;
  uhh2::Event::Handle<ReconstructionTstarHypothesis> h_tstartstar_hyp;
  bool is_HOTVR;   // For toggling softdropmass or "normal" mass depending on Jet collection
};
