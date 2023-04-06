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
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"


class TstarTstarGenDiscriminator : uhh2::AnalysisModule{

public:
  explicit TstarTstarGenDiscriminator(uhh2::Context&, bool do_only_ttbar=false);
  virtual bool process(uhh2::Event&) override;

 private:
  uhh2::Event::Handle<std::vector<ReconstructionTstarHypothesis>> h_tstartstar_hyp_vector;
  uhh2::Event::Handle<ReconstructionTstarHypothesis> h_tstartstar_hyp;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;

  bool is_tgtg, is_tgtgamma, only_ttbar;

};

class TstarTstarGenMatcher : uhh2::AnalysisModule{

public:
  explicit TstarTstarGenMatcher(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

 private:
  std::unique_ptr<uhh2::TTbarSemiLepMatchableSelection> TTbarSemiLepMatchable_selection;

  uhh2::Event::Handle<std::vector<ReconstructionTstarHypothesis>> h_tstartstar_hyp_vector;
  uhh2::Event::Handle<ReconstructionTstarHypothesis> h_tstartstar_hyp;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;

  bool is_tgtg, is_tgtgamma;

};
