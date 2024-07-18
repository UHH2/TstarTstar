#pragma once

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/NtupleObjects.h>
#include <UHH2/core/include/LorentzVector.h>
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/Utils.h"

class PuppiCHSMatching : public uhh2::AnalysisModule {

public:
  explicit PuppiCHSMatching(uhh2::Context& ctx, std::string chs_jet_handle);
  virtual bool process(uhh2::Event&) override;

private:
  uhh2::Event::Handle< std::vector<Jet> > h_CHSjets;
  uhh2::Event::Handle< std::vector<Jet> > h_CHS_matched_;
};
