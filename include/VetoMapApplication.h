#pragma once

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>

class VetoMapApplicator: public uhh2::AnalysisModule {
public:

  explicit VetoMapApplicator(uhh2::Context& ctx);
  virtual bool process(uhh2::Event&) override;

private:

    // cleaners to apply veto map IDs
    std::unique_ptr<AnalysisModule> AK4cleaner;
    std::unique_ptr<AnalysisModule> HOTVRcleaner;

    // matching, to repeat it
    std::unique_ptr<AnalysisModule> PuppiCHSMatcher;

    // object handles
    uhh2::Event::Handle<bool> h_is_btagevent;
    uhh2::Event::Handle<double> h_ST_AK4;
    uhh2::Event::Handle<double> h_ST_HOTVR;
    uhh2::Event::Handle<std::vector<Jet>> h_CHS_matched;

};
