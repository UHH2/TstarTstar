#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ReconstructionHypothesis.h"

#include "UHH2/TstarTstar/include/ReconstructionTstarHypothesis.h"

namespace uhh2{

  /**  \brief Example class for booking and filling histograms
   * 
   * NOTE: This class uses the 'hist' method to retrieve histograms.
   * This requires a string lookup and is therefore slow if you have
   * many histograms. Therefore, it is recommended to use histogram
   * pointers as member data instead, like in 'common/include/ElectronHists.h'.
   */
  class TstarTstarRecoTstarHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    TstarTstarRecoTstarHists(uhh2::Context & ctx, const std::string & dirname);
    
    virtual void fill(const uhh2::Event & ev) override;
    virtual ~TstarTstarRecoTstarHists();
    void fill_ttbarhyps(const Event & event, const ReconstructionHypothesis & hyp);
    
  private:
    
    uhh2::Event::Handle<ReconstructionTstarHypothesis> h_tstartstar_hyp;
    uhh2::Event::Handle<std::vector<ReconstructionTstarHypothesis>> h_tstartstar_hyp_vector;
    
    bool is_tgtg;
    bool is_tgtgamma;

  };

}
