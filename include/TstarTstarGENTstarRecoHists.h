#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"
#include "UHH2/TstarTstar/include/ReconstructionTstarHypothesis.h"

namespace uhh2{

/**  \brief Example class for booking and filling histograms
 *
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class TstarTstarGENTstarRecoHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    TstarTstarGENTstarRecoHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~TstarTstarGENTstarRecoHists();

  private:

    uhh2::Event::Handle<ReconstructionTstarHypothesis> h_tstartstar_hyp_gHOTVR;
    uhh2::Event::Handle<ReconstructionTstarHypothesis> h_tstartstar_hyp_gAK4;
};

}
