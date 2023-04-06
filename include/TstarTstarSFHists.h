#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"

namespace uhh2{

/**  \brief Example class for booking and filling histograms
 *
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class TstarTstarSFHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    TstarTstarSFHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~TstarTstarSFHists();

  private:

    uhh2::Event::Handle<float> h_weight_sfmu_id;
    uhh2::Event::Handle<float> h_weight_sfmu_id_down;
    uhh2::Event::Handle<float> h_weight_sfmu_id_up;

    uhh2::Event::Handle<float> h_weight_sfmu_isolation;
    uhh2::Event::Handle<float> h_weight_sfmu_isolation_down;
    uhh2::Event::Handle<float> h_weight_sfmu_isolation_up;

    uhh2::Event::Handle<float> h_weight_sfele_id;
    uhh2::Event::Handle<float> h_weight_sfele_id_down;
    uhh2::Event::Handle<float> h_weight_sfele_id_up;

};

}
