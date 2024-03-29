#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>
#include "UHH2/TstarTstar/include/ReconstructionTstarHypothesis.h"


namespace uhh2{

/**  \brief Example class for booking and filling histograms
 *
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class TstarTstarModelHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    explicit TstarTstarModelHists(uhh2::Context & ctx, const std::string & dirname);
    virtual void fill(const uhh2::Event & ev) override;

protected:
    bool is_mc;
    virtual ~TstarTstarModelHists();

    uhh2::Event::Handle<TTbarGen> h_ttbargen;
};



}
