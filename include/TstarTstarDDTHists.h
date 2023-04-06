#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"

namespace uhh2{

/**  \brief Example class for booking and filling histograms
 *
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class TstarTstarDDTHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    TstarTstarDDTHists(uhh2::Context & ctx, const std::string & dirname, const std::vector<TString> points);

    virtual void fill(const uhh2::Event & ev, std::vector<double> taggerScores);
    virtual void fill(const uhh2::Event & ev) override {};
    virtual ~TstarTstarDDTHists();

  private:

    // ST handle
    uhh2::Event::Handle<double> h_ST;

    // points to plot for
    std::vector<TString> points_;

};

}
