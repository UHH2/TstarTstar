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
class TstarTstarPDFNormHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    TstarTstarPDFNormHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~TstarTstarPDFNormHists();

  private:

    bool needsOtherMCweightHandling;

    // hist names for PDF
    std::string hist_names[100];

    uhh2::Event::Handle<float> h_murmuf_upup;
    uhh2::Event::Handle<float> h_murmuf_upnone;
    uhh2::Event::Handle<float> h_murmuf_noneup;
    uhh2::Event::Handle<float> h_murmuf_nonedown;
    uhh2::Event::Handle<float> h_murmuf_downnone;
    uhh2::Event::Handle<float> h_murmuf_downdown;

};

}
