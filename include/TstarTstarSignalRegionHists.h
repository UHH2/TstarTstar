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
class TstarTstarSignalRegionHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    TstarTstarSignalRegionHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~TstarTstarSignalRegionHists();

  private:

    // ST handle
    uhh2::Event::Handle<double> h_ST_HOTVR;
    bool is_MC;
    bool needsOtherMCweightHandling;

    // histogram containing decorrelation uncertainty
    TH1* decorrelationUncertainty;

    // debug bool
    bool debug = false;

    // hist names for PDF
    std::string hist_names[100];

    std::vector<TString> variations;

    // weight handles
    uhh2::Event::Handle<float> h_weight_ttbar;

    uhh2::Event::Handle<float> h_weight_puNominal;
    uhh2::Event::Handle<float> h_weight_puUp;
    uhh2::Event::Handle<float> h_weight_puDown;

    uhh2::Event::Handle<float> h_prefiringWeightNominal;
    uhh2::Event::Handle<float> h_prefiringWeightUp;
    uhh2::Event::Handle<float> h_prefiringWeightDown;

    uhh2::Event::Handle<float> h_weight_btagdiscNominal;

    uhh2::Event::Handle<float> h_weight_btagdisc_lfUp;
    uhh2::Event::Handle<float> h_weight_btagdisc_lfDown;

    uhh2::Event::Handle<float> h_weight_btagdisc_hfUp;
    uhh2::Event::Handle<float> h_weight_btagdisc_hfDown;

    uhh2::Event::Handle<float> h_weight_btagdisc_hfstats1Up;
    uhh2::Event::Handle<float> h_weight_btagdisc_hfstats1Down;

    uhh2::Event::Handle<float> h_weight_btagdisc_hfstats2Up;
    uhh2::Event::Handle<float> h_weight_btagdisc_hfstats2Down;

    uhh2::Event::Handle<float> h_weight_btagdisc_lfstats1Up;
    uhh2::Event::Handle<float> h_weight_btagdisc_lfstats1Down;

    uhh2::Event::Handle<float> h_weight_btagdisc_lfstats2Up;
    uhh2::Event::Handle<float> h_weight_btagdisc_lfstats2Down;

    uhh2::Event::Handle<float> h_weight_btagdisc_cferr1Up;
    uhh2::Event::Handle<float> h_weight_btagdisc_cferr1Down;

    uhh2::Event::Handle<float> h_weight_btagdisc_cferr2Up;
    uhh2::Event::Handle<float> h_weight_btagdisc_cferr2Down;

    uhh2::Event::Handle<float> h_weight_sfelec_idNominal;
    uhh2::Event::Handle<float> h_weight_sfelec_idUp;
    uhh2::Event::Handle<float> h_weight_sfelec_idDown;

    uhh2::Event::Handle<float> h_weight_sfelec_triggerNominal;
    uhh2::Event::Handle<float> h_weight_sfelec_triggerUp;
    uhh2::Event::Handle<float> h_weight_sfelec_triggerDown;

    uhh2::Event::Handle<float> h_weight_sfelec_recoNominal;
    uhh2::Event::Handle<float> h_weight_sfelec_recoUp;
    uhh2::Event::Handle<float> h_weight_sfelec_recoDown;

    uhh2::Event::Handle<float> h_weight_sfmu_idNominal;
    uhh2::Event::Handle<float> h_weight_sfmu_idUp;
    uhh2::Event::Handle<float> h_weight_sfmu_idDown;

    uhh2::Event::Handle<float> h_weight_sfmu_isoNominal;
    uhh2::Event::Handle<float> h_weight_sfmu_isoUp;
    uhh2::Event::Handle<float> h_weight_sfmu_isoDown;

    uhh2::Event::Handle<float> h_weight_sfmu_triggerNominal;
    uhh2::Event::Handle<float> h_weight_sfmu_triggerUp;
    uhh2::Event::Handle<float> h_weight_sfmu_triggerDown;

    uhh2::Event::Handle<float> h_murmuf_upup;
    uhh2::Event::Handle<float> h_murmuf_upnone;
    uhh2::Event::Handle<float> h_murmuf_noneup;
    uhh2::Event::Handle<float> h_murmuf_nonedown;
    uhh2::Event::Handle<float> h_murmuf_downnone;
    uhh2::Event::Handle<float> h_murmuf_downdown;

};

}
