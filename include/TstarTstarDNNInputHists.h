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
class TstarTstarDNNInputHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    TstarTstarDNNInputHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~TstarTstarDNNInputHists();

  private:
    uhh2::Event::Handle<std::vector<double>> h_DNN_Inputs;

    TH1D *DNN_2D_1, *DNN_2D_2, *DNN_2D_3, *DNN_2D_4, *DNN_2D_5, *DNN_2D_6, *DNN_2D_7, *DNN_2D_8, *DNN_2D_9, *DNN_2D_10;
    TH1D *DNN_2D_11, *DNN_2D_12, *DNN_2D_13, *DNN_2D_14, *DNN_2D_15, *DNN_2D_16, *DNN_2D_17, *DNN_2D_18, *DNN_2D_19, *DNN_2D_20;
    TH1D *DNN_2D_21, *DNN_2D_22, *DNN_2D_23, *DNN_2D_24, *DNN_2D_25, *DNN_2D_26, *DNN_2D_27, *DNN_2D_28, *DNN_2D_29, *DNN_2D_30;
    TH1D *DNN_2D_31, *DNN_2D_32, *DNN_2D_33, *DNN_2D_34, *DNN_2D_35, *DNN_2D_36, *DNN_2D_37, *DNN_2D_38, *DNN_2D_39, *DNN_2D_40;
    TH1D *DNN_2D_41, *DNN_2D_42, *DNN_2D_43, *DNN_2D_44, *DNN_2D_45, *DNN_2D_46, *DNN_2D_47, *DNN_2D_48, *DNN_2D_49, *DNN_2D_50;
    TH1D *DNN_2D_51, *DNN_2D_52, *DNN_2D_53, *DNN_2D_54;

};

}