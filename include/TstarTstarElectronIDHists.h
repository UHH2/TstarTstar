#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include <UHH2/common/include/TTbarGen.h>
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"

namespace uhh2{

/**  \brief Example class for booking and filling histograms
 *
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class TstarTstarElectronIDHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    explicit TstarTstarElectronIDHists(uhh2::Context & ctx, const std::string & dirname);
    virtual void fill(const uhh2::Event & ev) override;

  protected:
    bool is_mc;
    virtual ~TstarTstarElectronIDHists();

    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    std::unique_ptr<Selection> twodcut_sel;

};

}
