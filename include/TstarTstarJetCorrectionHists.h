#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"

namespace uhh2{

    class TstarTstarJetCorrectionHists: public uhh2::Hists {
    public:
        explicit TstarTstarJetCorrectionHists(uhh2::Context & ctx, const std::string & dirname);
        virtual void fill(const uhh2::Event & ev) override;

    protected:
        virtual ~TstarTstarJetCorrectionHists();

        double matching_radius;
    };

}
