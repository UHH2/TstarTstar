#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/common/include/ReconstructionHypothesis.h"
#include "UHH2/common/include/TopJetIds.h"


#include <string>
#include <vector>
#include <unordered_map>

namespace uhh2 {
  class TwoDCut : public Selection {
  public:
    explicit TwoDCut(float min_deltaR, float min_pTrel): min_deltaR_(min_deltaR), min_pTrel_(min_pTrel) {}
    virtual bool passes(const Event&) override;

  private:
    float min_deltaR_, min_pTrel_;
  };
  ////

  class TTbarSemiLepMatchableSelection: public uhh2::Selection{
  public:
    TTbarSemiLepMatchableSelection();
    ~TTbarSemiLepMatchableSelection(){};
    virtual bool passes(const uhh2::Event & event);
    std::pair<bool,double> check_reco(const ReconstructionHypothesis hyp);//compares match between reconstructed hypothesis vs gen tops and products of their decays filled in passes()
  private:
    GenParticle Wlep, Whad, blep, bhad, thad, tlep, lepton, neutrino, Whadd1,Whadd2;
  };




  class METCut : public Selection {

  public:
    explicit METCut(float, float max_met=1e7);
    virtual bool passes(const Event&) override;

  private:
    float min_met_, max_met_;
  };
  ////

  class STCut : public Selection {

  public:
    explicit STCut(float, float max_st=1e7);
    virtual bool passes(const Event&) override;

  private:
    float min_st_, max_st_;
  };
  ////

  class TopTagEventSelection: public Selection {
  public:
    explicit TopTagEventSelection(const TopJetId& tjet_id=CMSTopTag());
    virtual bool passes(const Event&) override;

  private:
    TopJetId topjetID_;
  };
  /////

  class Eled0dzCut {
  public:
    Eled0dzCut(float max_d0_, float max_dz_): max_d0(max_d0_), max_dz(max_dz_){}

      bool operator()(const Electron & ele, const uhh2::Event & event) const{
        const auto & pv = event.pvs->at(0);
        double pv_d0 = ele.gsfTrack_dxy_vertex(pv.x(), pv.y());
        double pv_dz = ele.gsfTrack_dz_vertex(pv.x(), pv.y(), pv.z());
        return (pv_d0 < max_d0 && pv_dz < max_dz);
      }

  private:
    float max_d0, max_dz;
  };

}
