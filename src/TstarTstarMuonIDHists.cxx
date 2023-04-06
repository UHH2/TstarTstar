#include "UHH2/TstarTstar/include/TstarTstarMuonIDHists.h"
#include "UHH2/core/include/Event.h"
#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/TTbarGen.h>
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

namespace {
  // invariant mass of a lorentzVector, but safe for timelike / spacelike vectors
  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }
}


TstarTstarMuonIDHists::TstarTstarMuonIDHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
  twodcut_sel.reset(new TwoDCut(0.4, 25.0));  // The same as in Z'->ttbar semileptonic


  book<TH1F>("count_beforettbargen", "N mu", 1, 0, 1);
  book<TH1F>("count_beforettbargen_after", "N mu", 1, 0, 1);

  book<TH1F>("count_ttbargen", "N mu", 1, 0, 1);
  book<TH1F>("count_ttbargen_after", "N mu", 1, 0, 1);
  book<TH1F>("sel_matched_to_gen", "matched?", 1, 0, 1);

  is_mc = ctx.get("dataset_type") == "MC";
}


void TstarTstarMuonIDHists::fill(const Event & event){

  assert(event.genparticles);

  double dR_cut_low = -1;
  double dR_cut_high = -1;
  double do2Dcut = false;
  double doAddMETcut = true;

  // ###### Lepton-2Dcut ######
  bool pass_2D = false;
  if(event.muons->size()==1) {
    for(auto& muo : *event.muons){
      float    dRmin, pTrel;
      std::tie(dRmin, pTrel) = drmin_pTrel(muo, *event.jets);
      muo.set_tag(Muon::twodcut_dRmin, dRmin);
      muo.set_tag(Muon::twodcut_pTrel, pTrel);
    }
    for(auto& ele : *event.electrons){
      float    dRmin, pTrel;
      std::tie(dRmin, pTrel) = drmin_pTrel(ele, *event.jets);
      ele.set_tag(Electron::twodcut_dRmin, dRmin);
      ele.set_tag(Electron::twodcut_pTrel, pTrel);
    }
    const bool pass_twodcut = twodcut_sel->passes(event);
    if(event.muons->at(0).pt()>120) pass_2D = pass_twodcut;
    else pass_2D = true;
  }

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  double dR_cut = 99999;

  hist("count_beforettbargen")->Fill(0.5, weight);
  if(event.muons->size() == 1) {
    double dR = 9999;
    for (const auto & jet : *event.jets) {
      double dR_tmp = deltaR(jet, event.muons->at(0));
      if(dR_tmp < dR) dR = dR_tmp;
    }
    if(event.muons->at(0).pt() > 120) dR_cut = dR_cut_high;
    else dR_cut = dR_cut_low;
    if((dR > dR_cut) && (!do2Dcut || (do2Dcut && pass_2D)) && (!doAddMETcut || (doAddMETcut && (event.met->pt() > 40)))) hist("count_beforettbargen_after")->Fill(0.5, weight);
  }

  // get ttbargen
  TTbarGen ttbargen = event.get(h_ttbargen);
  if(!ttbargen.IsSemiLeptonicDecay()) return;

  hist("count_ttbargen")->Fill(0.5, weight);
  if(event.muons->size() == 1) {
    double dR = 9999;
    for (const auto & jet : *event.jets) {
      double dR_tmp = deltaR(jet, event.muons->at(0));
      if(dR_tmp < dR) dR = dR_tmp;
    }
    if(event.muons->at(0).pt() > 120) dR_cut = dR_cut_high;
    else dR_cut = dR_cut_low;
    if ((dR > dR_cut) && (!do2Dcut || (do2Dcut && pass_2D)) && (!doAddMETcut || (doAddMETcut && (event.met->pt() > 40)))) {
      hist("count_ttbargen_after")->Fill(0.5, weight); // counting selection efficiency
      if(deltaR(ttbargen.ChargedLepton(), event.muons->at(0)) < 0.2) hist("sel_matched_to_gen")->Fill(1., weight);
      else hist("sel_matched_to_gen")->Fill(0., weight);
    }
  }



}

TstarTstarMuonIDHists::~TstarTstarMuonIDHists(){}
