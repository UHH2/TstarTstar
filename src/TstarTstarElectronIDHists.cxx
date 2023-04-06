#include "UHH2/TstarTstar/include/TstarTstarElectronIDHists.h"
#include "UHH2/core/include/Event.h"
#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/TTbarGen.h>
#include "UHH2/common/include/ElectronIds.h"
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


TstarTstarElectronIDHists::TstarTstarElectronIDHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
  twodcut_sel.reset(new TwoDCut(0.4, 25.0));  // The same as in Z'->ttbar semileptonic

  book<TH1F>("count_after_ttbargen", "nothing", 1, 0, 1);
  book<TH1F>("count_after_reco_clean", "nothing", 1, 0, 1);
  book<TH1F>("count_after_reco_clean_low", "nothing", 1, 0, 1);
  book<TH1F>("count_after_reco_clean_high", "nothing", 1, 0, 1);


  book<TH1F>("count_after_mvaEleID_Fall17_iso_V2_wp90", "nothing", 1, 0, 1);
  book<TH1F>("count_after_mvaEleID_Fall17_noIso_V2_wp90", "nothing", 1, 0, 1);
  book<TH1F>("count_after_mvaEleID_Fall17_iso_V2_wp80", "nothing", 1, 0, 1);
  book<TH1F>("count_after_mvaEleID_Fall17_noIso_V2_wp80", "nothing", 1, 0, 1);
  book<TH1F>("count_after_heepElectronID_HEEPV70", "nothing", 1, 0, 1);
  book<TH1F>("count_after_cutBasedElectronID_Fall17_94X_V2_loose", "nothing", 1, 0, 1);
  book<TH1F>("count_after_cutBasedElectronID_Fall17_94X_V2_medium", "nothing", 1, 0, 1);
  book<TH1F>("count_after_cutBasedElectronID_Fall17_94X_V2_tight", "nothing", 1, 0, 1);

  book<TH1F>("count_after_mvaEleID_Fall17_iso_V2_wp90_low", "nothing", 1, 0, 1);
  book<TH1F>("count_after_mvaEleID_Fall17_noIso_V2_wp90_low", "nothing", 1, 0, 1);
  book<TH1F>("count_after_mvaEleID_Fall17_iso_V2_wp80_low", "nothing", 1, 0, 1);
  book<TH1F>("count_after_mvaEleID_Fall17_noIso_V2_wp80_low", "nothing", 1, 0, 1);
  book<TH1F>("count_after_heepElectronID_HEEPV70_low", "nothing", 1, 0, 1);
  book<TH1F>("count_after_cutBasedElectronID_Fall17_94X_V2_loose_low", "nothing", 1, 0, 1);
  book<TH1F>("count_after_cutBasedElectronID_Fall17_94X_V2_medium_low", "nothing", 1, 0, 1);
  book<TH1F>("count_after_cutBasedElectronID_Fall17_94X_V2_tight_low", "nothing", 1, 0, 1);

  book<TH1F>("count_after_mvaEleID_Fall17_iso_V2_wp90_high", "nothing", 1, 0, 1);
  book<TH1F>("count_after_mvaEleID_Fall17_noIso_V2_wp90_high", "nothing", 1, 0, 1);
  book<TH1F>("count_after_mvaEleID_Fall17_iso_V2_wp80_high", "nothing", 1, 0, 1);
  book<TH1F>("count_after_mvaEleID_Fall17_noIso_V2_wp80_high", "nothing", 1, 0, 1);
  book<TH1F>("count_after_heepElectronID_HEEPV70_high", "nothing", 1, 0, 1);
  book<TH1F>("count_after_cutBasedElectronID_Fall17_94X_V2_loose_high", "nothing", 1, 0, 1);
  book<TH1F>("count_after_cutBasedElectronID_Fall17_94X_V2_medium_high", "nothing", 1, 0, 1);
  book<TH1F>("count_after_cutBasedElectronID_Fall17_94X_V2_tight_high", "nothing", 1, 0, 1);

  is_mc = ctx.get("dataset_type") == "MC";
}


void TstarTstarElectronIDHists::fill(const Event & event){

  double weight = event.weight;

  assert(event.genparticles);

  // get ttbargen
  TTbarGen ttbargen = event.get(h_ttbargen);
  if(!ttbargen.IsSemiLeptonicDecay()) return; // only semileptonic
  if(abs(ttbargen.ChargedLepton().pdgId()) != 11) return; // only electron

  hist("count_after_ttbargen")->Fill(0.5, weight);

  // only look at events that have any electrons above 40 GeV
  if(event.electrons->size()==0) return;
  if(event.electrons->at(0).pt() < 40) return;
  if(event.electrons->at(0).eta() > 2.4) return;

  // match leading electron to GEN electron
  if(deltaR(ttbargen.ChargedLepton(), event.electrons->at(0)) > 0.2) return;

  hist("count_after_reco_clean")->Fill(0.5, weight);

  if(  ElectronTagID(Electron::mvaEleID_Fall17_iso_V2_wp90)(event.electrons->at(0), event)) hist("count_after_mvaEleID_Fall17_iso_V2_wp90")->Fill(0.5, weight);
  if(  ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp90)(event.electrons->at(0), event)) hist("count_after_mvaEleID_Fall17_noIso_V2_wp90")->Fill(0.5, weight);
  if(  ElectronTagID(Electron::mvaEleID_Fall17_iso_V2_wp80)(event.electrons->at(0), event)) hist("count_after_mvaEleID_Fall17_iso_V2_wp80")->Fill(0.5, weight);
  if(  ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp80)(event.electrons->at(0), event)) hist("count_after_mvaEleID_Fall17_noIso_V2_wp80")->Fill(0.5, weight);
  if(  ElectronTagID(Electron::heepElectronID_HEEPV70)(event.electrons->at(0), event)) hist("count_after_heepElectronID_HEEPV70")->Fill(0.5, weight);
  if(  ElectronTagID(Electron::cutBasedElectronID_Fall17_94X_V2_loose)(event.electrons->at(0), event)) hist("count_after_cutBasedElectronID_Fall17_94X_V2_loose")->Fill(0.5, weight);
  if(  ElectronTagID(Electron::cutBasedElectronID_Fall17_94X_V2_medium)(event.electrons->at(0), event)) hist("count_after_cutBasedElectronID_Fall17_94X_V2_medium")->Fill(0.5, weight);
  if(  ElectronTagID(Electron::cutBasedElectronID_Fall17_94X_V2_tight)(event.electrons->at(0), event)) hist("count_after_cutBasedElectronID_Fall17_94X_V2_tight")->Fill(0.5, weight);  

  if(event.electrons->at(0).pt() > 120) {

    hist("count_after_reco_clean_high")->Fill(0.5, weight);

    if(  ElectronTagID(Electron::mvaEleID_Fall17_iso_V2_wp90)(event.electrons->at(0), event)) hist("count_after_mvaEleID_Fall17_iso_V2_wp90_high")->Fill(0.5, weight);
    if(  ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp90)(event.electrons->at(0), event)) hist("count_after_mvaEleID_Fall17_noIso_V2_wp90_high")->Fill(0.5, weight);
    if(  ElectronTagID(Electron::mvaEleID_Fall17_iso_V2_wp80)(event.electrons->at(0), event)) hist("count_after_mvaEleID_Fall17_iso_V2_wp80_high")->Fill(0.5, weight);
    if(  ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp80)(event.electrons->at(0), event)) hist("count_after_mvaEleID_Fall17_noIso_V2_wp80_high")->Fill(0.5, weight);
    if(  ElectronTagID(Electron::heepElectronID_HEEPV70)(event.electrons->at(0), event)) hist("count_after_heepElectronID_HEEPV70_high")->Fill(0.5, weight);
    if(  ElectronTagID(Electron::cutBasedElectronID_Fall17_94X_V2_loose)(event.electrons->at(0), event)) hist("count_after_cutBasedElectronID_Fall17_94X_V2_loose_high")->Fill(0.5, weight);
    if(  ElectronTagID(Electron::cutBasedElectronID_Fall17_94X_V2_medium)(event.electrons->at(0), event)) hist("count_after_cutBasedElectronID_Fall17_94X_V2_medium_high")->Fill(0.5, weight);
    if(  ElectronTagID(Electron::cutBasedElectronID_Fall17_94X_V2_tight)(event.electrons->at(0), event)) hist("count_after_cutBasedElectronID_Fall17_94X_V2_tight_high")->Fill(0.5, weight);  

  } else {

    hist("count_after_reco_clean_low")->Fill(0.5, weight);

    if(  ElectronTagID(Electron::mvaEleID_Fall17_iso_V2_wp90)(event.electrons->at(0), event)) hist("count_after_mvaEleID_Fall17_iso_V2_wp90_low")->Fill(0.5, weight);
    if(  ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp90)(event.electrons->at(0), event)) hist("count_after_mvaEleID_Fall17_noIso_V2_wp90_low")->Fill(0.5, weight);
    if(  ElectronTagID(Electron::mvaEleID_Fall17_iso_V2_wp80)(event.electrons->at(0), event)) hist("count_after_mvaEleID_Fall17_iso_V2_wp80_low")->Fill(0.5, weight);
    if(  ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp80)(event.electrons->at(0), event)) hist("count_after_mvaEleID_Fall17_noIso_V2_wp80_low")->Fill(0.5, weight);
    if(  ElectronTagID(Electron::heepElectronID_HEEPV70)(event.electrons->at(0), event)) hist("count_after_heepElectronID_HEEPV70_low")->Fill(0.5, weight);
    if(  ElectronTagID(Electron::cutBasedElectronID_Fall17_94X_V2_loose)(event.electrons->at(0), event)) hist("count_after_cutBasedElectronID_Fall17_94X_V2_loose_low")->Fill(0.5, weight);
    if(  ElectronTagID(Electron::cutBasedElectronID_Fall17_94X_V2_medium)(event.electrons->at(0), event)) hist("count_after_cutBasedElectronID_Fall17_94X_V2_medium_low")->Fill(0.5, weight);
    if(  ElectronTagID(Electron::cutBasedElectronID_Fall17_94X_V2_tight)(event.electrons->at(0), event)) hist("count_after_cutBasedElectronID_Fall17_94X_V2_tight_low")->Fill(0.5, weight);  

  }

}

TstarTstarElectronIDHists::~TstarTstarElectronIDHists(){}
