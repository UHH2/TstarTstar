#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/PhotonIds.h"
#include <UHH2/common/include/MuonIds.h>
#include "UHH2/common/include/MCWeight.h"
#include <UHH2/common/include/TriggerSelection.h>

#include "UHH2/TstarTstar/include/B2G_TrimMass_Hists.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

/**
 *
 * Simple quick analysis module to produce trigger efficiencies.
 *
 */
class B2G_TrimMass_Trigger_Studies: public AnalysisModule {
public:

  explicit B2G_TrimMass_Trigger_Studies(Context & ctx);
  virtual bool process(Event & event) override;

private:


  // module for cleaning
  std::unique_ptr<CommonModules> common; // ak4 (not used)
  unique_ptr<TopJetCleaner> AK8Cleaner; // ak8


  // ##### triggers #####
  unique_ptr<Selection> trg_mu_ref_1, trg_mu_ref_2;
  unique_ptr<Selection> trg_Ele15_HT450;
  unique_ptr<Selection> trg_Ele35, trg_Ele115, trg_Pho200, trg_Cross;

  // ##### Histograms classes #####
  // will all be objects of the B2G_TrimMass_Hists class, containting multiple histograms each time!
  std::unique_ptr<Hists> h_before;
  std::unique_ptr<Hists> h_afterAll, h_afterWithoutCross;

};


B2G_TrimMass_Trigger_Studies::B2G_TrimMass_Trigger_Studies(Context & ctx){

  // creating the objects
  trg_mu_ref_1.reset(new TriggerSelection("HLT_Mu50_v*"));
  trg_mu_ref_2.reset(new TriggerSelection("HLT_IsoMu24_v*"));

  trg_Ele35.reset(new TriggerSelection("HLT_Ele35_WPTight_Gsf_v*"));
  trg_Ele115.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
  trg_Pho200.reset(new TriggerSelection("HLT_Photon200_v*"));

  trg_Cross.reset(new TriggerSelection("HLT_Ele15_IsoVVVL_PFHT450_v*"));

  // Cleaners
  common.reset(new CommonModules());
  MuonId muID = MuonID(Muon::CutBasedIdLoose);
  double muon_pt(30.);
  ElectronId eleID = ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp90);
  double electron_pt(30.);
  common->set_electron_id(AndId<Electron>(PtEtaSCCut(electron_pt, 2.4), eleID));
  common->set_muon_id(AndId<Muon>(PtEtaCut(muon_pt, 2.4), muID));
  common->set_jet_id(AndId<Jet>(PtEtaCut(30., 2.4), JetPFID(JetPFID::WP_TIGHT_PUPPI)));
  common->init(ctx);
  AK8Cleaner.reset(new TopJetCleaner(ctx, PtEtaCut(200.0, 2.4)));

  // histograms
  h_before.reset(new B2G_TrimMass_Hists(ctx, "LeptonCross_before"));

  h_afterAll.reset(new B2G_TrimMass_Hists(ctx, "LeptonCross_afterAll"));
  h_afterWithoutCross.reset(new B2G_TrimMass_Hists(ctx, "LeptonCross_afterWithoutCross"));

}


bool B2G_TrimMass_Trigger_Studies::process(Event & event) {

  // applying cleaners
  if(!(common->process(event))) return false;
  if(!(AK8Cleaner->process(event))) return false;

  // getting trigger booleans
  bool pass_mu_ref = trg_mu_ref_1->passes(event) || trg_mu_ref_2->passes(event);;
  bool pass_Ele35 = trg_Ele35->passes(event);
  bool pass_Ele115 = trg_Ele115->passes(event);
  bool pass_Pho200 = trg_Pho200->passes(event);
  bool pass_Cross = trg_Cross->passes(event);

  // single mu reference selection
  if(!pass_mu_ref) return false;

  // requiring a muon and an electron
  if(event.muons->size() != 1) return false;
  if(event.electrons->size() != 1) return false;

  // histogram before triggers
  h_before->fill(event);

  // histograms for various trigger (combinations)
  if(pass_Ele35 || pass_Ele115 || pass_Pho200) h_afterWithoutCross->fill(event);
  if(pass_Ele35 || pass_Ele115 || pass_Pho200 || pass_Cross) h_afterAll->fill(event);


  return false;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarPreselectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(B2G_TrimMass_Trigger_Studies)

}
