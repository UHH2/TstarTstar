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

#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarElectronIDHists.h"
#include "UHH2/TstarTstar/include/TstarTstarMuonIDHists.h"
#include "UHH2/TstarTstar/include/TstarTstarCustomIds.h"

#include "UHH2/HOTVR/include/HOTVRJetCorrectionModule.h"



using namespace std;
using namespace uhh2;

namespace uhh2 {

/** bla
 *
 * blub
 *
 */
class TstarTstarGenStudies: public AnalysisModule {
public:

  explicit TstarTstarGenStudies(Context & ctx);
  virtual bool process(Event & event) override;

private:

  // ##### Histograms #####
  // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_GenHists;
  std::unique_ptr<Hists> h_ElectronIDHists_beforeCut;
  std::unique_ptr<Hists> h_ElectronIDHists_afterCut;
  std::unique_ptr<Hists> h_MuonIDHists_beforeCut;
  std::unique_ptr<Hists> h_MuonIDHists_afterCut;

  unique_ptr<Selection> met_sel;

  uhh2::Event::Handle<TTbarGen> h_ttbargen;

  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<AnalysisModule> HOTVRCorr;
  std::unique_ptr<TopJetCleaner> HOTVRcleaner;


};


TstarTstarGenStudies::TstarTstarGenStudies(Context & ctx){

  for(auto & kv : ctx.get_all()){
    cout << " " << kv.first << " = " << kv.second << endl;
  }

  h_GenHists.reset(new TstarTstarGenHists(ctx, "GenHists"));
  h_ElectronIDHists_beforeCut.reset(new TstarTstarElectronIDHists(ctx, "ElectronIDHists_beforeCut"));
  h_ElectronIDHists_afterCut.reset(new TstarTstarElectronIDHists(ctx, "ElectronIDHists_afterCut"));
  h_MuonIDHists_beforeCut.reset(new TstarTstarMuonIDHists(ctx, "MuonIDHists_beforeCut"));
  h_MuonIDHists_afterCut.reset(new TstarTstarMuonIDHists(ctx, "MuonIDHists_afterCut"));

  common.reset(new CommonModules());

  // electron ID
  ElectronId eleID_lowpt = ElectronTagID(Electron::mvaEleID_Fall17_iso_V2_wp90);
  ElectronId eleID_highpt = ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp90);
  //ElectronId eleID_highpt = AndId<Electron>(ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp90), Electron_MINIIso(0.1, "pf-weight"));
  double electron_pt_lowpt(55.); // was 30
  double electron_pt_highpt(120.);
  common->set_electron_id(OrId<Electron>( AndId<Electron>(PtEtaSCCut(electron_pt_lowpt, 2.4), eleID_lowpt, EleMaxPtCut(120.)),  AndId<Electron>(PtEtaSCCut(electron_pt_highpt, 2.4), eleID_highpt)));

  MuonId muID_lowpt = AndId<Muon>(MuonID(Muon::CutBasedIdTight), MuonID(Muon::PFIsoTight));
  MuonId muID_highpt = MuonID(Muon::CutBasedIdTight);
  double muon_pt_lowpt(27.);
  double muon_pt_highpt(60.);
  common->set_muon_id(OrId<Muon>( AndId<Muon>(PtEtaCut(muon_pt_lowpt, 2.4), muID_lowpt, MuMaxPtCut(60.)),  AndId<Muon>(PtEtaCut(muon_pt_highpt, 2.4), muID_highpt)));

  common->switch_jetlepcleaner();
  common->switch_jetPtSorter();
  double jet_pt(30.);
  common->set_jet_id(AndId<Jet>(PtEtaCut(jet_pt, 2.5), JetPFID(JetPFID::WP_TIGHT_PUPPI)));

  // HOTVR jets
  HOTVRCorr.reset(new HOTVRJetCorrectionModule(ctx)); // crashes
  TopJetId topjetID = AndId<TopJet>(HOTVRTopTag(), Tau32Groomed(0.56)); // Top Tag that is used later
  HOTVRcleaner.reset(new TopJetCleaner(ctx, PtEtaCut(150.0, 2.5)));

  common->init(ctx);

  met_sel.reset(new METCut  (50.,1e9));

  ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

}


bool TstarTstarGenStudies::process(Event & event) {

  if(!(common->process(event))) return false;
  if(!(HOTVRCorr->process(event))) return false;
  if(!(HOTVRcleaner->process(event))) return false;

  // apply selection cuts
  if(event.jets->size() < 4) return false;
  if(event.topjets->size() < 1) return false;
  bool pass_btagcut = false;
  for (const auto & jet: *event.jets){
    if(jet.btag_DeepCSV() > 0.2219) pass_btagcut = true;
  }
  if(!pass_btagcut) return false;
  if(!met_sel->passes(event)) return false;

  double st = 0.;
  for(const auto & jet : *event.topjets) st += jet.pt();
  for(const auto & lepton : *event.electrons) st += lepton.pt();
  for(const auto & lepton : *event.muons) st += lepton.pt();
  st += event.met->pt();
  if(st < 500) return false;


  ttgenprod->process(event);
  TTbarGen ttbargen = event.get(h_ttbargen);

  h_ElectronIDHists_beforeCut->fill(event);
  h_MuonIDHists_beforeCut->fill(event);

  if(!ttbargen.IsSemiLeptonicDecay()) return false;

  h_GenHists->fill(event);

  if(abs(ttbargen.ChargedLepton().pdgId()) == 11) {
    if(ttbargen.ChargedLepton().pt() > 30 && abs(ttbargen.ChargedLepton().eta()) < 2.4) h_ElectronIDHists_afterCut->fill(event);
  }
  else if(abs(ttbargen.ChargedLepton().pdgId()) == 13) {
    if(ttbargen.ChargedLepton().pt() > 30 && abs(ttbargen.ChargedLepton().eta()) < 2.4) h_MuonIDHists_afterCut->fill(event);
  }


  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarPreselectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarGenStudies)

}
