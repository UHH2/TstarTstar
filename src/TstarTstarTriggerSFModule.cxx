#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>

#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/CommonModules.h>
#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/TriggerSelection.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/TopPtReweight.h>

#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/LuminosityHists.h>
#include <UHH2/common/include/JetHists.h>
#include "UHH2/TstarTstar/include/CorrectionFactor.h"
#include "UHH2/TstarTstar/include/ElecTriggerSF.h"

/*
*******************************************************************
**************** TO DO ********************************************
*******************************************************************
- b tagging SF
*******************************************************************
*******************************************************************
*/

class TstarTstarTriggerSFModule : public AnalysisModule {

public:
  explicit TstarTstarTriggerSFModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

protected:
  enum lepton { muon, elec };
  lepton channel_;

  // selections

  std::unique_ptr<uhh2::Selection> trigger_mu_A;
  std::unique_ptr<uhh2::Selection> trigger_mu_B;
  std::unique_ptr<uhh2::Selection> trigger_el_A;
  std::unique_ptr<uhh2::Selection> trigger_el_B;
  std::unique_ptr<uhh2::Selection> trigger_el_C;

  Event::Handle<bool>h_recsel;
  Event::Handle<double>h_pt;
  Event::Handle<double>h_eta;
  Event::Handle<double>h_weight;
  Event::Handle<double>h_weight_SFpt;
  Event::Handle<double>h_weight_SFeta;
  Event::Handle<double>h_weight_SFetapt;
  Event::Handle<double>h_weight_SFetaptUP;
  Event::Handle<double>h_weight_SFetaptDOWN;
  Event::Handle<bool>h_passed;
  Event::Handle<bool>h_passed_elec;
  Event::Handle<bool>h_passed_photon;
  Event::Handle<int>h_run;
  Event::Handle<int>h_lumi;
  Event::Handle<int>h_eventnr;

  std::unique_ptr<uhh2::AnalysisModule> ele_id_SF, ele_trigger_SFpt, ele_trigger_SFeta, ele_reco_SF;
  std::unique_ptr<uhh2::AnalysisModule> ele_trigger_SFetapt, ele_trigger_SFetaptUP, ele_trigger_SFetaptDOWN;
  std::unique_ptr<uhh2::AnalysisModule> muo_tight_noniso_SF, muo_trigger_SF, muo_trigger_SF_B;

  bool debug = false;
  bool isMC; //define here to use it in "process" part
  bool year_16;
  bool year_17;
  bool year_18;
  Year year;

  bool data_is2017B = false;
  uhh2::Event::Handle<bool> h_MC_isfake2017B;

  std::unique_ptr<Hists> h_pass, h_all;
};

TstarTstarTriggerSFModule::TstarTstarTriggerSFModule(uhh2::Context& ctx){

  if(debug) cout << "Get Year ... " << endl;
  year_16 = false;
  year_17 = false;
  year_18 = false;
  year = extract_year(ctx);

  if(year == Year::is2016v3 or year == Year::isUL16preVFP or year == Year::isUL16postVFP) year_16 = true;
  else if(year == Year::is2017v2 or year == Year::isUL17) year_17 = true;
  else if(year == Year::is2018 or year == Year::isUL18) year_18 = true;
  else throw runtime_error("In PostSelectionModule: This Event is not from 2016v3, 2017v2 or 2018!");

  isMC = (ctx.get("dataset_type") == "MC");

  if(debug) cout << "Declare Output ... " << endl;
  ctx.undeclare_all_event_output();
  h_eta = ctx.declare_event_output<double>("eta");
  h_pt = ctx.declare_event_output<double>("pt");
  h_weight = ctx.declare_event_output<double>("weight");
  h_passed = ctx.declare_event_output<bool>("passed");
  h_run = ctx.declare_event_output<int>("run");
  h_lumi = ctx.declare_event_output<int>("lumi");
  h_eventnr = ctx.declare_event_output<int>("eventnr");

  if(debug) cout << "Define Trigger ... " << endl;
  if(year_16)      trigger_el_A = uhh2::make_unique<TriggerSelection>("HLT_Ele27_WPTight_Gsf_v*");
  else if(year_17) trigger_el_A = uhh2::make_unique<TriggerSelection>("HLT_Ele35_WPTight_Gsf_v*");
  else if(year_18) trigger_el_A = uhh2::make_unique<TriggerSelection>("HLT_Ele32_WPTight_Gsf_v*");
  trigger_el_B = uhh2::make_unique<TriggerSelection>("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*");
  if(year_16) trigger_el_C = uhh2::make_unique<TriggerSelection>("HLT_Photon175_v*");
  else        trigger_el_C = uhh2::make_unique<TriggerSelection>("HLT_Photon200_v*");

  h_pass.reset(new ElectronHists(ctx, "pass_Elec"));
  h_all.reset(new ElectronHists(ctx, "all_Elec"));

  h_MC_isfake2017B = ctx.get_handle<bool>("MC_isfake2017B");

  if(!isMC) data_is2017B = (ctx.get("dataset_version").find("SingleElectron_RunB_UL17") != std::string::npos) || (ctx.get("dataset_version").find("SingleMuon_RunB_UL17") != std::string::npos) || (ctx.get("dataset_version").find("SinglePhoton_RunB_UL17") != std::string::npos);

}


bool TstarTstarTriggerSFModule::process(uhh2::Event& event){
  if(debug) cout << " ------------------------------------------- " << endl;
  if(debug) cout << " ----------------- NewEvent ---------------- " << endl;
  if(debug) cout << " ------------------------------------------- " << endl;

  // fill hists with all events
  h_all->fill(event);

  // HERE FILL PT AND ETA HISTS FOR PASSING AND NOT PASSING ELEC TRIGGER
  if(debug) cout << "Start Fill ... " << endl;
  bool passed_elec_trigger = false;
  if(year_16) passed_elec_trigger = (trigger_el_A->passes(event) || trigger_el_B->passes(event) || trigger_el_C->passes(event));
  if(year_17){
    // for MC event.run=1
    if(data_is2017B || event.get(h_MC_isfake2017B)) passed_elec_trigger = (trigger_el_A->passes(event) || trigger_el_C->passes(event));
    else                                            passed_elec_trigger = (trigger_el_A->passes(event) || trigger_el_B->passes(event) || trigger_el_C->passes(event));
  }
  if(year_18)  passed_elec_trigger = (trigger_el_A->passes(event) || trigger_el_B->passes(event) || trigger_el_C->passes(event));

  if(debug) cout << "Set after Trigger pass... " << endl;
  event.set(h_pt, event.electrons->at(0).pt());
  event.set(h_eta, event.electrons->at(0).eta());
  event.set(h_weight, event.weight);
  event.set(h_passed, passed_elec_trigger);

  int run = 0;
  int lumi = 0;
  int eventnr = 0;

  if(!isMC){
    run = event.run;
    lumi = event.luminosityBlock;
    eventnr = event.event;
  }

  if(debug) cout << "Set run info ... " << endl;
  event.set(h_run, run);
  event.set(h_lumi, lumi);
  event.set(h_eventnr, eventnr);

  if(passed_elec_trigger){
    // fill pass histograms
    h_pass->fill(event);
  }

  // if(debug) cout << "Event done ... " << endl;
  if(debug) cout << "Event done ... " << endl;
  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarTriggerSFModule)
