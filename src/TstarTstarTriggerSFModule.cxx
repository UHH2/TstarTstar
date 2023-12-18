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


class TstarTstarTriggerSFModule : public AnalysisModule {

public:
  explicit TstarTstarTriggerSFModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

protected:
  enum lepton { muon, elec };
  lepton channel_;

  Event::Handle<bool> h_trigger_decision;
  Event::Handle<bool> h_trigger_decision_ele;

  Event::Handle<double>h_pt;
  Event::Handle<double>h_eta;
  Event::Handle<double>h_weight;
  Event::Handle<bool>h_passed;
  Event::Handle<int>h_run;
  Event::Handle<int>h_lumi;
  Event::Handle<int>h_eventnr;

  bool debug = false;
  bool isMC; //define here to use it in "process" part

  std::unique_ptr<Hists> h_pass, h_all;
};

TstarTstarTriggerSFModule::TstarTstarTriggerSFModule(uhh2::Context& ctx){

  isMC = (ctx.get("dataset_type") == "MC");

  h_trigger_decision = ctx.get_handle<bool>("trigger_decision");
  h_trigger_decision_ele = ctx.get_handle<bool>("trigger_decision_ele");

  if(debug) cout << "Declare Output ... " << endl;
  ctx.undeclare_all_event_output();
  h_eta = ctx.declare_event_output<double>("eta");
  h_pt = ctx.declare_event_output<double>("pt");
  h_weight = ctx.declare_event_output<double>("weight");
  h_passed = ctx.declare_event_output<bool>("passed");
  h_run = ctx.declare_event_output<int>("run");
  h_lumi = ctx.declare_event_output<int>("lumi");
  h_eventnr = ctx.declare_event_output<int>("eventnr");

  h_pass.reset(new ElectronHists(ctx, "pass_Elec"));
  h_all.reset(new ElectronHists(ctx, "all_Elec"));

}


bool TstarTstarTriggerSFModule::process(uhh2::Event& event){
  if(debug) cout << " ------------------------------------------- " << endl;
  if(debug) cout << " ----------------- NewEvent ---------------- " << endl;
  if(debug) cout << " ------------------------------------------- " << endl;

  if (!event.get(h_trigger_decision)) return false;

  // fill hists with all events
  h_all->fill(event);

  // HERE FILL PT AND ETA HISTS FOR PASSING AND NOT PASSING ELEC TRIGGER
  if(debug) cout << "Start Fill ... " << endl;
  

  if(debug) cout << "Set after Trigger pass... " << endl;
  event.set(h_pt, event.electrons->at(0).pt());
  event.set(h_eta, event.electrons->at(0).eta());
  event.set(h_weight, event.weight);
  event.set(h_passed, event.get(h_trigger_decision_ele));

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

  if( event.get(h_trigger_decision_ele) ){
    // fill pass histograms
    h_pass->fill(event);
  }

  // if(debug) cout << "Event done ... " << endl;
  if(debug) cout << "Event done ... " << endl;
  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarTriggerSFModule)
