#include <iostream>
#include <memory>
#include <string>

// UHH2 stuff
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/PhotonIds.h"
#include <UHH2/common/include/MuonIds.h>
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/LumiSelection.h"


// TstarTstar custom stuff
#include "UHH2/TstarTstar/include/TstarTstarCustomIds.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/TstarTstar/include/TstarTstarSFHists.h"
#include "UHH2/TstarTstar/include/TstarTstarScaleFactors.h"
#include "UHH2/TstarTstar/include/TstarTstarElectronIDHists.h"
#include "UHH2/TstarTstar/include/TstarTstarPDFNormHists.h"

// namespace definition
using namespace std;
using namespace uhh2;

namespace uhh2 {

/** \brief Module for the T*T*->ttbar gg/gamma preselection
 *
 * Corrects all objects via CommonModules and applies some loose cuts
 *
 */
class TstarTstarPreselectionModule: public AnalysisModule {
public:

  explicit TstarTstarPreselectionModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  // ##### Modules #####
  std::unique_ptr<Selection> lumi_selection;
  std::unique_ptr<AndSelection> metfilters_selection;
  std::vector<std::unique_ptr<AnalysisModule>> modules;
  std::unique_ptr<AnalysisModule> MCScaleVariations;

  // GEN stuff (used for top pt reweighting)
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;

  // trigger selections
  // these can be filled with different things per year, see below
  std::unique_ptr<Selection> trg_ele_low;
  std::unique_ptr<Selection> trg_ele_high;
  std::unique_ptr<Selection> trg_pho;

  std::unique_ptr<Selection> trg_mu_low_1;
  std::unique_ptr<Selection> trg_mu_low_2;
  std::unique_ptr<Selection> trg_mu_high_1;
  std::unique_ptr<Selection> trg_mu_high_2;
  std::unique_ptr<Selection> trg_mu_high_3;

  // ##### Histograms #####
  // total hists
  std::unique_ptr<Hists> h_nocuts,     h_common,         h_lepsel,        h_jetsel,       h_ST;

  // electron channel
  std::unique_ptr<Hists> h_lepsel_ele,          h_jetsel_ele,           h_ST_ele;
  std::unique_ptr<Hists> h_lepsel_ele_lowpt,    h_jetsel_ele_lowpt,     h_ST_ele_lowpt;
  std::unique_ptr<Hists> h_lepsel_ele_highpt,   h_jetsel_ele_highpt,    h_ST_ele_highpt;

  // muon channel
  std::unique_ptr<Hists> h_lepsel_mu,           h_jetsel_mu,            h_ST_mu;
  std::unique_ptr<Hists> h_lepsel_mu_lowpt,     h_jetsel_mu_lowpt,      h_ST_mu_lowpt;
  std::unique_ptr<Hists> h_lepsel_mu_highpt,    h_jetsel_mu_highpt,     h_ST_mu_highpt;

  // luminosity histograms
  std::unique_ptr<LuminosityHists> lumihist_common, lumihist_lepsel, lumihist_jetsel;

  // GEN histograms
  std::unique_ptr<Hists> h_nocuts_gen, h_common_gen,     h_lepsel_gen,    h_jetsel_gen,   h_ST_gen;
  std::unique_ptr<Hists> h_afterSelection_gen, h_afterSelection_genmatch;

  // histograms for electron ID
  std::unique_ptr<Hists> h_electronIDhists;

  // hist PDF NORM
  std::unique_ptr<Hists> h_PDFnorm;

  // ##### Handles #####
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<bool> h_is_muevt;
  uhh2::Event::Handle<bool> h_is_highpt;
  uhh2::Event::Handle<double> h_evt_weight;
  uhh2::Event::Handle<bool> h_trigger_decision;
  uhh2::Event::Handle<bool> h_trigger_decision_ele;
  uhh2::Event::Handle<bool> h_MC_isfake2017B;
  uhh2::Event::Handle<bool> h_MC_isfake2016B;

  // ##### other needed definitions #####
  // these are mainly booleans set in the constructor and used in the process method
  TString year;
  bool debug = false;
  bool is_MC;
  bool data_isMu = false;
  bool data_isEG = false;
  bool data_isEle = false;
  bool data_isPho = false;
  bool data_is2017B = false;
  bool data_is2016B = false;
  bool isTriggerSFMeasurement = false;


};


TstarTstarPreselectionModule::TstarTstarPreselectionModule(Context & ctx){

  // setting debugging flag from xml file
  if(ctx.get("debug", "<not set>") == "true") debug = true;

  // debugging status message
  if(debug) {
    cout << "Hello World from TstarTstarPreselectionModule!" << endl;
    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }
  }

  // ###### 0. Fetching required information ######
  // year of sample
  year = ctx.get("year", "<not set>");
  if(year == "<not set>"){
    if(ctx.get("dataset_version").find("2016") != std::string::npos) year = "2016";
    else if(ctx.get("dataset_version").find("2017") != std::string::npos) year = "2017";
    else if(ctx.get("dataset_version").find("2018") != std::string::npos) year = "2018";
    else if(ctx.get("dataset_version").find("UL16preVFP") != std::string::npos) year = "UL16preVFP";
    else if(ctx.get("dataset_version").find("UL16postVFP") != std::string::npos) year = "UL16postVFP";
    else if(ctx.get("dataset_version").find("UL17") != std::string::npos) year = "UL17";
    else if(ctx.get("dataset_version").find("UL18") != std::string::npos) year = "UL18";
    else throw "No year found in dataset name!";
  }
  if(debug) cout << "Year is " << year << "." << endl;

  // check if we are MC or real data
  is_MC = ctx.get("dataset_type") == "MC";

  // if its a data sample, check what kind
  if(!is_MC) {
    data_isMu = (ctx.get("dataset_version").find("SingleMuon") != std::string::npos);
    data_isEle = (ctx.get("dataset_version").find("SingleElectron") != std::string::npos);
    data_isPho = (ctx.get("dataset_version").find("SinglePhoton") != std::string::npos);
    data_isEG = (ctx.get("dataset_version").find("EGamma") != std::string::npos);

    if(data_isMu) std::cout << "This data sample is a muon sample." << std::endl;
    if(data_isPho) std::cout << "This data sample is a photon sample" << std::endl;
    if(data_isEle) std::cout << "This data sample is an electron sample." << std::endl;
    if(data_isEG) std::cout << "This data sample is an egamma sample" << std::endl;
  }

  // check if this is specific run
  // needed as some triggers were missing
  if(!is_MC) data_is2016B = ctx.get("dataset_version").find("RunB_UL16preVFP") != std::string::npos;
  if(!is_MC) data_is2017B = ctx.get("dataset_version").find("RunB_UL17") != std::string::npos;

  // ###### 1. setting up modules ######
  if(debug) cout << "Setting up modules" << endl;

  // TRIGGER
  // Muon channel
  std::cout << "Setting up muon triggers" << std::endl;

  if(year == "2016" || year == "UL16preVFP" || year == "UL16postVFP") {
    trg_mu_low_1.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    trg_mu_low_2.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
    trg_mu_high_1.reset(new TriggerSelection("HLT_Mu50_v*"));
    if(!data_is2016B) trg_mu_high_2.reset(new TriggerSelection("HLT_TkMu50_v*"));
  }
  else if(year == "2017" || year == "UL17") {
    trg_mu_low_1.reset(new TriggerSelection("HLT_IsoMu27_v*"));
    trg_mu_high_1.reset(new TriggerSelection("HLT_Mu50_v*"));
    trg_mu_high_2.reset(new TriggerSelection("HLT_TkMu100_v*"));
    trg_mu_high_3.reset(new TriggerSelection("HLT_OldMu100_v*"));
  }
  else if(year == "2018" || year == "UL18") {
    trg_mu_low_1.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    trg_mu_high_1.reset(new TriggerSelection("HLT_Mu50_v*"));
    trg_mu_high_2.reset(new TriggerSelection("HLT_TkMu100_v*"));
    trg_mu_high_3.reset(new TriggerSelection("HLT_OldMu100_v*"));
  }

  // Electron channel
  std::cout << "Setting up ele triggers" << std::endl;

  if(year == "2016" || year == "UL16preVFP" || year == "UL16postVFP") {
    trg_ele_low.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));
    trg_ele_high.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
    trg_pho.reset(new TriggerSelection("HLT_Photon175_v*"));
  }
  else if(year == "2017" || year == "UL17") {
    trg_ele_low.reset(new TriggerSelection("HLT_Ele35_WPTight_Gsf_v*"));
    trg_ele_high.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
    trg_pho.reset(new TriggerSelection("HLT_Photon200_v*"));
  }
  else if(year == "2018" || year == "UL18") {
    trg_ele_low.reset(new TriggerSelection("HLT_Ele32_WPTight_Gsf_v*"));
    trg_ele_high.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
    trg_pho.reset(new TriggerSelection("HLT_Photon200_v*"));
  }

  // CommonModules
  // replacing common modules (which was originally here) by manual calls of the appropriate things
  // what should be done here:
  // - good run selection (for data only) based on lumi_file defined in xml input
  lumi_selection.reset(new LumiSelection(ctx));
  // - apply MET filters (this is done early in common modules as well)
  metfilters_selection.reset(new AndSelection(ctx, "metfilters"));
  metfilters_selection->add<TriggerSelection>("goodVertices"                      ,"Flag_goodVertices");
  metfilters_selection->add<TriggerSelection>("globalSuperTightHalo2016Filter"    ,"Flag_globalSuperTightHalo2016Filter");
  metfilters_selection->add<TriggerSelection>("HBHENoiseFilter"                   ,"Flag_HBHENoiseFilter");
  metfilters_selection->add<TriggerSelection>("HBHENoiseIsoFilter"                ,"Flag_HBHENoiseIsoFilter");
  metfilters_selection->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter","Flag_EcalDeadCellTriggerPrimitiveFilter");
  metfilters_selection->add<TriggerSelection>("BadPFMuonFilter"                   ,"Flag_BadPFMuonFilter");
  metfilters_selection->add<TriggerSelection>("BadPFMuonDzFilter"                 ,"Flag_BadPFMuonDzFilter");
  metfilters_selection->add<TriggerSelection>("eeBadScFilter"                     ,"Flag_eeBadScFilter");
  if ( year == "UL17" || year == "UL18" ) {
    metfilters_selection->add<TriggerSelection>("Flag_ecalBadCalibFilter", "Flag_ecalBadCalibFilter");
  }
  PrimaryVertexId pvid=StandardPrimaryVertexId();
  metfilters_selection->add<NPVSelection>("1 good PV",1,-1,pvid);
  // - Primary vertex cleaner (remove all non-good PVs from the list of primary vertices)
  modules.emplace_back(new PrimaryVertexCleaner(pvid));
  if(is_MC){
    // - MCLumiWeight (for MC only; only has an effect if "use_sframe_weight" is set to false)
    modules.emplace_back(new MCLumiWeight(ctx));
    // - MCPileupReweight (for MC only)
    modules.emplace_back(new MCPileupReweight(ctx, "central"));
  } 

  // this will set the mc scale variation weights
  MCScaleVariations.reset(new MCScaleVariation(ctx) );


  // - Electron Cleaner
  ElectronId eleID_lowpt = ElectronTagID(Electron::mvaEleID_Fall17_iso_V2_wp90);
  ElectronId eleID_highpt = ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wp90);
  ElectronId total_eleID = OrId<Electron>( AndId<Electron>(PtEtaCut(40., 2.4), eleID_lowpt, EleMaxPtCut(120.)),  AndId<Electron>(PtEtaCut(120., 2.4), eleID_highpt));
  modules.emplace_back(new ElectronCleaner(total_eleID));

  // - Muon Cleaner
  MuonId muID_lowpt = AndId<Muon>(MuonID(Muon::CutBasedIdTight), MuonID(Muon::PFIsoTight));
  MuonId muID_highpt = MuonID(Muon::CutBasedIdGlobalHighPt);
  MuonId total_muID = OrId<Muon>( AndId<Muon>(PtEtaCut(30., 2.4), muID_lowpt, MuMaxPtCut(55.)), AndId<Muon>(PtEtaCut(55, 2.4), muID_highpt));
  modules.emplace_back(new MuonCleaner(total_muID));

  // find ttbar for GEN
  if(is_MC){
    // Prepare GEN
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
    if(debug) cout << "TTbarGenProducer done" << endl;
  }

  // ###### 3. setting up histograms ######
  // general
  h_nocuts.reset(new TstarTstarHists(ctx, "NoCuts"));
  h_common.reset(new TstarTstarHists(ctx, "AfterCommon"));
  h_lepsel.reset(new TstarTstarHists(ctx, "AfterLep"));
  h_jetsel.reset(new TstarTstarHists(ctx, "AfterJets"));
  h_ST.reset(new TstarTstarHists(ctx, "AfterST"));

  // electron channel
  h_lepsel_ele.reset(new TstarTstarHists(ctx, "AfterLepSel_ele"));
  h_jetsel_ele.reset(new TstarTstarHists(ctx, "AfterJets_ele"));
  h_ST_ele.reset(new TstarTstarHists(ctx, "AfterST_ele"));

  h_lepsel_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterLepSel_ele_lowpt"));
  h_jetsel_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterJets_ele_lowpt"));
  h_ST_ele_lowpt.reset(new TstarTstarHists(ctx, "AfterST_ele_lowpt"));

  h_lepsel_ele_highpt.reset(new TstarTstarHists(ctx, "AfterLepSel_ele_highpt"));
  h_jetsel_ele_highpt.reset(new TstarTstarHists(ctx, "AfterJets_ele_highpt"));
  h_ST_ele_highpt.reset(new TstarTstarHists(ctx, "AfterST_ele_highpt"));

  // muon channel
  h_lepsel_mu.reset(new TstarTstarHists(ctx, "AfterLepSel_mu"));
  h_jetsel_mu.reset(new TstarTstarHists(ctx, "AfterJets_mu"));
  h_ST_mu.reset(new TstarTstarHists(ctx, "AfterST_mu"));

  h_lepsel_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterLepSel_mu_lowpt"));
  h_jetsel_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterJets_mu_lowpt"));
  h_ST_mu_lowpt.reset(new TstarTstarHists(ctx, "AfterST_mu_lowpt"));

  h_lepsel_mu_highpt.reset(new TstarTstarHists(ctx, "AfterLepSel_mu_highpt"));
  h_jetsel_mu_highpt.reset(new TstarTstarHists(ctx, "AfterJets_mu_highpt"));
  h_ST_mu_highpt.reset(new TstarTstarHists(ctx, "AfterST_mu_highpt"));

  if(debug) cout << "normal hists donr done" << endl;

  // GEN hists
  h_nocuts_gen.reset(new TstarTstarGenHists(ctx, "NoCuts_gen"));
  h_common_gen.reset(new TstarTstarGenHists(ctx, "AfterCommon_gen"));
  h_lepsel_gen.reset(new TstarTstarGenHists(ctx, "AfterLep_gen"));
  h_jetsel_gen.reset(new TstarTstarGenHists(ctx, "AfterJets_gen"));
  h_ST_gen.reset(new TstarTstarGenHists(ctx, "AfterST_gen"));
  h_afterSelection_gen.reset(new TstarTstarGenHists(ctx, "AfterSel_gen"));
  h_afterSelection_genmatch.reset(new TstarTstarGenRecoMatchedHists(ctx, "AfterSel_genmatch"));

  if(debug) cout << "GEN hists donr done" << endl;

  // Lumi hists
  lumihist_common.reset(new LuminosityHists(ctx, "lumihist_AfterCommon"));
  lumihist_lepsel.reset(new LuminosityHists(ctx, "lumihist_AfterLep"));
  lumihist_jetsel.reset(new LuminosityHists(ctx, "lumihist_AfterJets"));

  if(debug) cout << "lumi hists donr done" << endl;

  // other histograms  
  h_electronIDhists.reset(new TstarTstarElectronIDHists(ctx, "ElectronIDHists"));
  if(is_MC) h_PDFnorm.reset(new TstarTstarPDFNormHists(ctx, "PDFNorm"));

  // ###### 4. init handles ######
  h_is_muevt = ctx.declare_event_output<bool>("is_muevt");
  h_is_highpt = ctx.declare_event_output<bool>("is_highpt");
  h_evt_weight = ctx.declare_event_output<double>("evt_weight");
  h_MC_isfake2017B = ctx.declare_event_output<bool>("MC_isfake2017B");
  h_MC_isfake2016B = ctx.declare_event_output<bool>("MC_isfake2016B");
  h_trigger_decision = ctx.declare_event_output<bool>("trigger_decision");
  h_trigger_decision_ele = ctx.declare_event_output<bool>("trigger_decision_ele");

  // ###### 4. other MISC stuff ######
  // handle trigger SF measurement, which changes some things
  string IsTriggerSFMeasurement = ctx.get("IsTriggerSFMeasurement", "False");
  if(IsTriggerSFMeasurement == "True") isTriggerSFMeasurement = true;
  if(isTriggerSFMeasurement) std::cout << "Changing the lepton requirement as this is Electron trigger SF measurement run!" << std::endl;

  // init srand
  // used later to set random MC "runs"
  std::srand(std::time(nullptr));

}


bool TstarTstarPreselectionModule::process(Event & event) {

  // debug status messages
  if(debug){
    cout << "TstarTstarPreselectionModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
    cout<<"N muons = "<<event.muons->size()<<", N electrons = "<<event.electrons->size()<<", N photons = "<<event.photons->size()<<endl;
    for(auto& muo : *event.muons) cout<<"BEFORE Muon (pt,eta): "<<muo.pt()<<", "<<muo.eta()<<endl;
  }

  // Fill ttbar object for gen
  if(is_MC){
    ttgenprod->process(event);
  }

  // hists before anything happened
  h_nocuts->fill(event);
  h_nocuts_gen->fill(event);
  if(is_MC) h_electronIDhists->fill(event);
  if(debug) cout<<"Filled hists without any cuts, and GEN with no cuts"<<endl;

  // ###### common modules, corrections & cleaning ######
  if(event.isRealData) if(!lumi_selection->passes(event)) return false;
  if(!metfilters_selection->passes(event)) return false;
  for(auto & m : modules){
    m->process(event);
  }
  if(debug) cout<<"common modules done"<<endl;

  // hists before selection
  h_common->fill(event);
  h_common_gen->fill(event);
  lumihist_common->fill(event);
  MCScaleVariations->process(event); // writing MC weights
  if(is_MC) h_PDFnorm->fill(event);
  if(debug) cout<<"Filled hists after cleaning"<<endl;

  // setting fake MC "runs"
  bool MC_isfake2017B = false;
  bool MC_isfake2016B = false;

  if(year == "UL17" && is_MC) {
    // in UL17 we need to fake 11.6% (run B part) to a different trigger

    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    if(r < 0.116) MC_isfake2017B = true;

    event.set(h_MC_isfake2017B, MC_isfake2017B);
    event.set(h_MC_isfake2016B, false);
    if(debug && MC_isfake2017B) std::cout << "This MC event was set to be part of a fake 2017B run. r was " << r << "." << std::endl;

  } else if(year == "UL16preVFP" && is_MC) {
    // in UL16preVFP we need to fake 14.29% (run B part) to a different trigger

    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    if(r < 0.1429) MC_isfake2016B = true;

    event.set(h_MC_isfake2017B, false);
    event.set(h_MC_isfake2016B, MC_isfake2016B);
    if(debug && MC_isfake2016B) std::cout << "This MC event was set to be part of a fake 2017B run. r was " << r << "." << std::endl;

  } else {
    event.set(h_MC_isfake2017B, false);
    event.set(h_MC_isfake2016B, false);
  }

  // #################
  // ### Selection ###
  // #################

  // ###### Lepton selection ######
  bool pass_lep1 = (((event.muons->size() == 1) || (event.electrons->size() == 1)) && (event.electrons->size()+event.muons->size()) == 1);
  if(isTriggerSFMeasurement) pass_lep1 = ((event.muons->size() == 1) && (event.electrons->size() == 1)); // this ensures orthogonal dataset
  if(!pass_lep1) return false;

  // now that we have leptons, we can set the event channel:
  // setting muevt handle
  if(event.muons->size() == 1) event.set(h_is_muevt, true);
  else event.set(h_is_muevt, false);

  // trigger SF measurement will be in "muon mode"
  if(isTriggerSFMeasurement) event.set(h_is_muevt, true);

  // set is_highpt
  if(event.get(h_is_muevt)){
    if(event.muons->at(0).pt()<=55.) event.set(h_is_highpt, false);
    else event.set(h_is_highpt, true);
  }
  else {
    if(event.electrons->at(0).pt()<=120.) event.set(h_is_highpt, false);
    else event.set(h_is_highpt, true);
  }

  // ###### Trigger selection ######
  // important: triggers are not yet "applied" here, just calculated for plotting!
  // to be able to plot trigger efficiencies later, events that don't pass are not thrown away yet
  bool pass_trigger = false;
  bool pass_trigger_SingleMu_lowpt = false;
  bool pass_trigger_SingleMu_highpt = false;
  bool pass_trigger_SingleEle = false;

  // using muon triggers if we are a muon event
  // also, if data, only for muon samples
  if ( event.get(h_is_muevt) && (is_MC || data_isMu ) ) {
    if(debug) std::cout << "Entered muon trigger logic" << std::endl;

    if(year == "2016" || year == "UL16preVFP" || year == "UL16postVFP") {
      pass_trigger_SingleMu_lowpt = (trg_mu_low_1->passes(event) || trg_mu_low_2->passes(event));
      if(data_is2016B || event.get(h_MC_isfake2016B)) pass_trigger_SingleMu_highpt = trg_mu_high_1->passes(event);
      else pass_trigger_SingleMu_highpt = (trg_mu_high_1->passes(event) || trg_mu_high_2->passes(event));
    }
    else if(year == "2017" || year == "UL17") {
      pass_trigger_SingleMu_lowpt = trg_mu_low_1->passes(event);
      if(data_is2017B || event.get(h_MC_isfake2017B)) pass_trigger_SingleMu_highpt = trg_mu_high_1->passes(event);
      else pass_trigger_SingleMu_highpt = (trg_mu_high_1->passes(event) || trg_mu_high_2->passes(event) || trg_mu_high_3->passes(event));
    }
    else if(year == "2018" || year == "UL18") {
      pass_trigger_SingleMu_lowpt = trg_mu_low_1->passes(event);
      pass_trigger_SingleMu_highpt = (trg_mu_high_1->passes(event) || trg_mu_high_2->passes(event) || trg_mu_high_3->passes(event));
    }

    if(event.get(h_is_highpt)) pass_trigger = pass_trigger_SingleMu_highpt;
    else pass_trigger = pass_trigger_SingleMu_lowpt;

    if(debug) std::cout << "Passed muon trigger logic: " << pass_trigger << std::endl;
  }
  
  if (( !event.get(h_is_muevt) && (is_MC || data_isEG || data_isEle || data_isPho ) ) || isTriggerSFMeasurement ){

    // making sure that this is exclusive!!!
    if (( event.get(h_is_muevt) && (is_MC || data_isMu )) && !isTriggerSFMeasurement) throw std::runtime_error("Trigger selection not exclusive, and this is not a SF measurement!");

    if(debug) std::cout << "Entered electron trigger logic" << std::endl;

    if(year == "2016" || year == "UL16preVFP" || year == "UL16postVFP") {

      if(is_MC || isTriggerSFMeasurement) {
        pass_trigger_SingleEle = (trg_ele_low->passes(event) || trg_ele_high->passes(event) || trg_pho->passes(event));
      } else {
        if (data_isPho) pass_trigger_SingleEle = (!trg_ele_low->passes(event) && !trg_ele_high->passes(event) && trg_pho->passes(event));
        else pass_trigger_SingleEle = (trg_ele_low->passes(event) || trg_ele_high->passes(event));
      }

    }
    else if(year == "2017" || year == "UL17") {

      if(is_MC || isTriggerSFMeasurement) {
        if (event.get(h_MC_isfake2017B) || (isTriggerSFMeasurement && data_is2017B) ) pass_trigger_SingleEle = (trg_ele_low->passes(event) || trg_pho->passes(event));
        else pass_trigger_SingleEle = (trg_ele_low->passes(event) || trg_ele_high->passes(event) || trg_pho->passes(event));
      } else {
        if (data_is2017B) {
          if(data_isPho) pass_trigger_SingleEle = (!trg_ele_low->passes(event) && trg_pho->passes(event));
          else pass_trigger_SingleEle = trg_ele_low->passes(event);
        } else {
          if(data_isPho) pass_trigger_SingleEle = (!trg_ele_low->passes(event) && !trg_ele_high->passes(event) && trg_pho->passes(event));
          else pass_trigger_SingleEle = (trg_ele_low->passes(event) || trg_ele_high->passes(event));
        }
      }

    }
    else if(year == "2018" || year == "UL18") {
      pass_trigger_SingleEle = (trg_ele_low->passes(event) || trg_ele_high->passes(event) || trg_pho->passes(event));
    }

    if (!isTriggerSFMeasurement) pass_trigger = pass_trigger_SingleEle;

    if(debug) std::cout << "Passed electron trigger logic: " << pass_trigger << std::endl;
  }

  // setting the trigger_decision
  event.set(h_trigger_decision, pass_trigger);
  event.set(h_trigger_decision_ele, pass_trigger_SingleEle);

  // hists
  // these are after our "first real" selection step, the lepton selection
  if(pass_trigger) {
    h_lepsel->fill(event);
    h_lepsel_gen->fill(event);
    lumihist_lepsel->fill(event);
    if(event.get(h_is_muevt)){
      h_lepsel_mu->fill(event);
      if(event.get(h_is_highpt)) h_lepsel_mu_highpt->fill(event);
      else h_lepsel_mu_lowpt->fill(event);
    }
    else {
      h_lepsel_ele->fill(event);
      if(event.get(h_is_highpt)) h_lepsel_ele_highpt->fill(event);
      else h_lepsel_ele_lowpt->fill(event);
    }
    if(debug) cout << "Filled hists after lepsel" << endl;
  }

  // ###### jet selection ######
  // just a loose preselection without any pt cut on the jets
  bool pass_njet = (event.jets->size()>2);
  if(isTriggerSFMeasurement) pass_njet = (event.jets->size()>1); // only need to require 2 jets for trigger SF measurement
  if(!pass_njet) return false;

  // hists
  if(pass_trigger) {
    h_jetsel->fill(event);
    h_jetsel_gen->fill(event);
    lumihist_jetsel->fill(event);
    if(event.get(h_is_muevt)){
      h_jetsel_mu->fill(event);
      if(event.get(h_is_highpt)) h_jetsel_mu_highpt->fill(event);
      else h_jetsel_mu_lowpt->fill(event);
    }
    else {
      h_jetsel_ele->fill(event);
      if(event.get(h_is_highpt)) h_jetsel_ele_highpt->fill(event);
      else h_jetsel_ele_lowpt->fill(event);
    }
    if(debug) cout << "Filled hists after jet sel" << endl;
  }

  // ##### ST cut to suppress file size substantially! #####
  // only tmp st calculation, as we'll calculate a proper one later after the full jet sel
  double st = 0.;
  for(const auto & lepton : *event.electrons)   st += lepton.pt();
  for(const auto & lepton : *event.muons)       st += lepton.pt();
  st += event.met->pt();
  for(const auto & jet : *event.topjets)        st += jet.pt();

  // st cut
  // not applied for trigger sf measurement
  if(st < 450 && !isTriggerSFMeasurement) return false;

  // hists
  if(pass_trigger) {
    h_ST->fill(event);
    h_ST_gen->fill(event);
    lumihist_jetsel->fill(event);
    if(event.get(h_is_muevt)){
      h_ST_mu->fill(event);
      if(event.get(h_is_highpt)) h_ST_mu_highpt->fill(event);
      else h_ST_mu_lowpt->fill(event);
    }
    else {
      h_ST_ele->fill(event);
      if(event.get(h_is_highpt)) h_ST_ele_highpt->fill(event);
      else h_ST_ele_lowpt->fill(event);
    }
    if(debug) cout << "Filled hists after ST sel" << endl;
  }

  // some gen check hists
  h_afterSelection_gen->fill(event);
  h_afterSelection_genmatch->fill(event);

  // outputting event weight for following modules
  event.set(h_evt_weight, event.weight);

  if(debug) cout << "########### Done with preselection! ###########" << endl << endl;
  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarPreselectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarPreselectionModule)

}
