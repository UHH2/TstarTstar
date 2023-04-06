#include <iostream>
#include <memory>

// TODO clean. Do i need all those?
// UHH2 stuff
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/PhotonIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/MCWeight.h"
#include <UHH2/common/include/DetectorCleaning.h>
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/LeptonScaleFactors.h"
#include "UHH2/common/include/TopPtReweight.h"
#include "UHH2/common/include/MCWeight.h"

// TstarTstar stuff
#include "UHH2/TstarTstar/include/ModuleBASE.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/TstarTstar/include/TstarTstarScaleFactors.h"
#include "UHH2/TstarTstar/include/ElecTriggerSF.h"
#include "UHH2/TstarTstar/include/TstarTstarJetCorrectionHists.h"
#include "UHH2/TstarTstar/include/JetMETCorrections.h"

// HOTVR stuff
#include "UHH2/HOTVR/include/HOTVRIds.h"
#include "UHH2/HOTVR/include/HadronicTop.h"
#include "UHH2/HOTVR/include/HOTVRScaleFactor.h"
#include "UHH2/HOTVR/include/HOTVRJetCorrector.h"
#include "UHH2/HOTVR/include/HOTVRJetCorrectionModule.h"

// namespace definitions
using namespace std;
using namespace uhh2;

namespace uhh2 {

// quick method to calculate inv_mass
float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

/** \brief Module for the full selection in T*T*->ttbar gg search
 *
 *
 */
class TstarTstarSelectionModule: public AnalysisModule {
public:

    explicit TstarTstarSelectionModule(Context & ctx);
    virtual bool process(Event & event) override;

private:

  // ###### Modules ######
  // general
  std::unique_ptr<JetMETCorrections> jetcorrections;
  std::unique_ptr<AnalysisModule> AK4cleaner;
  std::unique_ptr<AnalysisModule> reco_primlep;
  std::unique_ptr<AnalysisModule> ttgenprod;

  // HOTVR-related
  std::unique_ptr<AnalysisModule> HOTVRCorr;
  std::unique_ptr<HOTVRJetLeptonCleaner> HOTVRjlc;
  std::unique_ptr<TopJetCleaner> HOTVRcleaner;

  // scale factors
  std::unique_ptr<AnalysisModule> ScaleFactor_btagging;
  std::unique_ptr<AnalysisModule> ScaleFactor_NNLO;
  std::unique_ptr<AnalysisModule> MCScaleVariations;
  std::unique_ptr<AnalysisModule> TopPtReweighting;

  // lepton SFs
  // muon
  std::unique_ptr<AnalysisModule> sf_muon_iso;
  std::unique_ptr<AnalysisModule> sf_muon_ID_lowpt;
  std::unique_ptr<AnalysisModule> sf_muon_ID_highpt;
  std::unique_ptr<AnalysisModule> sf_muon_trigger_lowpt;
  std::unique_ptr<AnalysisModule> sf_muon_trigger_highpt;
  std::unique_ptr<AnalysisModule> sf_muon_iso_DUMMY;
  std::unique_ptr<AnalysisModule> sf_muon_ID_DUMMY;
  std::unique_ptr<AnalysisModule> sf_muon_trigger_DUMMY;

  // electron
  std::unique_ptr<AnalysisModule> sf_ele_ID_lowpt;
  std::unique_ptr<AnalysisModule> sf_ele_ID_highpt;
  std::unique_ptr<AnalysisModule> sf_ele_reco;
  std::unique_ptr<AnalysisModule> sf_ele_ID_DUMMY;
  std::unique_ptr<AnalysisModule> sf_ele_reco_DUMMY;
  std::unique_ptr<AnalysisModule> sf_ele_trigger; // no dummy needed as only applied when 1 electron present

  // selections
  unique_ptr<Selection> met_sel;
  unique_ptr<Selection> twodcut_sel;
  unique_ptr<HEMCleanerSelection> HEMCleaner;
  unique_ptr<HEMCleanerMCScale> HEMCleanerMCScaler;

  // ###### Histograms ######

  // set of nominal ones
  std::unique_ptr<Hists> h_beginSel,                 h_prefiring,                 h_jetcorrections,            h_jetcleaning,            h_triggerSF,             h_leptonSF,             h_HEMcut,             h_METcut,             h_AK4cut;
  std::unique_ptr<Hists> h_HOTVRcut,                 h_bcorrections,              h_byield,                 h_btagcut,              h_2Dcut,              h_STcut,              h_theorycorrections;

  // set of muon ones
  std::unique_ptr<Hists> h_beginSel_mu,              h_prefiring_mu,              h_jetcorrections_mu,         h_jetcleaning_mu,         h_triggerSF_mu,           h_leptonSF_mu,           h_HEMcut_mu,           h_METcut_mu,           h_AK4cut_mu;
  std::unique_ptr<Hists> h_HOTVRcut_mu,              h_bcorrections_mu,           h_byield_mu,              h_btagcut_mu,            h_2Dcut_mu,            h_STcut_mu,            h_theorycorrections_mu;
  std::unique_ptr<Hists> h_beginSel_mu_lowpt,        h_prefiring_mu_lowpt,        h_jetcorrections_mu_lowpt,   h_jetcleaning_mu_lowpt,   h_triggerSF_mu_lowpt,     h_leptonSF_mu_lowpt,     h_HEMcut_mu_lowpt,     h_METcut_mu_lowpt,     h_AK4cut_mu_lowpt;
  std::unique_ptr<Hists> h_HOTVRcut_mu_lowpt,        h_bcorrections_mu_lowpt,     h_byield_mu_lowpt,        h_btagcut_mu_lowpt,      h_2Dcut_mu_lowpt,      h_STcut_mu_lowpt,      h_theorycorrections_mu_lowpt;
  std::unique_ptr<Hists> h_beginSel_mu_highpt,       h_prefiring_mu_highpt,       h_jetcorrections_mu_highpt,  h_jetcleaning_mu_highpt,  h_triggerSF_mu_highpt,    h_leptonSF_mu_highpt,    h_HEMcut_mu_highpt,    h_METcut_mu_highpt,    h_AK4cut_mu_highpt;
  std::unique_ptr<Hists> h_HOTVRcut_mu_highpt,       h_bcorrections_mu_highpt,    h_byield_mu_highpt,       h_btagcut_mu_highpt,     h_2Dcut_mu_highpt,     h_STcut_mu_highpt,     h_theorycorrections_mu_highpt;

  // set of electron ones
  std::unique_ptr<Hists> h_beginSel_ele,             h_prefiring_ele,             h_jetcorrections_ele,        h_jetcleaning_ele,        h_triggerSF_ele,          h_leptonSF_ele,          h_HEMcut_ele,          h_METcut_ele,          h_AK4cut_ele;
  std::unique_ptr<Hists> h_HOTVRcut_ele,             h_bcorrections_ele,          h_byield_ele,             h_btagcut_ele,           h_2Dcut_ele,           h_STcut_ele,           h_theorycorrections_ele;
  std::unique_ptr<Hists> h_beginSel_ele_lowpt,       h_prefiring_ele_lowpt,       h_jetcorrections_ele_lowpt,  h_jetcleaning_ele_lowpt,  h_triggerSF_ele_lowpt,    h_leptonSF_ele_lowpt,    h_HEMcut_ele_lowpt,    h_METcut_ele_lowpt,    h_AK4cut_ele_lowpt;
  std::unique_ptr<Hists> h_HOTVRcut_ele_lowpt,       h_bcorrections_ele_lowpt,    h_byield_ele_lowpt,       h_btagcut_ele_lowpt,     h_2Dcut_ele_lowpt,     h_STcut_ele_lowpt,     h_theorycorrections_ele_lowpt;
  std::unique_ptr<Hists> h_beginSel_ele_highpt,      h_prefiring_ele_highpt,      h_jetcorrections_ele_highpt, h_jetcleaning_ele_highpt, h_triggerSF_ele_highpt,   h_leptonSF_ele_highpt,   h_HEMcut_ele_highpt,   h_METcut_ele_highpt,   h_AK4cut_ele_highpt;
  std::unique_ptr<Hists> h_HOTVRcut_ele_highpt,      h_bcorrections_ele_highpt,   h_byield_ele_highpt,      h_btagcut_ele_highpt,    h_2Dcut_ele_highpt,    h_STcut_ele_highpt,    h_theorycorrections_ele_highpt;

  // no-btag CR histograms for all steps after btagcut
  std::unique_ptr<Hists> h_nob_btagcut,              h_nob_2Dcut,                 h_nob_STcut,              h_nob_theorycorrections;
  std::unique_ptr<Hists> h_nob_btagcut_mu,           h_nob_2Dcut_mu,              h_nob_STcut_mu,           h_nob_theorycorrections_mu;
  std::unique_ptr<Hists> h_nob_btagcut_mu_lowpt,     h_nob_2Dcut_mu_lowpt,        h_nob_STcut_mu_lowpt,     h_nob_theorycorrections_mu_lowpt;
  std::unique_ptr<Hists> h_nob_btagcut_mu_highpt,    h_nob_2Dcut_mu_highpt,       h_nob_STcut_mu_highpt,    h_nob_theorycorrections_mu_highpt;
  std::unique_ptr<Hists> h_nob_btagcut_ele,          h_nob_2Dcut_ele,             h_nob_STcut_ele,          h_nob_theorycorrections_ele;
  std::unique_ptr<Hists> h_nob_btagcut_ele_lowpt,    h_nob_2Dcut_ele_lowpt,       h_nob_STcut_ele_lowpt,    h_nob_theorycorrections_ele_lowpt;
  std::unique_ptr<Hists> h_nob_btagcut_ele_highpt,   h_nob_2Dcut_ele_highpt,      h_nob_STcut_ele_highpt,   h_nob_theorycorrections_ele_highpt;

  // ttbar control region
  std::unique_ptr<Hists> h_theorycorrections_ttbar;

  // set of jet correction hists
  std::unique_ptr<Hists> h_jetCorr_beginSel,          h_jetCorr_prefiring,        h_jetCorr_jetcorrections,          h_jetCorr_jetcleaning,         h_jetCorr_theorycorrections;
  std::unique_ptr<Hists> h_jetCorr_beginSel_mu,       h_jetCorr_prefiring_mu,     h_jetCorr_jetcorrections_mu,       h_jetCorr_jetcleaning_mu,      h_jetCorr_theorycorrections_mu;
  std::unique_ptr<Hists> h_jetCorr_beginSel_ele,      h_jetCorr_prefiring_ele,    h_jetCorr_jetcorrections_ele,      h_jetCorr_jetcleaning_ele,     h_jetCorr_theorycorrections_ele;

  // ###### Handles ######
  uhh2::Event::Handle<bool> h_trigger_decision;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<bool> h_is_muevt;
  uhh2::Event::Handle<bool> h_is_highpt;
  uhh2::Event::Handle<double> h_evt_weight;
  uhh2::Event::Handle<double> h_ST_AK4;
  uhh2::Event::Handle<double> h_ST_HOTVR;
  uhh2::Event::Handle<bool> h_is_btagevent;
  uhh2::Event::Handle<bool> h_is_ttbarCR;
  uhh2::Event::Handle<bool> h_MC_isfake2017B;
  uhh2::Event::Handle<bool> h_MC_isfake2016B;


  // ###### other needed definitions ######
  bool debug = false;
  bool is_MC;
  bool data_isMu = false;
  bool data_isEG = false;
  bool data_isEle = false;
  bool data_isPho = false;
  bool isTriggerSFMeasurement = false;
  TString year;
  TString Prefiring_direction;

  // event yield factors, needed for b-tagging SFs.
  TH2D *eventYieldFactors;

};


TstarTstarSelectionModule::TstarTstarSelectionModule(Context & ctx) {

  // setting debug from xml file
  if(ctx.get("debug", "<not set>") == "true") debug = true;

  // debug message
  if(debug) {
    cout << "Hello World from TstarTstarSelectionModule!" << endl;
    for(auto & kv : ctx.get_all()){cout << " " << kv.first << " = " << kv.second << endl;}
  }

  // ###### 0. Setting variables ######
  // MC or real data
  is_MC = ctx.get("dataset_type") == "MC";

  // getting SF file path
  string SF_path = ctx.get("SF_path");

  // checking if this is trigger file measurement
  string IsTriggerSFMeasurement = ctx.get("IsTriggerSFMeasurement", "False");
  if(IsTriggerSFMeasurement == "True") isTriggerSFMeasurement = true;

  // year of samples
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

  // ###### 1. Set up modules ######

  // common modules was replaced here
  // I'll use the jet correction module from Christopher
  // this will do the following
  // - Jet corrections
  // - Jet resolution smearing
  // - PF ID Jet cleaning
  // - Jet Lepton Cleaning
  jetcorrections.reset(new JetMETCorrections());
  jetcorrections->switch_met_xy_correction(false); // I do not want to apply these
  jetcorrections->init(ctx);  
  
  // - JetCleaner
  const JetId jetID = PtEtaCut(30, 2.5); // only this needed, PFID is done within jet correction module
  AK4cleaner.reset(new JetCleaner(ctx, jetID));

  // correcting and cleaning HOTVR
  HOTVRjlc.reset(new HOTVRJetLeptonCleaner());
  HOTVRCorr.reset(new HOTVRJetCorrectionModule(ctx));
  HOTVRcleaner.reset(new TopJetCleaner(ctx, AndId<Jet>(PtEtaCut(200, 2.5), JetPFID(JetPFID::WP_TIGHT_PUPPI)) ));

  // primary lepton
  reco_primlep.reset(new PrimaryLepton(ctx));

  // lepton scale factors
  // electrons
  if(debug) cout << "Setting up electron scale." << endl;
  sf_ele_reco.reset( new uhh2::ElectronRecoScaleFactors(ctx) );
  sf_ele_reco_DUMMY.reset(new uhh2::ElectronRecoScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true) );
  sf_ele_ID_lowpt.reset(new uhh2::ElectronIdScaleFactors(ctx, Electron::mvaEleID_Fall17_iso_V2_wp90) );
  sf_ele_ID_highpt.reset( new uhh2::ElectronIdScaleFactors(ctx, Electron::mvaEleID_Fall17_noIso_V2_wp90) );
  sf_ele_ID_DUMMY.reset( new uhh2::ElectronIdScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true) );
  if(!isTriggerSFMeasurement) sf_ele_trigger.reset( new uhh2::ElecTriggerSF(ctx, "central", "eta_ptbins", year) );

  // muons
  if(debug) cout << "Setting up muon scale." << endl;
  sf_muon_ID_lowpt.reset( new uhh2::MuonIdScaleFactors(ctx, Muon::CutBasedIdTight) );
  sf_muon_ID_highpt.reset( new uhh2::MuonIdScaleFactors(ctx, Muon::CutBasedIdGlobalHighPt) );
  sf_muon_ID_DUMMY.reset( new uhh2::MuonIdScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true) );
  sf_muon_iso.reset( new uhh2::MuonIsoScaleFactors(ctx, Muon::PFIsoTight, Muon::CutBasedIdTight) ); // only for low pt
  sf_muon_iso_DUMMY.reset( new uhh2::MuonIsoScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true) );
  sf_muon_trigger_lowpt.reset( new uhh2::MuonTriggerScaleFactors(ctx, false) );
  sf_muon_trigger_highpt.reset( new uhh2::MuonTriggerScaleFactors(ctx, true, false) );
  sf_muon_trigger_DUMMY.reset( new uhh2::MuonTriggerScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true) );

  // b-tagging scale factors
  if(debug) cout << "Setting up btagging scale." << endl;
  ScaleFactor_btagging.reset(new MCBTagDiscriminantReweighting(ctx, BTag::algo::DEEPJET)); // should be enough like this
  if(is_MC) { // TODO put this into a module at some point
    TFile *f = new TFile("/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/btagYieldSFs_"+year+".root");
    TString sample_string = "";
    if(ctx.get("dataset_version").find("TT") != std::string::npos) sample_string = "TTbar";
    else if(ctx.get("dataset_version").find("ST") != std::string::npos) sample_string = "ST";
    else if(ctx.get("dataset_version").find("WJets") != std::string::npos) sample_string = "WJets";
    else if(ctx.get("dataset_version").find("QCD") != std::string::npos) sample_string = "QCD";
    else if(ctx.get("dataset_version").find("Diboson") != std::string::npos) sample_string = "VV";
    else if(ctx.get("dataset_version").find("DY") != std::string::npos) sample_string = "DYJets";
    else if(ctx.get("dataset_version").find("Tstar") != std::string::npos) sample_string = "TstarTstar";
    if(debug) std::cout << "Apply 2D b-taggin yield SFs for " << sample_string << std::endl;

    if(sample_string != "") eventYieldFactors = (TH2D*)f->Get(sample_string);
    else throw std::runtime_error("Error: can not determine sample type for btagging yield SFs.");
  }

  // ttbar on GEN
  if(is_MC) ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));

  // corrections
  ScaleFactor_NNLO.reset(new NLOCorrections(ctx));
  TopPtReweighting.reset( new TopPtReweight(ctx, 0.0615, -0.0005, "ttbargen", "weight_ttbar", true) );

  Prefiring_direction = ctx.get("Sys_prefiring", "nominal");
  MCScaleVariations.reset(new MCScaleVariation(ctx) );

  // ###### 2. set up selections ######
  if(debug) cout << "Setting up Selections." << endl;

  // MET selection
  met_sel.reset(new METCut  (50.,1e9));

  // 2D cut
  twodcut_sel.reset(new TwoDCut(0.4, 25.0));  // The same as in Z'->ttbar semileptonic

  // HEM cut
  HEMCleaner.reset(new HEMCleanerSelection(ctx, "jets", "topjets"));
  HEMCleanerMCScaler.reset(new HEMCleanerMCScale(ctx, "jets", "topjets"));


  // ###### 3. Set up Hists ######
  if(debug) cout << "Setting up Hists." << endl;

  h_beginSel.reset(new TstarTstarHists(ctx, "beginSel"));
  h_beginSel_mu.reset(new TstarTstarHists(ctx, "beginSel_mu"));
  h_beginSel_mu_lowpt.reset(new TstarTstarHists(ctx, "beginSel_mu_lowpt"));
  h_beginSel_mu_highpt.reset(new TstarTstarHists(ctx, "beginSel_mu_highpt"));
  h_beginSel_ele.reset(new TstarTstarHists(ctx, "beginSel_ele"));
  h_beginSel_ele_lowpt.reset(new TstarTstarHists(ctx, "beginSel_ele_lowpt"));
  h_beginSel_ele_highpt.reset(new TstarTstarHists(ctx, "beginSel_ele_highpt"));

  h_prefiring.reset(new TstarTstarHists(ctx, "prefiring"));
  h_prefiring_mu.reset(new TstarTstarHists(ctx, "prefiring_mu"));
  h_prefiring_mu_lowpt.reset(new TstarTstarHists(ctx, "prefiring_mu_lowpt"));
  h_prefiring_mu_highpt.reset(new TstarTstarHists(ctx, "prefiring_mu_highpt"));
  h_prefiring_ele.reset(new TstarTstarHists(ctx, "prefiring_ele"));
  h_prefiring_ele_lowpt.reset(new TstarTstarHists(ctx, "prefiring_ele_lowpt"));
  h_prefiring_ele_highpt.reset(new TstarTstarHists(ctx, "prefiring_ele_highpt"));

  h_jetcorrections.reset(new TstarTstarHists(ctx, "jetcorrections"));
  h_jetcorrections_mu.reset(new TstarTstarHists(ctx, "jetcorrections_mu"));
  h_jetcorrections_mu_lowpt.reset(new TstarTstarHists(ctx, "jetcorrections_mu_lowpt"));
  h_jetcorrections_mu_highpt.reset(new TstarTstarHists(ctx, "jetcorrections_mu_highpt"));
  h_jetcorrections_ele.reset(new TstarTstarHists(ctx, "jetcorrections_ele"));
  h_jetcorrections_ele_lowpt.reset(new TstarTstarHists(ctx, "jetcorrections_ele_lowpt"));
  h_jetcorrections_ele_highpt.reset(new TstarTstarHists(ctx, "jetcorrections_ele_highpt"));

  h_jetcleaning.reset(new TstarTstarHists(ctx, "jetcleaning"));
  h_jetcleaning_mu.reset(new TstarTstarHists(ctx, "jetcleaning_mu"));
  h_jetcleaning_mu_lowpt.reset(new TstarTstarHists(ctx, "jetcleaning_mu_lowpt"));
  h_jetcleaning_mu_highpt.reset(new TstarTstarHists(ctx, "jetcleaning_mu_highpt"));
  h_jetcleaning_ele.reset(new TstarTstarHists(ctx, "jetcleaning_ele"));
  h_jetcleaning_ele_lowpt.reset(new TstarTstarHists(ctx, "jetcleaning_ele_lowpt"));
  h_jetcleaning_ele_highpt.reset(new TstarTstarHists(ctx, "jetcleaning_ele_highpt"));

  h_triggerSF.reset(new TstarTstarHists(ctx, "triggerSF"));
  h_triggerSF_mu.reset(new TstarTstarHists(ctx, "triggerSF_mu"));
  h_triggerSF_mu_lowpt.reset(new TstarTstarHists(ctx, "triggerSF_mu_lowpt"));
  h_triggerSF_mu_highpt.reset(new TstarTstarHists(ctx, "triggerSF_mu_highpt"));
  h_triggerSF_ele.reset(new TstarTstarHists(ctx, "triggerSF_ele"));
  h_triggerSF_ele_lowpt.reset(new TstarTstarHists(ctx, "triggerSF_ele_lowpt"));
  h_triggerSF_ele_highpt.reset(new TstarTstarHists(ctx, "triggerSF_ele_highpt"));

  h_leptonSF.reset(new TstarTstarHists(ctx, "leptonSF"));
  h_leptonSF_mu.reset(new TstarTstarHists(ctx, "leptonSF_mu"));
  h_leptonSF_mu_lowpt.reset(new TstarTstarHists(ctx, "leptonSF_mu_lowpt"));
  h_leptonSF_mu_highpt.reset(new TstarTstarHists(ctx, "leptonSF_mu_highpt"));
  h_leptonSF_ele.reset(new TstarTstarHists(ctx, "leptonSF_ele"));
  h_leptonSF_ele_lowpt.reset(new TstarTstarHists(ctx, "leptonSF_ele_lowpt"));
  h_leptonSF_ele_highpt.reset(new TstarTstarHists(ctx, "leptonSF_ele_highpt"));

  h_HEMcut.reset(new TstarTstarHists(ctx, "HEMcut"));
  h_HEMcut_mu.reset(new TstarTstarHists(ctx, "HEMcut_mu"));
  h_HEMcut_mu_lowpt.reset(new TstarTstarHists(ctx, "HEMcut_mu_lowpt"));
  h_HEMcut_mu_highpt.reset(new TstarTstarHists(ctx, "HEMcut_mu_highpt"));
  h_HEMcut_ele.reset(new TstarTstarHists(ctx, "HEMcut_ele"));
  h_HEMcut_ele_lowpt.reset(new TstarTstarHists(ctx, "HEMcut_ele_lowpt"));
  h_HEMcut_ele_highpt.reset(new TstarTstarHists(ctx, "HEMcut_ele_highpt"));

  h_METcut.reset(new TstarTstarHists(ctx, "METcut"));
  h_METcut_mu.reset(new TstarTstarHists(ctx, "METcut_mu"));
  h_METcut_mu_lowpt.reset(new TstarTstarHists(ctx, "METcut_mu_lowpt"));
  h_METcut_mu_highpt.reset(new TstarTstarHists(ctx, "METcut_mu_highpt"));
  h_METcut_ele.reset(new TstarTstarHists(ctx, "METcut_ele"));
  h_METcut_ele_lowpt.reset(new TstarTstarHists(ctx, "METcut_ele_lowpt"));
  h_METcut_ele_highpt.reset(new TstarTstarHists(ctx, "METcut_ele_highpt"));

  h_AK4cut.reset(new TstarTstarHists(ctx, "AK4cut"));
  h_AK4cut_mu.reset(new TstarTstarHists(ctx, "AK4cut_mu"));
  h_AK4cut_mu_lowpt.reset(new TstarTstarHists(ctx, "AK4cut_mu_lowpt"));
  h_AK4cut_mu_highpt.reset(new TstarTstarHists(ctx, "AK4cut_mu_highpt"));
  h_AK4cut_ele.reset(new TstarTstarHists(ctx, "AK4cut_ele"));
  h_AK4cut_ele_lowpt.reset(new TstarTstarHists(ctx, "AK4cut_ele_lowpt"));
  h_AK4cut_ele_highpt.reset(new TstarTstarHists(ctx, "AK4cut_ele_highpt"));

  h_HOTVRcut.reset(new TstarTstarHists(ctx, "HOTVRcut"));
  h_HOTVRcut_mu.reset(new TstarTstarHists(ctx, "HOTVRcut_mu"));
  h_HOTVRcut_mu_lowpt.reset(new TstarTstarHists(ctx, "HOTVRcut_mu_lowpt"));
  h_HOTVRcut_mu_highpt.reset(new TstarTstarHists(ctx, "HOTVRcut_mu_highpt"));
  h_HOTVRcut_ele.reset(new TstarTstarHists(ctx, "HOTVRcut_ele"));
  h_HOTVRcut_ele_lowpt.reset(new TstarTstarHists(ctx, "HOTVRcut_ele_lowpt"));
  h_HOTVRcut_ele_highpt.reset(new TstarTstarHists(ctx, "HOTVRcut_ele_highpt"));

  h_bcorrections.reset(new TstarTstarHists(ctx, "bcorrections"));
  h_bcorrections_mu.reset(new TstarTstarHists(ctx, "bcorrections_mu"));
  h_bcorrections_mu_lowpt.reset(new TstarTstarHists(ctx, "bcorrections_mu_lowpt"));
  h_bcorrections_mu_highpt.reset(new TstarTstarHists(ctx, "bcorrections_mu_highpt"));
  h_bcorrections_ele.reset(new TstarTstarHists(ctx, "bcorrections_ele"));
  h_bcorrections_ele_lowpt.reset(new TstarTstarHists(ctx, "bcorrections_ele_lowpt"));
  h_bcorrections_ele_highpt.reset(new TstarTstarHists(ctx, "bcorrections_ele_highpt"));

  h_byield.reset(new TstarTstarHists(ctx, "byield"));
  h_byield_mu.reset(new TstarTstarHists(ctx, "byield_mu"));
  h_byield_mu_lowpt.reset(new TstarTstarHists(ctx, "byield_mu_lowpt"));
  h_byield_mu_highpt.reset(new TstarTstarHists(ctx, "byield_mu_highpt"));
  h_byield_ele.reset(new TstarTstarHists(ctx, "byield_ele"));
  h_byield_ele_lowpt.reset(new TstarTstarHists(ctx, "byield_ele_lowpt"));
  h_byield_ele_highpt.reset(new TstarTstarHists(ctx, "byield_ele_highpt"));

  h_btagcut.reset(new TstarTstarHists(ctx, "btagcut"));
  h_btagcut_mu.reset(new TstarTstarHists(ctx, "btagcut_mu"));
  h_btagcut_mu_lowpt.reset(new TstarTstarHists(ctx, "btagcut_mu_lowpt"));
  h_btagcut_mu_highpt.reset(new TstarTstarHists(ctx, "btagcut_mu_highpt"));
  h_btagcut_ele.reset(new TstarTstarHists(ctx, "btagcut_ele"));
  h_btagcut_ele_lowpt.reset(new TstarTstarHists(ctx, "btagcut_ele_lowpt"));
  h_btagcut_ele_highpt.reset(new TstarTstarHists(ctx, "btagcut_ele_highpt"));

  h_2Dcut.reset(new TstarTstarHists(ctx, "2Dcut"));
  h_2Dcut_mu.reset(new TstarTstarHists(ctx, "2Dcut_mu"));
  h_2Dcut_mu_lowpt.reset(new TstarTstarHists(ctx, "2Dcut_mu_lowpt"));
  h_2Dcut_mu_highpt.reset(new TstarTstarHists(ctx, "2Dcut_mu_highpt"));
  h_2Dcut_ele.reset(new TstarTstarHists(ctx, "2Dcut_ele"));
  h_2Dcut_ele_lowpt.reset(new TstarTstarHists(ctx, "2Dcut_ele_lowpt"));
  h_2Dcut_ele_highpt.reset(new TstarTstarHists(ctx, "2Dcut_ele_highpt"));

  h_STcut.reset(new TstarTstarHists(ctx, "STcut"));
  h_STcut_mu.reset(new TstarTstarHists(ctx, "STcut_mu"));
  h_STcut_mu_lowpt.reset(new TstarTstarHists(ctx, "STcut_mu_lowpt"));
  h_STcut_mu_highpt.reset(new TstarTstarHists(ctx, "STcut_mu_highpt"));
  h_STcut_ele.reset(new TstarTstarHists(ctx, "STcut_ele"));
  h_STcut_ele_lowpt.reset(new TstarTstarHists(ctx, "STcut_ele_lowpt"));
  h_STcut_ele_highpt.reset(new TstarTstarHists(ctx, "STcut_ele_highpt"));

  h_theorycorrections.reset(new TstarTstarHists(ctx, "theorycorrections"));
  h_theorycorrections_mu.reset(new TstarTstarHists(ctx, "theorycorrections_mu"));
  h_theorycorrections_mu_lowpt.reset(new TstarTstarHists(ctx, "theorycorrections_mu_lowpt"));
  h_theorycorrections_mu_highpt.reset(new TstarTstarHists(ctx, "theorycorrections_mu_highpt"));
  h_theorycorrections_ele.reset(new TstarTstarHists(ctx, "theorycorrections_ele"));
  h_theorycorrections_ele_lowpt.reset(new TstarTstarHists(ctx, "theorycorrections_ele_lowpt"));
  h_theorycorrections_ele_highpt.reset(new TstarTstarHists(ctx, "theorycorrections_ele_highpt"));
  h_theorycorrections_ttbar.reset(new TstarTstarHists(ctx, "theorycorrections_ttbar"));
  
  h_nob_btagcut.reset(new TstarTstarHists(ctx, "nob_btagcut"));
  h_nob_btagcut_mu.reset(new TstarTstarHists(ctx, "nob_btagcut_mu"));
  h_nob_btagcut_mu_lowpt.reset(new TstarTstarHists(ctx, "nob_btagcut_lowpt"));
  h_nob_btagcut_mu_highpt.reset(new TstarTstarHists(ctx, "nob_btagcut_mu_highpt"));
  h_nob_btagcut_ele.reset(new TstarTstarHists(ctx, "nob_btagcut_ele"));
  h_nob_btagcut_ele_lowpt.reset(new TstarTstarHists(ctx, "nob_btagcut_ele_lowpt"));
  h_nob_btagcut_ele_highpt.reset(new TstarTstarHists(ctx, "nob_btagcut_ele_highpt"));

  h_nob_2Dcut.reset(new TstarTstarHists(ctx, "nob_2Dcut"));
  h_nob_2Dcut_mu.reset(new TstarTstarHists(ctx, "nob_2Dcut_mu"));
  h_nob_2Dcut_mu_lowpt.reset(new TstarTstarHists(ctx, "nob_2Dcut_lowpt"));
  h_nob_2Dcut_mu_highpt.reset(new TstarTstarHists(ctx, "nob_2Dcut_mu_highpt"));
  h_nob_2Dcut_ele.reset(new TstarTstarHists(ctx, "nob_2Dcut_ele"));
  h_nob_2Dcut_ele_lowpt.reset(new TstarTstarHists(ctx, "nob_2Dcut_ele_lowpt"));
  h_nob_2Dcut_ele_highpt.reset(new TstarTstarHists(ctx, "nob_2Dcut_ele_highpt"));

  h_nob_STcut.reset(new TstarTstarHists(ctx, "nob_STcut"));
  h_nob_STcut_mu.reset(new TstarTstarHists(ctx, "nob_STcut_mu"));
  h_nob_STcut_mu_lowpt.reset(new TstarTstarHists(ctx, "nob_STcut_lowpt"));
  h_nob_STcut_mu_highpt.reset(new TstarTstarHists(ctx, "nob_STcut_mu_highpt"));
  h_nob_STcut_ele.reset(new TstarTstarHists(ctx, "nob_STcut_ele"));
  h_nob_STcut_ele_lowpt.reset(new TstarTstarHists(ctx, "nob_STcut_ele_lowpt"));
  h_nob_STcut_ele_highpt.reset(new TstarTstarHists(ctx, "nob_STcut_ele_highpt"));

  h_nob_theorycorrections.reset(new TstarTstarHists(ctx, "nob_theorycorrections"));
  h_nob_theorycorrections_mu.reset(new TstarTstarHists(ctx, "nob_theorycorrections_mu"));
  h_nob_theorycorrections_mu_lowpt.reset(new TstarTstarHists(ctx, "nob_theorycorrections_lowpt"));
  h_nob_theorycorrections_mu_highpt.reset(new TstarTstarHists(ctx, "nob_theorycorrections_mu_highpt"));
  h_nob_theorycorrections_ele.reset(new TstarTstarHists(ctx, "nob_theorycorrections_ele"));
  h_nob_theorycorrections_ele_lowpt.reset(new TstarTstarHists(ctx, "nob_theorycorrections_ele_lowpt"));
  h_nob_theorycorrections_ele_highpt.reset(new TstarTstarHists(ctx, "nob_theorycorrections_ele_highpt"));

  // jet correction hists
  h_jetCorr_beginSel.reset(new TstarTstarJetCorrectionHists(ctx, "jetCorr_beginSel"));
  h_jetCorr_prefiring.reset(new TstarTstarJetCorrectionHists(ctx, "jetCorr_prefiring"));
  h_jetCorr_jetcorrections.reset(new TstarTstarJetCorrectionHists(ctx, "jetCorr_jetcorrections"));
  h_jetCorr_jetcleaning.reset(new TstarTstarJetCorrectionHists(ctx, "jetCorr_jetcleaning"));
  h_jetCorr_theorycorrections.reset(new TstarTstarJetCorrectionHists(ctx, "jetCorr_theorycorrections"));

  h_jetCorr_beginSel_mu.reset(new TstarTstarJetCorrectionHists(ctx, "jetCorr_beginSel_mu"));
  h_jetCorr_prefiring_mu.reset(new TstarTstarJetCorrectionHists(ctx, "jetCorr_prefiring_mu"));
  h_jetCorr_jetcorrections_mu.reset(new TstarTstarJetCorrectionHists(ctx, "jetCorr_jetcorrections_mu"));
  h_jetCorr_jetcleaning_mu.reset(new TstarTstarJetCorrectionHists(ctx, "jetCorr_jetcleaning_mu"));
  h_jetCorr_theorycorrections_mu.reset(new TstarTstarJetCorrectionHists(ctx, "jetCorr_theorycorrections_mu"));

  h_jetCorr_beginSel_ele.reset(new TstarTstarJetCorrectionHists(ctx, "jetCorr_beginSel_ele"));
  h_jetCorr_prefiring_ele.reset(new TstarTstarJetCorrectionHists(ctx, "jetCorr_prefiring_ele"));
  h_jetCorr_jetcorrections_ele.reset(new TstarTstarJetCorrectionHists(ctx, "jetCorr_jetcorrections_ele"));
  h_jetCorr_jetcleaning_ele.reset(new TstarTstarJetCorrectionHists(ctx, "jetCorr_jetcleaning_ele"));
  h_jetCorr_theorycorrections_ele.reset(new TstarTstarJetCorrectionHists(ctx, "jetCorr_theorycorrections_ele"));

  // ###### 4. Init Handles ######
  h_trigger_decision = ctx.get_handle<bool>("trigger_decision");
  h_is_muevt = ctx.get_handle<bool>("is_muevt");
  h_is_highpt = ctx.get_handle<bool>("is_highpt");
  h_evt_weight = ctx.get_handle<double>("evt_weight");
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_MC_isfake2017B = ctx.get_handle<bool>("MC_isfake2017B");
  h_MC_isfake2016B = ctx.get_handle<bool>("MC_isfake2016B");

  h_ST_AK4 = ctx.declare_event_output<double>("ST_AK4");
  h_ST_HOTVR = ctx.declare_event_output<double>("ST_HOTVR");
  h_is_btagevent = ctx.declare_event_output<bool>("is_btagevent");
  h_is_ttbarCR = ctx.declare_event_output<bool>("is_ttbarCR");

}


bool TstarTstarSelectionModule::process(Event & event) {

  // debug status message
  if(debug) cout << "TstarTstarSelectionModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;

  if(debug && event.get(h_MC_isfake2017B)) cout << "This MC event was randomly assigned to 2017 Run B handling" << endl;
  if(debug && event.get(h_is_highpt)) cout << "The event is classified as having high pt" << std::endl;

  // reapply event weights from handle
  event.weight = event.get(h_evt_weight);
  if(debug) cout << "weights applied." << endl;

  // Fill ttgen object for correct matching check, etc
  if(is_MC) ttgenprod->process(event);
  if(debug) std::cout << "Filled ttgenprod" << endl;

  // set primary lepton
  reco_primlep->process(event);
  if(debug) std::cout << "Got primary lepton" << endl;

  // hists before anything happened
  if(debug) std::cout << "Fill Crosscheck hists" << endl;
  if(event.get(h_trigger_decision)) {
    h_beginSel->fill(event);
    if(is_MC) h_jetCorr_beginSel->fill(event);
    if(event.get(h_is_muevt)){
      h_beginSel_mu->fill(event);
      if(is_MC) h_jetCorr_beginSel_mu->fill(event);
      if(event.get(h_is_highpt)) h_beginSel_mu_highpt->fill(event);
      else h_beginSel_mu_lowpt->fill(event);
    } else {
      h_beginSel_ele->fill(event);
      if(is_MC) h_jetCorr_beginSel_ele->fill(event);
      if(event.get(h_is_highpt)) h_beginSel_ele_highpt->fill(event);
      else h_beginSel_ele_lowpt->fill(event);
    }
  }

  // Prefiring weights
  if (is_MC) {
     if (Prefiring_direction == "nominal") event.weight *= event.prefiringWeight;
     else if (Prefiring_direction == "up") event.weight *= event.prefiringWeightUp;
     else if (Prefiring_direction == "down") event.weight *= event.prefiringWeightDown;
  }
  // writing MC weights
  MCScaleVariations->process(event);

  // hists before anything happened
  if(debug) std::cout << "Fill Crosscheck hists" << endl;
  if(event.get(h_trigger_decision)) {
    h_prefiring->fill(event);
    if(is_MC) h_jetCorr_prefiring->fill(event);
    if(event.get(h_is_muevt)){
      h_prefiring_mu->fill(event);
      if(is_MC) h_jetCorr_prefiring_mu->fill(event);
      if(event.get(h_is_highpt)) h_prefiring_mu_highpt->fill(event);
      else h_prefiring_mu_lowpt->fill(event);
    } else {
      h_prefiring_ele->fill(event);
      if(is_MC) h_jetCorr_prefiring_ele->fill(event);
      if(event.get(h_is_highpt)) h_prefiring_ele_highpt->fill(event);
      else h_prefiring_ele_lowpt->fill(event);
    }
  }
  
  // JetMET corrections (and jet lepton cleaning)
  if(!jetcorrections->process(event)) return false;
  if(!(HOTVRjlc->process(event))) return false;
  if(!(HOTVRCorr->process(event))) return false;

  // hists after jet corrections
  if(debug) std::cout << "Fill correction hists" << endl;
  if(event.get(h_trigger_decision)) {
    h_jetcorrections->fill(event);
    if(is_MC) h_jetCorr_jetcorrections->fill(event);
    if(event.get(h_is_muevt)){
      h_jetcorrections_mu->fill(event);
      if(is_MC) h_jetCorr_jetcorrections_mu->fill(event);
      if(event.get(h_is_highpt)) h_jetcorrections_mu_highpt->fill(event);
      else h_jetcorrections_mu_lowpt->fill(event);
    } else {
      h_jetcorrections_ele->fill(event);
      if(is_MC) h_jetCorr_jetcorrections_ele->fill(event);
      if(event.get(h_is_highpt)) h_jetcorrections_ele_highpt->fill(event);
      else h_jetcorrections_ele_lowpt->fill(event);
    }
  }

  // Jet cleaning
  if(!(AK4cleaner->process(event))) return false;
  if(!(HOTVRcleaner->process(event))) return false;

  // st calculation, for usage in plotting later on
  double st = 0.;
  for(const auto & lepton : *event.electrons) st += lepton.pt();
  for(const auto & lepton : *event.muons) st += lepton.pt();
  st += event.met->pt();
  double stHOTVR = st;

  for(const auto & jet : *event.jets) st += jet.pt();
  for(const auto & jet : *event.topjets) stHOTVR += jet.pt();

  event.set(h_ST_AK4, st);
  event.set(h_ST_HOTVR, stHOTVR);

  // hists after jet cleaning
  if(debug) std::cout << "Fill jet cleaning hists" << endl;
  if(event.get(h_trigger_decision)) {
    h_jetcleaning->fill(event);
    if(is_MC) h_jetCorr_jetcleaning->fill(event);
    if(event.get(h_is_muevt)){
      h_jetcleaning_mu->fill(event);
      if(is_MC) h_jetCorr_jetcleaning_mu->fill(event);
      if(event.get(h_is_highpt)) h_jetcleaning_mu_highpt->fill(event);
      else h_jetcleaning_mu_lowpt->fill(event);
    } else {
      h_jetcleaning_ele->fill(event);
      if(is_MC) h_jetCorr_jetcleaning_ele->fill(event);
      if(event.get(h_is_highpt)) h_jetcleaning_ele_highpt->fill(event);
      else h_jetcleaning_ele_lowpt->fill(event);
    }
  }

  // trigger SFs
  // these are only applied if the event is triggered!
  if(debug) std::cout << "before muon trigger SFs" << std::endl;
  if(is_MC && event.get(h_is_muevt) && !isTriggerSFMeasurement && event.get(h_trigger_decision)){
    if(event.get(h_is_highpt)) sf_muon_trigger_highpt->process(event);
    else sf_muon_trigger_lowpt->process(event);
  } else {
    sf_muon_trigger_DUMMY->process(event);
  }
  if(debug) std::cout << "before electron trigger SFs" << std::endl;
  if(!isTriggerSFMeasurement) sf_ele_trigger->process(event); // this one is its own dummy!

  // hists after triggerSFs
  if(debug) std::cout << "Fill triggerSF hists" << endl;
  if(event.get(h_trigger_decision)) {
    h_triggerSF->fill(event);
    if(event.get(h_is_muevt)){
      h_triggerSF_mu->fill(event);
      if(event.get(h_is_highpt)) h_triggerSF_mu_highpt->fill(event);
      else h_triggerSF_mu_lowpt->fill(event);
    } else {
      h_triggerSF_ele->fill(event);
      if(event.get(h_is_highpt)) h_triggerSF_ele_highpt->fill(event);
      else h_triggerSF_ele_lowpt->fill(event);
    }
  }

  // lepton SFs
  if(event.get(h_is_muevt)){

    if(event.get(h_is_highpt)) {
      sf_muon_ID_highpt->process(event);
      sf_muon_iso_DUMMY->process(event);
    } else {
      sf_muon_ID_lowpt->process(event);
      sf_muon_iso->process(event);
    }

    if(isTriggerSFMeasurement) {
      sf_ele_reco->process(event);

      if(event.electrons->at(0).pt() > 120) {
        sf_ele_ID_highpt->process(event);
      } else {
        sf_ele_ID_lowpt->process(event);
      }
    } else {
      sf_ele_reco_DUMMY->process(event);
      sf_ele_ID_DUMMY->process(event);
    }

  } else {

    sf_ele_reco->process(event);

    if(event.get(h_is_highpt)) {
      sf_ele_ID_highpt->process(event);
    } else {
      sf_ele_ID_lowpt->process(event);
    }

    sf_muon_ID_DUMMY->process(event);
    sf_muon_iso_DUMMY->process(event);

  }
  if(debug) std::cout << "Done Lepton ID, ISO SFs" << endl;

  // hists after lepton SFs
  if(debug) std::cout << "Fill lepton SF hists" << endl;
  if(event.get(h_trigger_decision)) {
    h_leptonSF->fill(event);
    if(event.get(h_is_muevt)){
      h_leptonSF_mu->fill(event);
      if(event.get(h_is_highpt)) h_leptonSF_mu_highpt->fill(event);
      else h_leptonSF_mu_lowpt->fill(event);
    } else {
      h_leptonSF_ele->fill(event);
      if(event.get(h_is_highpt)) h_leptonSF_ele_highpt->fill(event);
      else h_leptonSF_ele_lowpt->fill(event);
    }
  }

  // exclude electron in overlap region
  if(!event.get(h_is_muevt) || isTriggerSFMeasurement) {
    if( (abs(event.electrons->at(0).eta()) >= 1.444) && (abs(event.electrons->at(0).eta()) <= 1.566)) return false;
  }

  // Fixing the HEM issue (only 2018)
  if(!HEMCleaner->passes(event)) return false;
  if(is_MC && (year == "2018" || year == "UL18")) HEMCleanerMCScaler->process(event);

  // hists after HEM & lepton overlap removal
  if(debug) std::cout << "Fill HEM hists" << endl;
  if(event.get(h_trigger_decision)) {
    h_HEMcut->fill(event);
    if(event.get(h_is_muevt)){
      h_HEMcut_mu->fill(event);
      if(event.get(h_is_highpt)) h_HEMcut_mu_highpt->fill(event);
      else h_HEMcut_mu_lowpt->fill(event);
    } else {
      h_HEMcut_ele->fill(event);
      if(event.get(h_is_highpt)) h_HEMcut_ele_highpt->fill(event);
      else h_HEMcut_ele_lowpt->fill(event);
    }
  }


  // MET selection
  bool pass_MET = met_sel->passes(event);
  if(!pass_MET) return false;

  // hists after MET cut
  if(debug) std::cout << "Fill MET cut hists" << endl;
  if(event.get(h_trigger_decision)) {
    h_METcut->fill(event);
    if(event.get(h_is_muevt)){
      h_METcut_mu->fill(event);
      if(event.get(h_is_highpt)) h_METcut_mu_highpt->fill(event);
      else h_METcut_mu_lowpt->fill(event);
    } else {
      h_METcut_ele->fill(event);
      if(event.get(h_is_highpt)) h_METcut_ele_highpt->fill(event);
      else h_METcut_ele_lowpt->fill(event);
    }
  }


  // AK4 jet selection
  bool pass_njet = (event.jets->size()>3);
  if(isTriggerSFMeasurement) pass_njet = (event.jets->size()>1); // only need to require 2 jets for trigger SF measurement
  if(!pass_njet) return false;

  // hists after AK4 cut
  if(debug) std::cout << "Fill AK4 cut hists" << endl;
  if(event.get(h_trigger_decision)) {
    h_AK4cut->fill(event);
    if(event.get(h_is_muevt)){
      h_AK4cut_mu->fill(event);
      if(event.get(h_is_highpt)) h_AK4cut_mu_highpt->fill(event);
      else h_AK4cut_mu_lowpt->fill(event);
    } else {
      h_AK4cut_ele->fill(event);
      if(event.get(h_is_highpt)) h_AK4cut_ele_highpt->fill(event);
      else h_AK4cut_ele_lowpt->fill(event);
    }
  }


  // ###### HOTVR jet selection ######
  bool pass_fat_njet = (event.topjets->size()>0);
  if(!pass_fat_njet && !isTriggerSFMeasurement) return false;

  // hists after HOTVR cut
  if(debug) std::cout << "Fill HOTVR cut hists" << endl;
  if(event.get(h_trigger_decision)) {
    h_HOTVRcut->fill(event);
    if(event.get(h_is_muevt)){
      h_HOTVRcut_mu->fill(event);
      if(event.get(h_is_highpt)) h_HOTVRcut_mu_highpt->fill(event);
      else h_HOTVRcut_mu_lowpt->fill(event);
    } else {
      h_HOTVRcut_ele->fill(event);
      if(event.get(h_is_highpt)) h_HOTVRcut_ele_highpt->fill(event);
      else h_HOTVRcut_ele_lowpt->fill(event);
    }
  }


  // b-tagging sfs
  ScaleFactor_btagging->process(event);

  // hists after b-tagging SFs
  if(debug) std::cout << "Fill btagging SF hists" << endl;
  if(event.get(h_trigger_decision)) {
    h_bcorrections->fill(event);
    if(event.get(h_is_muevt)){
      h_bcorrections_mu->fill(event);
      if(event.get(h_is_highpt)) h_bcorrections_mu_highpt->fill(event);
      else h_bcorrections_mu_lowpt->fill(event);
    } else {
      h_bcorrections_ele->fill(event);
      if(event.get(h_is_highpt)) h_bcorrections_ele_highpt->fill(event);
      else h_bcorrections_ele_lowpt->fill(event);
    }
  }

  // b-tagging yield correction
  // done as a function of AK4 HT and N(AK4)
  if(is_MC) {
    double ht = 0.;
    for(const auto & jet : *event.jets) ht += jet.pt();
    if(ht >= 4000.) ht = 3999.9;

    double btaggingYieldWeight = eventYieldFactors->GetBinContent( eventYieldFactors->GetXaxis()->FindBin(ht),  eventYieldFactors->GetYaxis()->FindBin(event.jets->size()) );
    event.weight *= btaggingYieldWeight;
  }

  // defining a ttbar CR, now that btags have been corrected
  TopJetId topjetID = AndId<TopJet>(HOTVRTopTag(), Tau32Groomed(0.56)); // Top Tag that is used later
  bool passHOTVRtoptag = false;
  for (const auto & jet: *event.topjets){
    if(topjetID(jet, event)) passHOTVRtoptag = true;
  }
  BTag bJetID_loose = BTag(BTag::algo::DEEPJET, BTag::wp::WP_LOOSE);
  int N_btag_loose = 0;
  for (const auto & jet: *event.jets){
    if(bJetID_loose(jet, event)) N_btag_loose++;
  }
  if(passHOTVRtoptag && N_btag_loose > 1) event.set(h_is_ttbarCR, true);
  else event.set(h_is_ttbarCR, false);

  // hists after b-tagging yield corrections
  if(debug) std::cout << "Fill btagging yield hists" << endl;
  if(event.get(h_trigger_decision)) {
    h_byield->fill(event);
    if(event.get(h_is_muevt)){
      h_byield_mu->fill(event);
      if(event.get(h_is_highpt)) h_byield_mu_highpt->fill(event);
      else h_byield_mu_lowpt->fill(event);
    } else {
      h_byield_ele->fill(event);
      if(event.get(h_is_highpt)) h_byield_ele_highpt->fill(event);
      else h_byield_ele_lowpt->fill(event);
    }
  }


  // btag selection
  BTag bJetID = BTag(BTag::algo::DEEPJET, BTag::wp::WP_MEDIUM);
  bool pass_btagcut = false;
  for (const auto & jet: *event.jets){
    if(bJetID(jet, event)) pass_btagcut = true;
  }
  event.set(h_is_btagevent, pass_btagcut); // not throwing events away, as we'll keep these for iur

  // hists after b-tagging cut
  if(debug) std::cout << "Fill btagging cut hists" << endl;
  if(event.get(h_trigger_decision)) {
    if(event.get(h_is_btagevent)) {
      h_btagcut->fill(event);
      if(event.get(h_is_muevt)){
        h_btagcut_mu->fill(event);
        if(event.get(h_is_highpt)) h_btagcut_mu_highpt->fill(event);
        else h_btagcut_mu_lowpt->fill(event);
      } else {
        h_btagcut_ele->fill(event);
        if(event.get(h_is_highpt)) h_btagcut_ele_highpt->fill(event);
        else h_btagcut_ele_lowpt->fill(event);
      }
    } else {
      h_nob_btagcut->fill(event);
      if(event.get(h_is_muevt)){
        h_nob_btagcut_mu->fill(event);
        if(event.get(h_is_highpt)) h_nob_btagcut_mu_highpt->fill(event);
        else h_nob_btagcut_mu_lowpt->fill(event);
      } else {
        h_nob_btagcut_ele->fill(event);
        if(event.get(h_is_highpt)) h_nob_btagcut_ele_highpt->fill(event);
        else h_nob_btagcut_ele_lowpt->fill(event);
      }
    }
  }


  // lepton-jet 2D cut
  if(debug) std::cout << "Start 2D cut" << endl;
  bool pass_2D = true;
  for(auto& muo : *event.muons){
    if(debug) cout<<"AFTER Muon (pt,eta): "<<muo.pt()<<", "<<muo.eta()<<endl;
    float    dRmin, pTrel;
    std::tie(dRmin, pTrel) = drmin_pTrel(muo, *event.jets);
    muo.set_tag(Muon::twodcut_dRmin, dRmin);
    muo.set_tag(Muon::twodcut_pTrel, pTrel);
  }
  for(auto& ele : *event.electrons){
    if(debug) cout<<"Electron (pt,eta): "<<ele.pt()<<", "<<ele.eta()<<endl;
    float    dRmin, pTrel;
    std::tie(dRmin, pTrel) = drmin_pTrel(ele, *event.jets);
    ele.set_tag(Electron::twodcut_dRmin, dRmin);
    ele.set_tag(Electron::twodcut_pTrel, pTrel);
  }
  const bool pass_twodcut = twodcut_sel->passes(event);
  if(event.get(h_is_highpt)) pass_2D = pass_twodcut;
  if(!pass_2D) return false;

  // hists after 2D cut
  if(debug) std::cout << "Fill 2D cut hists" << endl;
  if(event.get(h_trigger_decision)) {
    if(event.get(h_is_btagevent)) {
      h_2Dcut->fill(event);
      if(event.get(h_is_muevt)){
        h_2Dcut_mu->fill(event);
        if(event.get(h_is_highpt)) h_2Dcut_mu_highpt->fill(event);
        else h_2Dcut_mu_lowpt->fill(event);
      } else {
        h_2Dcut_ele->fill(event);
        if(event.get(h_is_highpt)) h_2Dcut_ele_highpt->fill(event);
        else h_2Dcut_ele_lowpt->fill(event);
      }
    } else {
      h_nob_2Dcut->fill(event);
      if(event.get(h_is_muevt)){
        h_nob_2Dcut_mu->fill(event);
        if(event.get(h_is_highpt)) h_nob_2Dcut_mu_highpt->fill(event);
        else h_nob_2Dcut_mu_lowpt->fill(event);
      } else {
        h_nob_2Dcut_ele->fill(event);
        if(event.get(h_is_highpt)) h_nob_2Dcut_ele_highpt->fill(event);
        else h_nob_2Dcut_ele_lowpt->fill(event);
      }
    }
  }


  // st cut
  if(stHOTVR < 500 && !isTriggerSFMeasurement) return false;

  // hists after ST cut
  if(debug) std::cout << "Fill ST cut hists" << endl;
  if(event.get(h_trigger_decision)) {
    if(event.get(h_is_btagevent)) {
      h_STcut->fill(event);
      if(event.get(h_is_muevt)){
        h_STcut_mu->fill(event);
        if(event.get(h_is_highpt)) h_STcut_mu_highpt->fill(event);
        else h_STcut_mu_lowpt->fill(event);
      } else {
        h_STcut_ele->fill(event);
        if(event.get(h_is_highpt)) h_STcut_ele_highpt->fill(event);
        else h_STcut_ele_lowpt->fill(event);
      }
    } else {
      h_nob_STcut->fill(event);
      if(event.get(h_is_muevt)){
        h_nob_STcut_mu->fill(event);
        if(event.get(h_is_highpt)) h_nob_STcut_mu_highpt->fill(event);
        else h_nob_STcut_mu_lowpt->fill(event);
      } else {
        h_nob_STcut_ele->fill(event);
        if(event.get(h_is_highpt)) h_nob_STcut_ele_highpt->fill(event);
        else h_nob_STcut_ele_lowpt->fill(event);
      }
    }
  }


  // NNLO corrections
  ScaleFactor_NNLO->process(event);
  TopPtReweighting->process(event); // is not applied, but calculated

  // hists after theory corrections
  if(debug) std::cout << "Fill theory correction hists" << endl;
  if(event.get(h_trigger_decision)) {
    if(event.get(h_is_btagevent)) {
      h_theorycorrections->fill(event);
      if(is_MC) h_jetCorr_theorycorrections->fill(event);
      if(event.get(h_is_ttbarCR)) {h_theorycorrections_ttbar->fill(event);}
      if(event.get(h_is_muevt)){
        h_theorycorrections_mu->fill(event);
        if(is_MC) h_jetCorr_theorycorrections_mu->fill(event);
        if(event.get(h_is_highpt)) h_theorycorrections_mu_highpt->fill(event);
        else h_theorycorrections_mu_lowpt->fill(event);
      } else {
        h_theorycorrections_ele->fill(event);
        if(is_MC) h_jetCorr_theorycorrections_ele->fill(event);
        if(event.get(h_is_highpt)) h_theorycorrections_ele_highpt->fill(event);
        else h_theorycorrections_ele_lowpt->fill(event);
      }
    } else {
      h_nob_theorycorrections->fill(event);
      if(event.get(h_is_muevt)){
        h_nob_theorycorrections_mu->fill(event);
        if(event.get(h_is_highpt)) h_nob_theorycorrections_mu_highpt->fill(event);
        else h_nob_theorycorrections_mu_lowpt->fill(event);
      } else {
        h_nob_theorycorrections_ele->fill(event);
        if(event.get(h_is_highpt)) h_nob_theorycorrections_ele_highpt->fill(event);
        else h_nob_theorycorrections_ele_lowpt->fill(event);
      }
    }
  }

  // save final weight for outputting later
  event.set(h_evt_weight, event.weight);
  return true;

}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarSelectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarSelectionModule)

}
