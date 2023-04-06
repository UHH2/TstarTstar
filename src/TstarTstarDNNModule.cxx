#include <iostream>
#include <memory>
#include <string>
#include <regex>

// UHH2 stuff
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/PhotonIds.h"
#include <UHH2/common/include/MuonIds.h>
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/MCWeight.h"

// TstarTstar
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarDNNHists.h"
#include "UHH2/TstarTstar/include/TstarTstarDNNInputHists.h"
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarAllGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/TstarTstar/include/ReconstructionTstarHypothesis.h"
#include "UHH2/TstarTstar/include/TstarTstarGenMatch.h"
#include "UHH2/TstarTstar/include/NeuralNetworkModules.h"
#include "UHH2/TstarTstar/include/TstarTstarSignalRegionHists.h"
#include "UHH2/TstarTstar/include/TstarTstarScaleFactors.h"
#include "UHH2/TstarTstar/include/TstarTstarDDTHists.h"


// other
#include "UHH2/HOTVR/include/HOTVRIds.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

/** \brief Module for the T*T*->ttbar gg MC based study
 *
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class TstarTstarDNNModule: public AnalysisModule {
public:

  explicit TstarTstarDNNModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  // ###### Modules ######
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;
  std::unique_ptr<NeuralNetworkIncluder> DNN_Includer;


  // ##### Histograms #####
  std::unique_ptr<Hists> h_topcheck, h_topcheck_reweighted, h_AfterDNNcut_06_UGLYFIX, h_NotDNNcut_06_UGLYFIX;
  std::unique_ptr<Hists> h_STreweighted, h_crosscheck, h_crosscheck_noeletriggerSF;

  std::unique_ptr<Hists> h_newTaggerSR, h_newTaggerCR, h_newTaggerCR_noelectriggerSF, h_newTagger_btagCR, h_newTagger_btagCR_DNNSR;

  std::unique_ptr<Hists> h_AfterDNNcut_02, h_AfterDNNcut_03, h_AfterDNNcut_04, h_AfterDNNcut_05, h_AfterDNNcut_06, h_AfterDNNcut_07, h_AfterDNNcut_08;
  std::unique_ptr<Hists> h_notDNNcut_02,   h_notDNNcut_03,   h_notDNNcut_04,   h_notDNNcut_05,   h_notDNNcut_06,   h_notDNNcut_07,   h_notDNNcut_08;

  std::unique_ptr<Hists> h_BtagControl_AfterDNNcut_02, h_BtagControl_AfterDNNcut_03, h_BtagControl_AfterDNNcut_04, h_BtagControl_AfterDNNcut_05, h_BtagControl_AfterDNNcut_06, h_BtagControl_AfterDNNcut_07, h_BtagControl_AfterDNNcut_08;
  std::unique_ptr<Hists> h_BtagControl_notDNNcut_02,   h_BtagControl_notDNNcut_03,   h_BtagControl_notDNNcut_04,   h_BtagControl_notDNNcut_05,   h_BtagControl_notDNNcut_06,   h_BtagControl_notDNNcut_07,   h_BtagControl_notDNNcut_08;

  std::unique_ptr<Hists> h_DNN, h_DNN_newTagger, h_DNN_newTagger_2, h_DNN_BtagControl, h_DNN_reweighted, h_DNN_reweighted_2;
  std::unique_ptr<Hists> h_DNN_lowST, h_DNN_highST, h_DNN_lowDNN, h_DNN_highDNN, h_DNN_highST_lowDNN, h_DNN_highST_highDNN;
  std::unique_ptr<Hists> h_AfterDNN, h_AfterDNN_lowST, h_AfterDNN_highST, h_AfterDNN_lowDNN, h_AfterDNN_highDNN, h_AfterDNN_highST_lowDNN, h_AfterDNN_highST_highDNN;

  std::unique_ptr<Hists> h_highLepton, h_highLepton_AfterDNNcut_06, h_highLepton_notDNNcut_06;

  std::unique_ptr<Hists> h_SignalRegion_total, h_SignalRegion_mu, h_SignalRegion_ele;
  std::unique_ptr<Hists> h_ValidationRegion_total, h_ValidationRegion_mu, h_ValidationRegion_ele;
  std::unique_ptr<Hists> h_ControlRegion_total, h_ControlRegion_mu, h_ControlRegion_ele;

  std::unique_ptr<Hists> h_SignalRegion_datadriven_FuncUp_total, h_SignalRegion_datadriven_FuncDown_total;
  std::unique_ptr<Hists> h_SignalRegion_datadriven_FuncUp_mu, h_SignalRegion_datadriven_FuncDown_mu;
  std::unique_ptr<Hists> h_SignalRegion_datadriven_FuncUp_ele, h_SignalRegion_datadriven_FuncDown_ele;
  std::unique_ptr<Hists> h_SignalRegion_datadriven_BtagUp_total, h_SignalRegion_datadriven_BtagDown_total;
  std::unique_ptr<Hists> h_SignalRegion_datadriven_BtagUp_mu, h_SignalRegion_datadriven_BtagDown_mu;
  std::unique_ptr<Hists> h_SignalRegion_datadriven_BtagUp_ele, h_SignalRegion_datadriven_BtagDown_ele;

  std::unique_ptr<Hists> h_ValidationRegion_datadriven_FuncUp_total, h_ValidationRegion_datadriven_FuncDown_total;
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_FuncUp_mu, h_ValidationRegion_datadriven_FuncDown_mu;
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_FuncUp_ele, h_ValidationRegion_datadriven_FuncDown_ele;
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_BtagUp_total, h_ValidationRegion_datadriven_BtagDown_total;
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_BtagUp_mu, h_ValidationRegion_datadriven_BtagDown_mu;
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_BtagUp_ele, h_ValidationRegion_datadriven_BtagDown_ele;


  std::unique_ptr<Hists> h_ttbarCR_v1, h_ttbarCR_v2, h_ttbarCR_v3;

  std::unique_ptr<uhh2::TstarTstarDDTHists> h_DDTtestHists;

  // ###### Control switches ######
  bool debug = false;
  bool do_masspoint = false;
  bool doAddInputs = false;


  // ###### Handles ######
  uhh2::Event::Handle<int> h_flag_toptagevent;
  uhh2::Event::Handle<int> h_flag_muonevent;
  uhh2::Event::Handle<double> h_evt_weight;
  uhh2::Event::Handle<bool> h_is_btagevent;


  uhh2::Event::Handle<double> h_ST_weight;

  uhh2::Event::Handle<double> h_DNN_output;
  uhh2::Event::Handle<bool> h_do_masspoint;
  uhh2::Event::Handle<double> h_ST;
  uhh2::Event::Handle<double> h_ST_HOTVR;
  uhh2::Event::Handle<bool> h_DoAddInputs;
  uhh2::Event::Handle<bool> h_is_mu;
  uhh2::Event::Handle<double> h_newTagger;
  uhh2::Event::Handle<TString> h_region;

  uhh2::Event::Handle<float> h_weight_sfelec_triggerNominal;


  // ###### other parameters ######
  bool is_MC;
  bool is_TTbar = false;
  bool is_Signal = false;

  bool is_datadriven_BG_run = false;

  TGraphAsymmErrors* bgest_purity;
  TF1* backgroundEstimationFunctionNominal_SR, * backgroundEstimationFunctionNominal_VR;
  TF1* backgroundEstimationFunctionFuncUp_SR, * backgroundEstimationFunctionFuncUp_VR;
  TF1* backgroundEstimationFunctionFuncDown_SR, * backgroundEstimationFunctionFuncDown_VR;
  TF1* backgroundEstimationFunctionBtagUp_SR, * backgroundEstimationFunctionBtagUp_VR;
  TF1* backgroundEstimationFunctionBtagDown_SR, * backgroundEstimationFunctionBtagDown_VR;
  std::vector<TF1*> DDTFunctions;
  TF1* BestDDTFunction;

  std::unique_ptr<uhh2::AnalysisModule> TstarTstarSpinSwitcher;

};


TstarTstarDNNModule::TstarTstarDNNModule(Context & ctx){

  // setting debug from xml file
  if(ctx.get("debug", "<not set>") == "true") debug = true;

  // debug message
  if(debug) {
    cout << "Hello World from TstarTstarDNNModule!" << endl;
    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }
  }

  // ###### 0. Setting variables ######
  // MC or real data
  is_MC = ctx.get("dataset_type") == "MC";

  if(!is_MC) is_datadriven_BG_run = ctx.get("use_data_for", "regular") == "background_extrapolation"; // get which running mode to use for data

  if(is_datadriven_BG_run) {
    TString background_estimation_purity_filepath = ctx.get("background_estimation_purity_file");
    TFile *f = new TFile(background_estimation_purity_filepath);
    if(!f) throw std::runtime_error("ERROR: cant open background estimation purity at "+ctx.get("background_estimation_purity_file"));
    bgest_purity = (TGraphAsymmErrors*)f->Get("purity");
  }

  if(is_MC) {
    is_TTbar = (ctx.get("dataset_version").find("TT") != std::string::npos);
    is_Signal = (ctx.get("dataset_version").find("Tstar") != std::string::npos);
  }

  // ###### 1. set up modules ######
  // primlep
  reco_primlep.reset(new PrimaryLepton(ctx));

  // DNN modules
  DNN_Includer.reset(new NeuralNetworkIncluder(ctx, do_masspoint));


  // 3. Set up Hists classes:
  h_topcheck.reset(new TstarTstarAllGenHists(ctx, "topcheck"));
  h_topcheck_reweighted.reset(new TstarTstarAllGenHists(ctx, "topcheck_reweighted"));

  h_crosscheck.reset(new TstarTstarHists(ctx, "crosscheck"));
  h_crosscheck_noeletriggerSF.reset(new TstarTstarHists(ctx, "crosscheck_noeletriggerSF"));
  h_STreweighted.reset(new TstarTstarHists(ctx, "STreweighted"));

  h_newTaggerSR.reset(new TstarTstarHists(ctx, "newTaggerSR"));
  h_newTaggerCR.reset(new TstarTstarHists(ctx, "newTaggerCR"));
  h_newTaggerCR_noelectriggerSF.reset(new TstarTstarHists(ctx, "newTaggerCR_noElecTriggerSF"));
  h_newTagger_btagCR.reset(new TstarTstarHists(ctx, "newTagger_btagCR"));
  h_newTagger_btagCR_DNNSR.reset(new TstarTstarHists(ctx, "newTagger_btagCR_ele"));

  /**
  h_AfterDNNcut_02.reset(new TstarTstarHists(ctx, "AfterDNNcut_02"));
  h_notDNNcut_02.reset(new TstarTstarHists(ctx, "notDNNcut_02"));
  h_AfterDNNcut_03.reset(new TstarTstarHists(ctx, "AfterDNNcut_03"));
  h_notDNNcut_03.reset(new TstarTstarHists(ctx, "notDNNcut_03"));
  h_AfterDNNcut_04.reset(new TstarTstarHists(ctx, "AfterDNNcut_04"));
  h_notDNNcut_04.reset(new TstarTstarHists(ctx, "notDNNcut_04"));
  h_AfterDNNcut_05.reset(new TstarTstarHists(ctx, "AfterDNNcut_05"));
  h_notDNNcut_05.reset(new TstarTstarHists(ctx, "notDNNcut_05"));
  h_AfterDNNcut_06.reset(new TstarTstarHists(ctx, "AfterDNNcut_06"));
  h_notDNNcut_06.reset(new TstarTstarHists(ctx, "notDNNcut_06"));
  h_AfterDNNcut_07.reset(new TstarTstarHists(ctx, "AfterDNNcut_07"));
  h_notDNNcut_07.reset(new TstarTstarHists(ctx, "notDNNcut_07"));
  h_AfterDNNcut_08.reset(new TstarTstarHists(ctx, "AfterDNNcut_08"));
  h_notDNNcut_08.reset(new TstarTstarHists(ctx, "notDNNcut_08"));
  **/

  h_AfterDNNcut_06_UGLYFIX.reset(new TstarTstarHists(ctx, "AfterDNNcut_06_UGLYFIX"));
  h_NotDNNcut_06_UGLYFIX.reset(new TstarTstarHists(ctx, "NotDNNcut_06_UGLYFIX"));

  /**
  h_BtagControl_AfterDNNcut_02.reset(new TstarTstarHists(ctx, "AfterDNNcut_02_BtagControl"));
  h_BtagControl_notDNNcut_02.reset(new TstarTstarHists(ctx, "BtagControl_notDNNcut_02"));
  h_BtagControl_AfterDNNcut_03.reset(new TstarTstarHists(ctx, "BtagControl_AfterDNNcut_03"));
  h_BtagControl_notDNNcut_03.reset(new TstarTstarHists(ctx, "BtagControl_notDNNcut_03"));
  h_BtagControl_AfterDNNcut_04.reset(new TstarTstarHists(ctx, "BtagControl_AfterDNNcut_04"));
  h_BtagControl_notDNNcut_04.reset(new TstarTstarHists(ctx, "BtagControl_notDNNcut_04"));
  h_BtagControl_AfterDNNcut_05.reset(new TstarTstarHists(ctx, "BtagControl_AfterDNNcut_05"));
  h_BtagControl_notDNNcut_05.reset(new TstarTstarHists(ctx, "BtagControl_notDNNcut_05"));
  h_BtagControl_AfterDNNcut_06.reset(new TstarTstarHists(ctx, "BtagControl_AfterDNNcut_06"));
  h_BtagControl_notDNNcut_06.reset(new TstarTstarHists(ctx, "BtagControl_notDNNcut_06"));
  h_BtagControl_AfterDNNcut_07.reset(new TstarTstarHists(ctx, "BtagControl_AfterDNNcut_07"));
  h_BtagControl_notDNNcut_07.reset(new TstarTstarHists(ctx, "BtagControl_notDNNcut_07"));
  h_BtagControl_AfterDNNcut_08.reset(new TstarTstarHists(ctx, "BtagControl_AfterDNNcut_08"));
  h_BtagControl_notDNNcut_08.reset(new TstarTstarHists(ctx, "BtagControl_notDNNcut_08"));
  **/

  h_highLepton.reset(new TstarTstarHists(ctx, "highLepton"));
  h_highLepton_AfterDNNcut_06.reset(new TstarTstarHists(ctx, "highLepton_AfterDNNcut_06"));
  h_highLepton_notDNNcut_06.reset(new TstarTstarHists(ctx, "highLepton_notDNNcut_06"));

  h_DNN.reset(new TstarTstarDNNHists(ctx, "DNN"));
  h_DNN_newTagger.reset(new TstarTstarDNNHists(ctx, "DNN_newTagger"));
  h_DNN_BtagControl.reset(new TstarTstarDNNHists(ctx, "DNN_BtagControl"));
  h_DNN_reweighted.reset(new TstarTstarDNNHists(ctx, "DNN_reweighted"));

  h_AfterDNN.reset(new TstarTstarHists(ctx, "AfterDNN"));

  h_SignalRegion_total.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_total"));
  h_SignalRegion_mu.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_mu"));
  h_SignalRegion_ele.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_ele"));

  h_ValidationRegion_total.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_total"));
  h_ValidationRegion_mu.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_mu"));
  h_ValidationRegion_ele.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_ele"));

  h_ControlRegion_total.reset(new TstarTstarSignalRegionHists(ctx, "ControlRegion_total"));
  h_ControlRegion_mu.reset(new TstarTstarSignalRegionHists(ctx, "ControlRegion_mu"));
  h_ControlRegion_ele.reset(new TstarTstarSignalRegionHists(ctx, "ControlRegion_ele"));

  h_SignalRegion_datadriven_FuncUp_total.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_FuncUp_total"));
  h_SignalRegion_datadriven_FuncDown_total.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_FuncDown_total"));
  h_ValidationRegion_datadriven_FuncUp_total.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_FuncUp_total")); // using SR hists here, althoigh it is not
  h_ValidationRegion_datadriven_FuncDown_total.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_FuncDown_total")); // using SR hists here, althoigh it is not

  h_SignalRegion_datadriven_FuncUp_mu.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_FuncUp_mu"));
  h_SignalRegion_datadriven_FuncDown_mu.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_FuncDown_mu"));
  h_ValidationRegion_datadriven_FuncUp_mu.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_FuncUp_mu")); // using SR hists here, althoigh it is not
  h_ValidationRegion_datadriven_FuncDown_mu.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_FuncDown_mu")); // using SR hists here, althoigh it is not

  h_SignalRegion_datadriven_FuncUp_ele.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_FuncUp_ele"));
  h_SignalRegion_datadriven_FuncDown_ele.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_FuncDown_ele"));
  h_ValidationRegion_datadriven_FuncUp_ele.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_FuncUp_ele")); // using SR hists here, althoigh it is not
  h_ValidationRegion_datadriven_FuncDown_ele.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_FuncDown_ele")); // using SR hists here, althoigh it is not

  h_SignalRegion_datadriven_BtagUp_total.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_BtagUp_total"));
  h_SignalRegion_datadriven_BtagDown_total.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_BtagDown_total"));
  h_ValidationRegion_datadriven_BtagUp_total.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_BtagUp_total")); // using SR hists here, althoigh it is not
  h_ValidationRegion_datadriven_BtagDown_total.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_BtagDown_total")); // using SR hists here, althoigh it is not

  h_SignalRegion_datadriven_BtagUp_mu.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_BtagUp_mu"));
  h_SignalRegion_datadriven_BtagDown_mu.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_BtagDown_mu"));
  h_ValidationRegion_datadriven_BtagUp_mu.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_BtagUp_mu")); // using SR hists here, althoigh it is not
  h_ValidationRegion_datadriven_BtagDown_mu.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_BtagDown_mu")); // using SR hists here, althoigh it is not

  h_SignalRegion_datadriven_BtagUp_ele.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_BtagUp_ele"));
  h_SignalRegion_datadriven_BtagDown_ele.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_BtagDown_ele"));
  h_ValidationRegion_datadriven_BtagUp_ele.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_BtagUp_ele")); // using SR hists here, althoigh it is not
  h_ValidationRegion_datadriven_BtagDown_ele.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_BtagDown_ele")); // using SR hists here, althoigh it is not

  h_ttbarCR_v1.reset(new TstarTstarHists(ctx, "ttbarCR_v1"));
  h_ttbarCR_v2.reset(new TstarTstarHists(ctx, "ttbarCR_v2"));
  h_ttbarCR_v3.reset(new TstarTstarHists(ctx, "ttbarCR_v3"));

  // ###### 4. init handles ######
  ctx.undeclare_all_event_output();
  h_evt_weight = ctx.get_handle<double>("evt_weight");
  h_flag_toptagevent = ctx.get_handle<int>("flag_toptagevent");
  h_DoAddInputs = ctx.declare_event_output<bool>("doAddInputs");
  h_newTagger = ctx.declare_event_output<double>("newTagger");
  h_is_btagevent = ctx.get_handle<bool>("is_btagevent");
  h_is_mu = ctx.get_handle<bool>("is_muevt");
  h_ST = ctx.get_handle<double>("ST");
  h_ST_HOTVR = ctx.get_handle<double>("STHOTVR");
  h_region = ctx.declare_event_output<TString>("region");

  h_ST_weight = ctx.declare_event_output<double>("ST_weight");

  h_DNN_output = ctx.get_handle<double>("DNN_output");

  if(is_datadriven_BG_run) {
    TString path = "/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/";
    TString filename_nominal_SR = "alphaFunction_HOTVR_SR.root";
    TString filename_nominal_VR = "alphaFunction_HOTVR_VR.root";
    TString filename_btagUp_SR = "alphaFunction_HOTVR_SR_btagging_totalUp.root";
    TString filename_btagDown_SR = "alphaFunction_HOTVR_SR_btagging_totalDown.root";
    TString filename_btagUp_VR = "alphaFunction_HOTVR_VR_btagging_totalUp.root";
    TString filename_btagDown_VR = "alphaFunction_HOTVR_VR_btagging_totalDown.root";

    // main file for signal region
    TFile *file_nominal_SR = new TFile(path+filename_nominal_SR);
    backgroundEstimationFunctionNominal_SR = (TF1*)file_nominal_SR->Get("fit_mean");
    backgroundEstimationFunctionFuncUp_SR = (TF1*)file_nominal_SR->Get("fit1");
    backgroundEstimationFunctionFuncDown_SR = (TF1*)file_nominal_SR->Get("fit2");

    // main file for validation region
    TFile *file_nominal_VR = new TFile(path+filename_nominal_VR);
    backgroundEstimationFunctionNominal_VR = (TF1*)file_nominal_VR->Get("fit_mean");
    backgroundEstimationFunctionFuncUp_VR = (TF1*)file_nominal_VR->Get("fit1");
    backgroundEstimationFunctionFuncDown_VR = (TF1*)file_nominal_VR->Get("fit2");

    // now getting the btag variations
    TFile *file_btagUp_SR = new TFile(path+filename_btagUp_SR);
    backgroundEstimationFunctionBtagUp_SR = (TF1*)file_btagUp_SR->Get("fit_mean");
    TFile *file_btagDown_SR = new TFile(path+filename_btagDown_SR);
    backgroundEstimationFunctionBtagDown_SR = (TF1*)file_btagDown_SR->Get("fit_mean");
    TFile *file_btagUp_VR = new TFile(path+filename_btagUp_VR);
    backgroundEstimationFunctionBtagUp_VR = (TF1*)file_btagUp_VR->Get("fit_mean");
    TFile *file_btagDown_VR = new TFile(path+filename_btagDown_VR);
    backgroundEstimationFunctionBtagDown_VR = (TF1*)file_btagDown_VR->Get("fit_mean");
  }

  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/";
  TString bestFunction = "0p3";
  TString filename_base = "DDTfunc_";
  
  // getting the best function
  TFile *f = new TFile(path+filename_base+bestFunction+".root");
  BestDDTFunction = (TF1*)f->Get("fit");

  // now lets load all the others
  std::vector<TString> points_to_check = {"0p1", "0p15", "0p2", "0p25", "0p3", "0p35", "0p4", "0p45", "0p5"};
  for (auto point : points_to_check) {
    TFile *f = new TFile(path+filename_base+point+".root");
    auto DDTFunction = (TF1*)f->Get("fit");
    DDTFunctions.push_back(DDTFunction);
  }

  h_DDTtestHists.reset( new TstarTstarDDTHists(ctx, "DDTHists", points_to_check) ) ;

  h_weight_sfelec_triggerNominal = ctx.get_handle<float>("weight_sfelec_trigger");

  TstarTstarSpinSwitcher.reset(new TstarTstarSpinScale(ctx, "/nfs/dust/cms/user/flabe/TstarTstar/Jupyter/reweight_factors/"));

}


bool TstarTstarDNNModule::process(Event & event) {

  // debug message
  if(debug){cout << endl << "TstarTstarDNNModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;}

  // fixing poitential empty handle
  event.set(h_region,"not set");

  // TODO remove me
  // if (event.get(h_ST_HOTVR) < 600) return false; // TODO REMOVE ME IN NEXT ITERATION
  // std::cout << "REMOVE TEMP FIX HERE!!!" << std::endl;

  // reapply weights
  event.weight = event.get(h_evt_weight);
  if(debug) cout << "weights applied." << endl;

  if(!(TstarTstarSpinSwitcher->process(event))) return false; // this will reweight, but only if set so in config file
  event.set(h_evt_weight, event.weight); // we'll use this later, so need to re-set

  // get ST weights
  double ST_weight = event.get(h_ST_weight);

  // set primary lepton
  reco_primlep->process(event);

  // setting addInputs
  event.set(h_DoAddInputs, doAddInputs);

  // placeholder as this flag seems to be wrong!!!
  BTag bJetID = BTag(BTag::algo::DEEPJET, BTag::wp::WP_MEDIUM);
  bool pass_btagcut = false;
  for (const auto & jet: *event.jets){
    if(bJetID(jet, event)) pass_btagcut = true;
  }

  h_crosscheck->fill(event);

  double reset_weight_etriggera = event.weight;
  if( event.get(h_weight_sfelec_triggerNominal) != 0) event.weight /= event.get(h_weight_sfelec_triggerNominal);
  h_crosscheck_noeletriggerSF->fill(event);
  event.weight = reset_weight_etriggera;
  
  event.weight = ST_weight * event.get(h_evt_weight);
  if(pass_btagcut) h_STreweighted->fill(event);
  event.weight = event.get(h_evt_weight);


  // ################
  // ### DNN Part ###
  // ################

  // including trained DNN model
  if(debug) cout << "Include DNN model" << endl;
  DNN_Includer->process(event);

  if(debug) cout << "Done DNN include" << endl;
  if(event.get(h_DNN_output) > 0.6) h_AfterDNNcut_06_UGLYFIX->fill(event);
  else h_NotDNNcut_06_UGLYFIX->fill(event);

  if(debug) cout << "Start filling hists" << endl;
  // hists
  if(pass_btagcut) {
    h_DNN->fill(event);
    h_topcheck->fill(event);
    if(is_TTbar) event.weight *= ST_weight;
    h_DNN_reweighted->fill(event);
    h_topcheck_reweighted->fill(event);
    event.weight = event.get(h_evt_weight);

    // some more plotting
    double DNNoutput = event.get(h_DNN_output);
    h_AfterDNN->fill(event);

  }

  /**
  else {
    h_DNN_BtagControl->fill(event);
    if(event.get(h_DNN_output) > 0.2) h_BtagControl_AfterDNNcut_02->fill(event);
    else h_BtagControl_notDNNcut_02->fill(event);
    if(event.get(h_DNN_output) > 0.3) h_BtagControl_AfterDNNcut_03->fill(event);
    else h_BtagControl_notDNNcut_03->fill(event);
    if(event.get(h_DNN_output) > 0.4) h_BtagControl_AfterDNNcut_04->fill(event);
    else h_BtagControl_notDNNcut_04->fill(event);
    if(event.get(h_DNN_output) > 0.5) h_BtagControl_AfterDNNcut_05->fill(event);
    else h_BtagControl_notDNNcut_05->fill(event);
    if(event.get(h_DNN_output) > 0.6) h_BtagControl_AfterDNNcut_06->fill(event);
    else h_BtagControl_notDNNcut_06->fill(event);
    if(event.get(h_DNN_output) > 0.7) h_BtagControl_AfterDNNcut_07->fill(event);
    else h_BtagControl_notDNNcut_07->fill(event);
    if(event.get(h_DNN_output) > 0.8) h_BtagControl_AfterDNNcut_08->fill(event);
    else h_BtagControl_notDNNcut_08->fill(event);
  }
  **/

  double newFunc = BestDDTFunction->Eval(event.get(h_ST));
  double newTagger = event.get(h_DNN_output) - ( 1 - newFunc );
  event.set(h_newTagger, newTagger);
  
  // lets store all the other tagger outputs in a vector
  std::vector<double> newTaggerResults;
  for (auto function : DDTFunctions) {
    double value = event.get(h_DNN_output) - ( 1 - function->Eval(event.get(h_ST)) );
    newTaggerResults.push_back(value);
  }

  // datadriven background estimation
  double transfer_weight_nominal_SR = 1;
  double transfer_weight_nominal_VR = 1;
  double transfer_weight_funcUp_SR = 1;
  double transfer_weight_funcUp_VR = 1;
  double transfer_weight_funcDown_SR = 1;
  double transfer_weight_funcDown_VR = 1;
  double transfer_weight_btagUp_SR = 1;
  double transfer_weight_btagUp_VR = 1;
  double transfer_weight_btagDown_SR = 1;
  double transfer_weight_btagDown_VR = 1;

  double weight_for_resetting = event.weight; 
  double purity_value = 1;
  if(is_datadriven_BG_run && !pass_btagcut) {

    if(debug) cout << "Doing datadriven BG estimation" << endl;
    // getting all the weight values
    transfer_weight_nominal_SR = backgroundEstimationFunctionNominal_SR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_nominal_VR = backgroundEstimationFunctionNominal_VR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_funcUp_SR = backgroundEstimationFunctionFuncUp_SR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_funcUp_VR = backgroundEstimationFunctionFuncUp_VR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_funcDown_SR = backgroundEstimationFunctionFuncDown_SR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_funcDown_VR = backgroundEstimationFunctionFuncDown_VR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_btagUp_SR = backgroundEstimationFunctionBtagUp_SR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_btagUp_VR = backgroundEstimationFunctionBtagUp_VR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_btagDown_SR = backgroundEstimationFunctionBtagDown_SR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_btagDown_VR = backgroundEstimationFunctionBtagDown_VR->Eval(event.get(h_ST_HOTVR));
    if(debug) cout << "Gotten all values" << endl;
    
    purity_value = bgest_purity->Eval(event.get(h_ST_HOTVR));
    if(debug) cout << "purity: " << purity_value << endl;

    // fill all the up and down-variations

    // function
    event.weight *= transfer_weight_funcUp_SR * purity_value;
    h_SignalRegion_datadriven_FuncUp_total->fill(event);
    if(event.get(h_is_mu)) h_SignalRegion_datadriven_FuncUp_mu->fill(event);
    else h_SignalRegion_datadriven_FuncUp_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_funcDown_SR * purity_value;
    h_SignalRegion_datadriven_FuncDown_total->fill(event);
    if(event.get(h_is_mu)) h_SignalRegion_datadriven_FuncDown_mu->fill(event);
    else h_SignalRegion_datadriven_FuncDown_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_funcUp_VR * purity_value;
    h_ValidationRegion_datadriven_FuncUp_total->fill(event);
    if(event.get(h_is_mu)) h_ValidationRegion_datadriven_FuncUp_mu->fill(event);
    else h_ValidationRegion_datadriven_FuncUp_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_funcDown_VR * purity_value;
    h_ValidationRegion_datadriven_FuncDown_total->fill(event);
    if(event.get(h_is_mu)) h_ValidationRegion_datadriven_FuncDown_mu->fill(event);
    else h_ValidationRegion_datadriven_FuncDown_ele->fill(event);
    event.weight = weight_for_resetting;

    // btag
    event.weight *= transfer_weight_btagUp_SR * purity_value;
    h_SignalRegion_datadriven_BtagUp_total->fill(event);
    if(event.get(h_is_mu)) h_SignalRegion_datadriven_BtagUp_mu->fill(event);
    else h_SignalRegion_datadriven_BtagUp_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_btagDown_SR * purity_value;
    h_SignalRegion_datadriven_BtagDown_total->fill(event);
    if(event.get(h_is_mu)) h_SignalRegion_datadriven_BtagDown_mu->fill(event);
    else h_SignalRegion_datadriven_BtagDown_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_btagUp_VR * purity_value;
    h_ValidationRegion_datadriven_BtagUp_total->fill(event);
    if(event.get(h_is_mu)) h_ValidationRegion_datadriven_BtagUp_mu->fill(event);
    else h_ValidationRegion_datadriven_BtagUp_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_btagDown_VR * purity_value;
    h_ValidationRegion_datadriven_BtagDown_total->fill(event);
    if(event.get(h_is_mu)) h_ValidationRegion_datadriven_BtagDown_mu->fill(event);
    else h_ValidationRegion_datadriven_BtagDown_ele->fill(event);
    event.weight = weight_for_resetting;

  }

  if(pass_btagcut || (is_datadriven_BG_run && !pass_btagcut)) {

    h_DNN_newTagger->fill(event);

    bool fillSR = false;
    bool fillCR = false;

    if(is_datadriven_BG_run && !pass_btagcut) {
      fillSR = true;
      fillCR = true;
      event.set(h_region, "datadrivenEst");
    } else if (!is_datadriven_BG_run) {
      if(newTagger > 0) {
        fillSR = true;
        fillCR = false;
        event.set(h_region, "SR");
      } else {
        fillSR = false;
        fillCR = true;
        event.set(h_region, "VR");
      }
    }
   
    if(fillSR) {
      if (is_datadriven_BG_run) event.weight *= transfer_weight_nominal_SR * purity_value;
      if(is_MC /* blinding */ || is_datadriven_BG_run) {
        if(debug) std::cout << "Filling SR" << std::endl;
        h_newTaggerSR->fill(event);
        h_SignalRegion_total->fill(event);
        if(event.get(h_is_mu)) h_SignalRegion_mu->fill(event);
        else h_SignalRegion_ele->fill(event);
      }
      if (is_datadriven_BG_run) event.weight = weight_for_resetting;
    }
    if(fillCR) {
      if (is_datadriven_BG_run) event.weight *= transfer_weight_nominal_VR * purity_value;
      h_newTaggerCR->fill(event);
      h_ValidationRegion_total->fill(event);
        if(event.get(h_is_mu)) h_ValidationRegion_mu->fill(event);
        else {

          double reset_weight_etrigger = event.weight;
          if( event.get(h_weight_sfelec_triggerNominal) != 0) event.weight /= event.get(h_weight_sfelec_triggerNominal);
          h_newTaggerCR_noelectriggerSF->fill(event);
          event.weight = reset_weight_etrigger;

          h_ValidationRegion_ele->fill(event);
        }
      if (is_datadriven_BG_run) event.weight = weight_for_resetting;

      // a subsest of this will be used as a (to-test) ttbar CR

      // HOTVR top tag
      TopJetId topjetID = AndId<TopJet>(HOTVRTopTag(), Tau32Groomed(0.56)); // Top Tag that is used later
      bool passHOTVRtoptag = false;
      for (const auto & jet: *event.topjets){
        if(topjetID(jet, event)) passHOTVRtoptag = true;
      }

      if(passHOTVRtoptag) {

        // tighter b-tag selection
        // ###### Btag Selection ######
        BTag bJetID_loose = BTag(BTag::algo::DEEPJET, BTag::wp::WP_LOOSE);
        BTag bJetID_tight = BTag(BTag::algo::DEEPJET, BTag::wp::WP_TIGHT);
        int N_btag_loose = 0;
        int N_btag_tight = 0;
        for (const auto & jet: *event.jets){
          if(bJetID_loose(jet, event)) N_btag_loose++;
          if(bJetID_tight(jet, event)) N_btag_tight++;
        }

        // try three different versions of this CR
        if(N_btag_loose >= 2) h_ttbarCR_v1->fill(event);
        if(N_btag_tight >= 1) h_ttbarCR_v2->fill(event);
        if(N_btag_tight >= 2) h_ttbarCR_v3->fill(event);

      }


    }

    // lets fill some other histogram class that will give us the output for each point to check
    h_DDTtestHists->fill(event, newTaggerResults);

  } else {
    h_newTagger_btagCR->fill(event);

    h_ControlRegion_total->fill(event);
    if(event.get(h_is_mu)) h_ControlRegion_mu->fill(event);
    else h_ControlRegion_ele->fill(event);

    if(newTagger > 0) {
      h_newTagger_btagCR_DNNSR->fill(event);
      event.set(h_region,"CR");
    }
  }

  // testing tighter pt cut
  if(debug) cout << "testing tighter pt cut" << endl;
  if(pass_btagcut) {
    if(event.muons->size() == 1) {
      h_highLepton->fill(event);
      if(event.get(h_DNN_output) > 0.6) h_highLepton_AfterDNNcut_06->fill(event);
      else h_highLepton_notDNNcut_06->fill(event);
    }
    else {
      if(event.electrons->at(0).pt() > 50 && event.met->pt() > 80) {
        h_highLepton->fill(event);
        if(event.get(h_DNN_output) > 0.6) h_highLepton_AfterDNNcut_06->fill(event);
        else h_highLepton_notDNNcut_06->fill(event);
      }
    }

  }

  if(debug){cout << "Done ##################################" << endl;}
  return true;

}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarDNNModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarDNNModule)

}
