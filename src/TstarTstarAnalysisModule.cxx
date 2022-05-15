#include <iostream>
#include <memory>
#include <string>

// UHH2 stuff
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/PhotonIds.h"
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/TriggerSelection.h>
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/MCWeight.h"

// TstarTstar stuff
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGENTstarRecoHists.h"
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

// other stuff
#include "UHH2/HOTVR/include/HOTVRIds.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

// quick method to calculate inv_mass
float inv_mass_3(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

/** \brief Module for the T*T*->ttbar gg MC based study
 *
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class TstarTstarAnalysisModule: public AnalysisModule {
public:

  explicit TstarTstarAnalysisModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  // ###### Modules ######
  unique_ptr<Selection> toptagevt_sel;
  unique_ptr<NeuralNetworkInputWriter> DNN_InputWriter;

  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;

  std::unique_ptr<TstarTstar_tgtg_TopTag_Reconstruction> Tstarreco_gHOTVR;
  std::unique_ptr<TstarTstar_tgtg_AK4_Reconstruction> Tstarreco_gAK4;
  std::unique_ptr<TstarTstar_Discrimination> TstarDiscriminator_gHOTVR;
  std::unique_ptr<TstarTstar_Discrimination> TstarDiscriminator_gAK4;
  uhh2::Event::Handle<ReconstructionTstarHypothesis> h_tstartstar_hyp;
  uhh2::Event::Handle<ReconstructionTstarHypothesis> h_tstartstar_hyp_gHOTVR;
  uhh2::Event::Handle<ReconstructionTstarHypothesis> h_tstartstar_hyp_gAK4;
  uhh2::Event::Handle<ReconstructionTstarHypothesis> h_tstartstar_GENhyp_gHOTVR;
  uhh2::Event::Handle<ReconstructionTstarHypothesis> h_tstartstar_GENhyp_gAK4;
  std::unique_ptr<TstarTstarGenDiscriminator> TstarGENDiscriminator_gHOTVR;
  std::unique_ptr<TstarTstarGenDiscriminator> TstarGENDiscriminator_gAK4;


  // ##### Histograms #####
  std::unique_ptr<Hists> h_main, h_main_ttag, h_main_nottag, h_main_mu, h_main_mu_lowpt, h_main_mu_highpt, h_main_ele, h_main_ele_lowpt, h_main_ele_highpt;
  std::unique_ptr<Hists> h_reco, h_reco_ttag, h_reco_nottag, h_reco_mu, h_reco_mu_lowpt, h_reco_mu_highpt, h_reco_ele, h_reco_ele_lowpt, h_reco_ele_highpt;
  std::unique_ptr<Hists> h_crosscheck, h_STreweighted;
  std::unique_ptr<Hists> h_main_gen;

  std::unique_ptr<Hists> h_DNN_Inputs, h_DNN_Inputs_reweighted;

  std::unique_ptr<Hists> h_GENTstarReco;

  // ###### Handles ######
  uhh2::Event::Handle<double> h_evt_weight;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
  uhh2::Event::Handle<double> h_ST;
  uhh2::Event::Handle<double> h_STHOTVR;
  uhh2::Event::Handle<LorentzVector> h_neutrino;
  uhh2::Event::Handle<bool> h_is_btagevent;

  // for reconstruction
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<int> h_flag_toptagevent;
  uhh2::Event::Handle<int> h_flag_muonevent;

  // for DNN output
  uhh2::Event::Handle<bool> h_do_masspoint;

  // ###### Control Switches ######
  bool debug = false;
  bool outputDNNvalues = true;
  bool do_masspoint = false;


  // ###### other needed definitions ######
  bool isTrigger;
  bool is_MC;

  bool is_TTbar;
  bool is_Signal;
  bool is_Data;

  TH1D* ST_ratio;
  uhh2::Event::Handle<double> h_ST_weight;
};


TstarTstarAnalysisModule::TstarTstarAnalysisModule(Context & ctx){

  // setting debug from xml file
  if(ctx.get("debug", "<not set>") == "true") debug = true;

  // debug messagt
  if(debug) {
    cout << "Hello World from TstarTstarAnalysisModule!" << endl;
    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }
  }

  // ###### 0. Setting Variables ######
  // MC or real data
  is_MC = ctx.get("dataset_type") == "MC";

  is_TTbar = (ctx.get("dataset_version").find("TT") != std::string::npos);
  is_Signal = (ctx.get("dataset_version").find("Tstar") != std::string::npos);


  // ###### 1. Set up modules ######
  // primary lepton
  reco_primlep.reset(new PrimaryLepton(ctx));

  // GEN things
  if(is_MC){
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  }

  // top tag definition
  TopJetId topjetID = AndId<TopJet>(HOTVRTopTag(), Tau32Groomed(0.56));
  toptagevt_sel.reset(new TopTagEventSelection(topjetID));;


  // ###### 3. Set up histograms ######
  // before Reconstruction
  h_crosscheck.reset(new TstarTstarHists(ctx, "crosscheck"));
  h_main.reset(new TstarTstarHists(ctx, "main"));
  h_main_ttag.reset(new TstarTstarHists(ctx, "main_ttag"));
  h_main_nottag.reset(new TstarTstarHists(ctx, "main_nottag"));
  h_main_mu.reset(new TstarTstarHists(ctx, "main_mu"));
  h_main_mu_lowpt.reset(new TstarTstarHists(ctx, "main_mu_lowpt"));
  h_main_mu_highpt.reset(new TstarTstarHists(ctx, "main_mu_highpt"));
  h_main_ele.reset(new TstarTstarHists(ctx, "main_ele"));
  h_main_ele_lowpt.reset(new TstarTstarHists(ctx, "main_ele_lowpt"));
  h_main_ele_highpt.reset(new TstarTstarHists(ctx, "main_ele_highpt"));
  h_main_gen.reset(new TstarTstarGenHists(ctx, "main_gen"));

  // after Reconstruction
  h_reco.reset(new TstarTstarHists(ctx, "reco"));
  h_reco_ttag.reset(new TstarTstarHists(ctx, "reco_ttag"));
  h_reco_nottag.reset(new TstarTstarHists(ctx, "reco_nottag"));
  h_reco_mu.reset(new TstarTstarHists(ctx, "reco_mu"));
  h_reco_mu_lowpt.reset(new TstarTstarHists(ctx, "reco_mu_lowpt"));
  h_reco_mu_highpt.reset(new TstarTstarHists(ctx, "reco_mu_highpt"));
  h_reco_ele.reset(new TstarTstarHists(ctx, "reco_ele"));
  h_reco_ele_lowpt.reset(new TstarTstarHists(ctx, "reco_ele_lowpt"));
  h_reco_ele_highpt.reset(new TstarTstarHists(ctx, "reco_ele_highpt"));

  h_STreweighted.reset(new TstarTstarHists(ctx, "STreweighted"));

  // GEN reco
  h_GENTstarReco.reset(new TstarTstarGENTstarRecoHists(ctx, "GENTstarReco"));

  // DNN hists
  h_DNN_Inputs.reset(new TstarTstarDNNInputHists(ctx, "DNN_Inputs"));
  h_DNN_Inputs_reweighted.reset(new TstarTstarDNNInputHists(ctx, "DNN_Inputs_reweighted"));

  // ###### 4. Init handles ######
  h_evt_weight = ctx.get_handle<double>("evt_weight");
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_flag_muonevent = ctx.declare_event_output<int>("flag_muonevent");
  h_flag_toptagevent = ctx.declare_event_output<int>("flag_toptagevent");
  h_neutrino = ctx.get_handle<LorentzVector>("neutrino");
  h_is_btagevent = ctx.get_handle<bool>("is_btagevent");

  if(is_MC) h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
  h_ST = ctx.get_handle<double>("ST");
  h_STHOTVR = ctx.get_handle<double>("STHOTVR");

  // DNN output
  if(outputDNNvalues){
    DNN_InputWriter.reset(new NeuralNetworkInputWriter(ctx));
    h_do_masspoint = ctx.get_handle<bool>("do_masspoint");
    h_ST_weight = ctx.declare_event_output<double>("ST_weight");
  }

  TFile *f = new TFile("/nfs/dust/cms/user/flabe/MLCorner/TstarNN/reweightingApproach/output/data/ST_weights.root");
  ST_ratio = (TH1D*)f->Get("ST_ratio");

  // TstarTstar reconstruction
  Tstarreco_gHOTVR.reset(new TstarTstar_tgtg_TopTag_Reconstruction(ctx, NeutrinoReconstruction, topjetID));
  Tstarreco_gAK4.reset(new TstarTstar_tgtg_AK4_Reconstruction(ctx, NeutrinoReconstruction, topjetID));
  TstarDiscriminator_gHOTVR.reset(new TstarTstar_Discrimination(ctx));
  TstarDiscriminator_gAK4.reset(new TstarTstar_Discrimination(ctx));
  h_tstartstar_hyp = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_Hyp");
  h_tstartstar_hyp_gHOTVR = ctx.declare_event_output<ReconstructionTstarHypothesis>("TstarTstar_Hyp_gHOTVR");
  h_tstartstar_hyp_gAK4 = ctx.declare_event_output<ReconstructionTstarHypothesis>("TstarTstar_Hyp_gAK4");
  TstarGENDiscriminator_gAK4.reset(new TstarTstarGenDiscriminator(ctx));
  TstarGENDiscriminator_gHOTVR.reset(new TstarTstarGenDiscriminator(ctx));
  h_tstartstar_GENhyp_gHOTVR = ctx.declare_event_output<ReconstructionTstarHypothesis>("TstarTstar_GENHyp_gHOTVR");
  h_tstartstar_GENhyp_gAK4 = ctx.declare_event_output<ReconstructionTstarHypothesis>("TstarTstar_GENHyp_gAK4");

}


bool TstarTstarAnalysisModule::process(Event & event) {

  // debug message
  if(debug){cout << endl << "TstarTstarAnalysisModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;}

  // reapply weights
  event.weight = event.get(h_evt_weight);
  if(debug) cout << "weights applied." << endl;

  // set lepton channel
  const bool muon_evt = (event.muons->size() == 1);
  event.set(h_flag_muonevent, int(muon_evt));

  h_crosscheck->fill(event);

  // set primlep
  reco_primlep->process(event);

  // check ttag
  const bool pass_ttag = toptagevt_sel->passes(event);
  event.set(h_flag_toptagevent, int(pass_ttag));

  // ttgen
  if(is_MC) ttgenprod->process(event);

  if(debug) std::cout << "Hists before everything" << std::endl;
  // fill hists before things have happened
  if(event.get(h_is_btagevent)) {
    h_main->fill(event);
    h_main_gen->fill(event);
    if(pass_ttag) h_main_ttag->fill(event);
    else h_main_nottag->fill(event);
    if(event.get(h_flag_muonevent)){
      h_main_mu->fill(event);
      if(event.get(h_primlep).pt()<60) h_main_mu_lowpt->fill(event);
      else h_main_mu_highpt->fill(event);
    }
    else {
      h_main_ele->fill(event);
      if(event.get(h_primlep).pt()<120) h_main_ele_lowpt->fill(event);
      else h_main_ele_highpt->fill(event);
    }
  }

  // ########################################
  // ########### DNN Preparation ############
  // ########################################

  if(debug) cout << "Start DNN stuff" << endl;

  // Filling output for DNN
  if(outputDNNvalues){
    event.set(h_do_masspoint, do_masspoint);
    if(debug) cout << "Write inputs" << endl;
    DNN_InputWriter->process(event);
    if(debug) cout << "plot inputs" << endl;
    if(event.get(h_is_btagevent)) h_DNN_Inputs->fill(event);
    double st = event.get(h_ST);
    double ST_weight = 1;
    if(is_TTbar){
      if(st > 4000) st = 3999;
      double value = ST_ratio->GetBinContent(ST_ratio->GetXaxis()->FindBin(st));
      if(value > 0) ST_weight = 1/value;
      else ST_weight = 0;
    }
    event.set(h_ST_weight, ST_weight);
    event.weight = event.get(h_ST_weight) * event.get(h_evt_weight);
    if(event.get(h_is_btagevent)) h_STreweighted->fill(event);
    if(event.get(h_is_btagevent)) h_DNN_Inputs_reweighted->fill(event);
    event.weight = event.get(h_evt_weight);
  }

  // ########################################
  // ###### TstarTstar Reconstruction #######
  // ########################################

  // gHOTVR case
  if(is_MC){ // blinding
    if(Tstarreco_gHOTVR->process(event)) {
      TstarDiscriminator_gHOTVR->process(event);
      event.set(h_tstartstar_hyp_gHOTVR, event.get(h_tstartstar_hyp));
      TstarGENDiscriminator_gHOTVR->process(event);
      event.set(h_tstartstar_GENhyp_gHOTVR, event.get(h_tstartstar_hyp));
    } else {
      event.set(h_tstartstar_hyp_gHOTVR, ReconstructionTstarHypothesis());
      event.set(h_tstartstar_GENhyp_gHOTVR, ReconstructionTstarHypothesis());
    }
  } else {
    event.set(h_tstartstar_hyp_gHOTVR, ReconstructionTstarHypothesis());
    event.set(h_tstartstar_GENhyp_gHOTVR, ReconstructionTstarHypothesis());
  }

  // gAK4 case
  if(is_MC){ // blinding
    if(Tstarreco_gAK4->process(event)) {
      TstarDiscriminator_gAK4->process(event);
      event.set(h_tstartstar_hyp_gAK4, event.get(h_tstartstar_hyp));
      TstarGENDiscriminator_gAK4->process(event);
      event.set(h_tstartstar_GENhyp_gAK4, event.get(h_tstartstar_hyp));
    } else {
      event.set(h_tstartstar_hyp_gAK4, ReconstructionTstarHypothesis());
      event.set(h_tstartstar_GENhyp_gAK4, ReconstructionTstarHypothesis());
    }
  } else {
    event.set(h_tstartstar_hyp_gAK4, ReconstructionTstarHypothesis());
    event.set(h_tstartstar_GENhyp_gAK4, ReconstructionTstarHypothesis());
  }

  // filling hists after reco
  if(event.get(h_is_btagevent)) {
    h_reco->fill(event);
    h_GENTstarReco->fill(event);
    if(pass_ttag) h_reco_ttag->fill(event);
    else h_reco_nottag->fill(event);
    if(event.get(h_flag_muonevent)){
      h_reco_mu->fill(event);
      if(event.get(h_primlep).pt()<60) h_reco_mu_lowpt->fill(event);
      else h_reco_mu_highpt->fill(event);
    }
    else {
      h_reco_ele->fill(event);
      if(event.get(h_primlep).pt()<120) h_reco_ele_lowpt->fill(event);
      else h_reco_ele_highpt->fill(event);
    }
  }

  if(debug){cout << "Done ##################################" << endl;}
  return true;

}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarAnalysisModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarAnalysisModule)

}
