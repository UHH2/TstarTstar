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
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/TstarTstar/include/ModuleBASE.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"

#include "UHH2/common/include/MCWeight.h"

using namespace std;
using namespace uhh2;

//namespace uhh2 {

/** \brief Module for the full selection in T*T*->ttbar gg/gamma search 
 *
 * All objects are expected to be corrected in PreSelection stage
 *
 */
// class TstarTstarSelectionModule: public AnalysisModule {
class TstarTstarSelectionModule: public ModuleBASE {
public:
    
    explicit TstarTstarSelectionModule(Context & ctx);
    virtual bool process(Event & event) override;
    void book_histograms(uhh2::Context&, vector<string>);
    void fill_histograms(uhh2::Event&, string, bool);
    void book_TTreeVars(uhh2::Context& ctx, vector<string>);
    void fill_TTreeVars(uhh2::Event & event, string var_name, float var_value);

    void fill_MLinput(uhh2::Event & event);

private:
    

    // Declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
    // to avoid memory leaks.

    unique_ptr<Selection> TTbarSemiLepMatchable_selection;
  
    unique_ptr<Selection> triggerPFHT_sel;
    unique_ptr<Selection> triggerSingleJet450_sel;
    unique_ptr<Selection> triggerSingleLeptonEle1_sel;
    unique_ptr<Selection> triggerSingleLeptonEle2_sel;
    unique_ptr<Selection> triggerSingleLeptonEle3_sel;
    unique_ptr<Selection> triggerSingleLeptonMu1_sel;
    unique_ptr<Selection> triggerSingleLeptonMu2_sel;
    unique_ptr<Selection> triggerSingleLeptonMu3_sel;
    unique_ptr<Selection> triggerSingleLeptonMu4_sel;
    unique_ptr<Selection> triggerHT1_sel, triggerHT2_sel, triggerHT3_sel, triggerHT4_sel, triggerHT5_sel,  triggerHT6_sel;
    unique_ptr<Selection>  met_sel, st_sel; //selections defined in UHH2/TstarTstar/include/TstarTstarSelections.h
    unique_ptr<Selection> topjet_selection;

    unique_ptr<AnalysisModule> LumiWeight_module;
    // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  //  std::unique_ptr<Hists> h_nocuts, h_common, h_lepsel, h_2dcut, h_semilepttbarmatch, h_nosemilepttbarmatch;// h_trigger;
  //    std::unique_ptr<Hists> h_semilepttbarmatch_gen;
  //    std::unique_ptr<Hists> h_semilepttbarmatch_genreco;
  bool debug = false;
  //  bool debug = true;

  //Tstar-Tstar reconstruction
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod; 
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;
  std::unique_ptr<uhh2::AnalysisModule> ttbar_reco;
  std::unique_ptr<ttbarChi2Discriminator> ttbar_discriminator;
  //  std::unique_ptr<ttbarCorrectMatchDiscriminator> ttbar_CorrectMatchDiscriminator;
  std::unique_ptr<CorrectMatchDiscriminator> ttbar_CorrectMatchDiscriminator;
  std::unique_ptr<TstarTstar_Reconstruction> TstarTstar_reco;
  std::unique_ptr<TstarTstar_tgluon_tgluon_Reconstruction> TstarTstar_tgluon_tgluon_reco;
  //  uhh2::Event::Handle<std::vector<ReconstructionHypothesis>> h_ttbar_hyps;
  uhh2::Event::Handle<bool> h_is_ttbar_reconstructed;
  uhh2::Event::Handle<ReconstructionHypothesis> h_recohyp;

  uhh2::Event::Handle<std::vector<ReconstructionTstarHypothesis>> h_tstartstar_hyps;
  uhh2::Event::Handle<ReconstructionTstarHypothesis> h_recohyp_tstartstar;

  std::unique_ptr<TstarTstarRecoTstarHists>  h_RecoPlots_After_TstarTstar_tgtg, h_RecoPlots_After_TstarTstar_tgtg_ttbarsemilep;

  bool is_tgtg, is_tgtgamma;
  std::vector<Event::Handle<float>> h_MLinput_values;
  std::vector<string> h_MLinput_names;
  Event::Handle<float> h_weight;

};

void TstarTstarSelectionModule::book_histograms(uhh2::Context& ctx, vector<string> tags){
  for(const auto & tag : tags){
    string  mytag = tag + "_RECO";
    book_HFolder(mytag, new TstarTstarHists(ctx,mytag));
    mytag = tag + "_GEN";
    book_HFolder(mytag, new TstarTstarGenHists(ctx,mytag));
    mytag = tag + "_GENRECO";
    book_HFolder(mytag, new TstarTstarGenRecoMatchedHists(ctx,mytag));
  }
}

  void TstarTstarSelectionModule::fill_histograms(uhh2::Event& event, string tag, bool pass_ttbarsemilep){
    //  for(const auto & tag : tags){
    string mytag = tag + "_RECO";
    HFolder(mytag)->fill(event);
    if(pass_ttbarsemilep){
      mytag = tag + "_GEN";
      HFolder(mytag)->fill(event);
      mytag = tag + "_GENRECO";
      HFolder(mytag)->fill(event);
    }
    //  }
}
  
void TstarTstarSelectionModule::book_TTreeVars(uhh2::Context& ctx, vector<string> vars){
  for(const auto & var : vars){
    string  myvar = "TstarTstarReco_" + var;
    Event::Handle<float> dummy;
    h_MLinput_values.push_back(dummy);
    h_MLinput_names.push_back(myvar);
    h_MLinput_values[h_MLinput_values.size()-1] = ctx.declare_event_output<float>(myvar);
  }
}

void TstarTstarSelectionModule::fill_TTreeVars(uhh2::Event& event,string var_name, float var_value){ //fill the variable with calculated value
  for(unsigned int iname = 0; iname<h_MLinput_names.size(); iname++){
    if(h_MLinput_names[iname]==var_name)  event.set(h_MLinput_values[iname],var_value);
  }
}

void TstarTstarSelectionModule::fill_MLinput(uhh2::Event & event){
  ReconstructionTstarHypothesis hyp_Tstar = event.get(h_recohyp_tstartstar);  
  ReconstructionHypothesis hyp_ttbar = hyp_Tstar.ttbar_hyp();

  //lepton
  float lep_pt = hyp_ttbar.lepton().pt();
  fill_TTreeVars(event, "TstarTstarReco_lepton_pt", lep_pt);
  float lep_eta = hyp_ttbar.lepton().eta();
  fill_TTreeVars(event, "TstarTstarReco_lepton_eta", lep_eta);
  float lep_phi = hyp_ttbar.lepton().phi();
  fill_TTreeVars(event, "TstarTstarReco_lepton_phi", lep_phi);

  //dR
  float dR_ttbar = deltaR(hyp_ttbar.toplep_v4(),hyp_ttbar.tophad_v4());
  fill_TTreeVars(event, "TstarTstarReco_dR_ttbar", dR_ttbar);
  float dR_gluon_toplep = deltaR(hyp_Tstar.gluon1_v4(),hyp_ttbar.toplep_v4());
  fill_TTreeVars(event, "TstarTstarReco_dR_gluon_toplep", dR_gluon_toplep);
  float dR_gluon_tophad = deltaR(hyp_Tstar.gluon1_v4(),hyp_ttbar.tophad_v4());
  fill_TTreeVars(event, "TstarTstarReco_dR_gluon_tophad", dR_gluon_tophad);
  float dR_tstartstar = deltaR(hyp_Tstar.tstarlep_v4(),hyp_Tstar.tstarhad_v4());
  fill_TTreeVars(event, "TstarTstarReco_dR_tstartstar", dR_tstartstar);

  //gluon
  float gluon_pt = hyp_Tstar.gluon1_v4().pt();
  float gluon_eta = hyp_Tstar.gluon1_v4().eta();
  float gluon_phi = hyp_Tstar.gluon1_v4().phi();
  fill_TTreeVars(event, "TstarTstarReco_gluon_pt", gluon_pt);
  fill_TTreeVars(event, "TstarTstarReco_gluon_eta", gluon_eta);
  fill_TTreeVars(event, "TstarTstarReco_gluon_phi", gluon_phi);

  //Tstar leptonic and hadronic
  float Tstarlep_pt = hyp_Tstar.tstarlep_v4().pt();
  float Tstarlep_eta = hyp_Tstar.tstarlep_v4().eta();
  float Tstarlep_phi = hyp_Tstar.tstarlep_v4().phi();
  float Tstarhad_pt = hyp_Tstar.tstarhad_v4().pt();
  float Tstarhad_eta = hyp_Tstar.tstarhad_v4().eta();
  float Tstarhad_phi = hyp_Tstar.tstarhad_v4().phi();

  fill_TTreeVars(event, "TstarTstarReco_Tstarlep_pt", Tstarlep_pt);
  fill_TTreeVars(event, "TstarTstarReco_Tstarlep_eta", Tstarlep_eta);
  fill_TTreeVars(event, "TstarTstarReco_Tstarlep_phi", Tstarlep_phi);
  fill_TTreeVars(event, "TstarTstarReco_Tstarhad_pt", Tstarhad_pt);
  fill_TTreeVars(event, "TstarTstarReco_Tstarhad_eta", Tstarhad_eta);
  fill_TTreeVars(event, "TstarTstarReco_Tstarhad_phi", Tstarhad_phi);


}


TstarTstarSelectionModule::TstarTstarSelectionModule(Context & ctx){
    
  is_tgtg = false; is_tgtgamma = false;
  if(ctx.get("channel") == "tgtg") is_tgtg = true;
  if(ctx.get("channel") == "tgtgamma") is_tgtgamma = true;
  ctx.undeclare_all_event_output();//store only what we really want later

  h_weight = ctx.declare_event_output<float>("eventweight");

  if(debug) {
    cout << "Hello World from TstarTstarSelectionModule!" << endl;  
    
    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }
    
   }

    // 1. setup common selection  module -> skipped here, because done at preselection
    // 1b. set up lumi rewitghting
    LumiWeight_module.reset(new MCLumiWeight(ctx));
    

    // 2. set up selections
    //Trigger selection
    //

    //HT triggers
    triggerHT1_sel.reset(new TriggerSelection("HLT_HT430to450_v*"));
    triggerHT2_sel.reset(new TriggerSelection("HLT_HT450to470_v*"));
    triggerHT3_sel.reset(new TriggerSelection("HLT_HT470to500_v*"));
    triggerHT4_sel.reset(new TriggerSelection("HLT_HT500to550_v*"));
    triggerHT5_sel.reset(new TriggerSelection("HLT_HT550to650_v*"));
    triggerHT6_sel.reset(new TriggerSelection("HLT_HT650_v*"));

    //PF HT trigger
    triggerPFHT_sel.reset(new TriggerSelection("HLT_PFHT900_v*"));

    //SingleJet trigger
    triggerSingleJet450_sel.reset(new TriggerSelection("HLT_PFJet450_v*"));

    //Lepton trigger
    triggerSingleLeptonEle1_sel.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
    triggerSingleLeptonEle2_sel.reset(new TriggerSelection("HLT_Ele25_eta2p1_WPTight_Gsf_v*"));
    triggerSingleLeptonEle3_sel.reset(new TriggerSelection("HLT_Ele32_eta2p1_WPTight_Gsf_v*"));


    triggerSingleLeptonMu1_sel.reset(new TriggerSelection("HLT_Mu50_v*"));
    triggerSingleLeptonMu2_sel.reset(new TriggerSelection("HLT_Mu55_v*"));
    triggerSingleLeptonMu3_sel.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    triggerSingleLeptonMu4_sel.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));

    //MET selection
    met_sel.reset(new METCut  (50.,1e6));

    //ST selection
    st_sel.reset(new STCut  (500.,1e6));

    //Ak8jet selection
    topjet_selection.reset(new NTopJetSelection(1, -1, TopJetId(PtEtaCut(100, 2.1))));

    //Match to semileptonic ttbar
    TTbarSemiLepMatchable_selection.reset(new TTbarSemiLepMatchableSelection());// for x-checks

    // 3. Set up Tstar-Tstar reconstruction
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
    h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

    const std::string ttbar_hyps_label("TTbarReconstruction");
    const std::string ttbar_chi2_label("Chi2");
    reco_primlep.reset(new PrimaryLepton(ctx));
    if(is_tgtg){
      ttbar_reco.reset(new HighMassSkipJetsTTbarReconstruction(ctx, NeutrinoReconstruction,ttbar_hyps_label,0));
      //      h_ttbar_hyps = ctx.get_handle<std::vector<ReconstructionHypothesis>>(ttbar_hyps_label);
      h_is_ttbar_reconstructed = ctx.get_handle< bool >("is_ttbar_reconstructed_chi2");
      h_recohyp = ctx.declare_event_output<ReconstructionHypothesis>(ttbar_hyps_label+"_best");

      const std::string tstartstar_hyps_label("TstarTstar_tgtg");
      h_tstartstar_hyps = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>(tstartstar_hyps_label);
      h_recohyp_tstartstar = ctx.declare_event_output<ReconstructionTstarHypothesis>(tstartstar_hyps_label+"_best");
    }
    if(is_tgtgamma){
      ttbar_reco.reset(new HighMassSkipJetsTTbarReconstruction(ctx, NeutrinoReconstruction,ttbar_hyps_label,1));
      //      h_ttbar_hyps = ctx.get_handle<std::vector<ReconstructionHypothesis>>(ttbar_hyps_label);
      h_is_ttbar_reconstructed = ctx.get_handle< bool >("is_ttbar_reconstructed_chi2");
      h_recohyp = ctx.declare_event_output<ReconstructionHypothesis>(ttbar_hyps_label+"_best");

      const std::string tstartstar_hyps_label("TstarTstar_tgtgamma");
      h_tstartstar_hyps = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>(tstartstar_hyps_label);
      h_recohyp_tstartstar = ctx.declare_event_output<ReconstructionTstarHypothesis>(tstartstar_hyps_label+"_best");
    }

    ttbar_discriminator.reset(new ttbarChi2Discriminator(ctx));
    ttbar_CorrectMatchDiscriminator.reset(new CorrectMatchDiscriminator(ctx,ttbar_hyps_label));

    TstarTstar_reco.reset(new TstarTstar_Reconstruction(ctx));
    TstarTstar_tgluon_tgluon_reco.reset(new TstarTstar_tgluon_tgluon_Reconstruction(ctx));

    // 4. Set up Hists
    vector<string> histogram_tags = {"PreSelection","PreSelection_mu","PreSelection_ele","AK8sel","AK8sel_mu","AK8sel_ele",
				     "MET","MET_mu","MET_ele","ST","ST_mu","ST_ele","triggerSingleLeptonMu",
				     "triggerSingleLeptonEle","triggerSingleJet_mu","triggerSingleJet_ele",
				     "triggerHT_mu","triggerHT_ele","triggerPFHT_mu","triggerPFHT_ele"};
    book_histograms(ctx, histogram_tags);
    h_RecoPlots_After_TstarTstar_tgtg.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_After_TstarTstar_tgtg"));
    h_RecoPlots_After_TstarTstar_tgtg_ttbarsemilep.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_After_TstarTstar_tgtg_SemiLepTTBarMatch"));

    // 5. Set up variables for ML
    vector<string> vars = {"lepton_pt","lepton_eta","lepton_phi",
			   "dR_ttbar","dR_gluon_toplep","dR_gluon_tophad","dR_tstartstar",
			   "gluon_pt","gluon_eta","gluon_phi",
			   "Tstarlep_pt","Tstarlep_eta","Tstarlep_phi",
			   "Tstarhad_pt","Tstarhad_eta","Tstarhad_phi"};
    book_TTreeVars(ctx, vars);
}


bool TstarTstarSelectionModule::process(Event & event) {

  event.set(h_is_ttbar_reconstructed, false);
  event.set(h_recohyp, ReconstructionHypothesis());
  event.set(h_recohyp_tstartstar,ReconstructionTstarHypothesis());
  
  if(debug)   
    cout << "TstarTstarSelectionModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
  LumiWeight_module->process(event);
  const bool pass_ttbarsemilep = TTbarSemiLepMatchable_selection->passes(event);
  fill_histograms(event, "PreSelection", pass_ttbarsemilep);
  if(event.muons->size() == 1) fill_histograms(event, "PreSelection_mu", pass_ttbarsemilep);
  if(event.electrons->size() == 1) fill_histograms(event, "PreSelection_ele", pass_ttbarsemilep);
 if(debug)   
   cout << "pass_ttbarsemilep" <<endl;
 if(event.jets->size()<1) return false;//FixMe: why this is happening?


 bool pass_ak8 = topjet_selection->passes(event);
 if(!pass_ak8) return false;

 fill_histograms(event, "AK8sel", pass_ttbarsemilep);
 if(event.muons->size() == 1) fill_histograms(event, "AK8sel_mu", pass_ttbarsemilep);
 if(event.electrons->size() == 1) fill_histograms(event, "AK8sel_ele", pass_ttbarsemilep);


 bool pass_MET =  met_sel->passes(event);
 if(!pass_MET) return false;

 fill_histograms(event, "MET", pass_ttbarsemilep);
 if(event.muons->size() == 1) fill_histograms(event, "MET_mu", pass_ttbarsemilep);
 if(event.electrons->size() == 1) fill_histograms(event, "MET_ele", pass_ttbarsemilep);

 bool pass_ST =  st_sel->passes(event);
 if(!pass_ST) return false;
 fill_histograms(event, "ST", pass_ttbarsemilep);
 if(event.muons->size() == 1) fill_histograms(event, "ST_mu", pass_ttbarsemilep);
 if(event.electrons->size() == 1) fill_histograms(event, "ST_ele", pass_ttbarsemilep);

 //Fill ttgen object for correct matching check, etc
  ttgenprod->process(event);

  bool pass_trigger_SingleJet = (triggerSingleJet450_sel->passes(event) && event.jets->at(0).pt()>500);
  if(pass_trigger_SingleJet){
    if(event.muons->size() == 1) fill_histograms(event, "triggerSingleJet_mu", pass_ttbarsemilep);
    if(event.electrons->size() == 1) fill_histograms(event, "triggerSingleJet_ele", pass_ttbarsemilep);
  }
 if(debug)   
   cout << "pass_trigger_SingleJet" <<endl;

  bool pass_trigger_SingleMu = (triggerSingleLeptonMu1_sel->passes(event) || triggerSingleLeptonMu2_sel->passes(event) 
				|| triggerSingleLeptonMu3_sel->passes(event) || triggerSingleLeptonMu4_sel->passes(event));

  if(pass_trigger_SingleMu){
    if((event.muons->size() == 1) && (event.muons->at(0).pt()>60))
      fill_histograms(event, "triggerSingleLeptonMu", pass_ttbarsemilep); //FixMe: each Muon trigger should have its own threshold
  }
 if(debug)   
   cout << "pass_trigger_SingleMu" <<endl;

 bool pass_trigger_SingleEle = (triggerSingleLeptonEle1_sel->passes(event) || triggerSingleLeptonEle2_sel->passes(event) || triggerSingleLeptonEle3_sel->passes(event));
 if(pass_trigger_SingleEle){
   if((event.electrons->size() == 1) && (event.electrons->at(0).pt()>120)) 
     fill_histograms(event, "triggerSingleLeptonEle", pass_ttbarsemilep); //FixMe: each Electron trigger should have its own threshold
 }
 if(debug)   
   cout << "pass_trigger_SingleEle" <<endl;

 double st_jets = 0.;
 std::vector<Jet>* jets = event.jets;
 for(const auto & jet : *jets) st_jets += jet.pt();

 bool pass_trigegr_HT = (triggerHT1_sel->passes(event) || triggerHT2_sel->passes(event) || triggerHT3_sel->passes(event) 
    || triggerHT4_sel->passes(event) || triggerHT5_sel->passes(event) || triggerHT6_sel->passes(event));
  if(pass_trigegr_HT){
    if((event.muons->size() == 1) && st_jets>650) fill_histograms(event, "triggerHT_mu", pass_ttbarsemilep); //FixME: each HT trigger should have different st_jets threshold
    if((event.electrons->size() == 1) && st_jets>650) fill_histograms(event, "triggerHT_ele", pass_ttbarsemilep);
  }
 if(debug)   
   cout << "pass_trigegr_HT" <<endl;

  bool pass_trigegr_PFHT = triggerPFHT_sel->passes(event);
  if(pass_trigegr_PFHT){
    if((event.muons->size() == 1) && st_jets>900) fill_histograms(event, "triggerPFHT_mu", pass_ttbarsemilep); //FixME: PFHT trigger might have different st_jets threshold
    if((event.electrons->size() == 1) && st_jets>900) fill_histograms(event, "triggerPFHT_ele", pass_ttbarsemilep);
  }
 if(debug)   
   cout << "pass_trigegr_PFHT"<<endl;


  reco_primlep->process(event);//set "primary lepton"
  if(debug) {cout << "Starting ttbar reconstruction... ";}\
  ttbar_reco->process(event);//reconstruct ttbar
  ttbar_CorrectMatchDiscriminator->process(event);//find matched to ttbar gen hypothesis
  if(debug) {cout << "Finished. Finding best Hypothesis..."<< endl;}   
  ttbar_discriminator->process(event);

  //make sure we fill all values with smth
  for(unsigned int iname = 0; iname<h_MLinput_names.size(); iname++)
    event.set(h_MLinput_values[iname],0);

  ReconstructionHypothesis hyp = event.get(h_recohyp);  
  if(debug) {cout << "Start TstarTstar reconstruction ..."<< endl;}
  bool pass_tgluon_tgluon_reco = false;
  if(event.get(h_is_ttbar_reconstructed)){
    if(is_tgtgamma){ // Tstar+Tstar -> t+g + t+gamma //FixME: the code is broken and most probably won't work. Sorry!
      cout<<"Nothing is done for tgtgamma at the moment"<<endl;
    }
    if(is_tgtg){ // Tstar+Tstar -> t+g + t+g
      pass_tgluon_tgluon_reco = TstarTstar_tgluon_tgluon_reco->process(event);
      if(pass_tgluon_tgluon_reco){
	h_RecoPlots_After_TstarTstar_tgtg->fill(event);
	fill_MLinput(event);
      }
      if(pass_tgluon_tgluon_reco && pass_ttbarsemilep) 
	h_RecoPlots_After_TstarTstar_tgtg_ttbarsemilep->fill(event);
    }
    // void TstarTstarSelectionModule::fill_TTreeVars(string var_name, float var_value){ 
  }
  else{
    if(debug){cout << "Event has no best hypothesis!" << endl;}
  }

  event.set(h_weight,event.weight);
  //  return true;
  return pass_tgluon_tgluon_reco;//store only events, which pass tgtg reconstruction
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarSelectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarSelectionModule)

//}
