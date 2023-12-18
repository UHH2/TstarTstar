#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/TstarTstar/include/ReconstructionTstarHypothesis.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"
#include "UHH2/common/include/JetIds.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

float inv_mass_2(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }


TstarTstarHists::TstarTstarHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  topjetID = AndId<TopJet>(HOTVRTopTag(), Tau32Groomed(0.56));

  is_MC = ctx.get("dataset_type") == "MC";

  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_ST_AK4 = ctx.get_handle<double>("ST_AK4");
  h_ST_HOTVR = ctx.get_handle<double>("ST_HOTVR");

  h_tstartstar_hyp_gHOTVR = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_Hyp_gHOTVR");
  h_tstartstar_hyp_gAK4 = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_Hyp_gAK4");

  // book all histograms here
  book<TH1F>("weight", "weight", 55, -1, 10);
  book<TH1F>("weight_2", "weight_2", 55, -10, 100);
  book<TH1F>("weight_3", "weight_3", 55, -100, 1000);
  book<TH1F>("weight_4", "weight_4", 55, -1000, 10000);
  book<TH1F>("weight_5", "weight_5", 55, -10000, 100000);

  book<TH1F>("prefiring_weight", "weight", 50, 0, 2);

  h_weight_sfelec_idNominal = ctx.get_handle<float>("weight_sfelec_id");
  h_weight_sfelec_recoNominal = ctx.get_handle<float>("weight_sfelec_reco");
  h_weight_sfelec_triggerNominal = ctx.get_handle<float>("weight_sfelec_trigger");
  book<TH1F>("electron_id_weight", "weight", 50, 0.75, 1.25);
  book<TH1F>("electron_reco_weight", "weight", 50, 0.75, 1.25);
  book<TH1F>("electron_trigger_weight", "weight", 50, 0.75, 1.25);

  // jets
  book<TH1F>("N_jets", "N_{AK4 jets}", 20, 0, 20);
  book<TH1F>("N_jets_btag_loose", "N_{AK4 jets, deepCSV > 0.22}", 20, 0, 20);
  book<TH1F>("N_jets_btag_medium", "N_{AK4 jets, deepCSV > 0.63}", 20, 0, 20);
  book<TH1F>("N_jets_btag_tight", "N_{AK4 jets, deepCSV > 0.90}", 20, 0, 20);
  book<TH1F>("DeepJetscore", "deepJet output", 20, 0, 1);
  book<TH1F>("DeepJetscore_firstJet", "deepJet output first jet", 20, 0, 1);

  book<TH1F>("N_HOTVRjets", "N_{HOTVR jets}", 20, 0, 20);
  book<TH1F>("N_toptagged_HOTVRjets", "N_{toptagged HOTVR jets}", 20, 0, 20);
  book<TH1F>("N_PU", "N_{PU}", 100, 0, 100);

  const int n_jet_pt_bins = 15;
  double jet_pt_bins[n_jet_pt_bins] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 750, 1000, 1500};
  book<TH1F>("pt_jet0", "p_{T}^{jet 0}", 100, 10, 2000);
  book<TH1F>("pt_jet1", "p_{T}^{jet 1}", 100, 10, 2000);
  book<TH1F>("pt_jet2", "p_{T}^{jet 2}", 100, 10, 2000);
  book<TH1F>("pt_jet3", "p_{T}^{jet 3}", 100, 10, 2000);
  book<TH1F>("pt_jet0_rebinned", "p_{T}^{jet 0}", n_jet_pt_bins - 1, jet_pt_bins);
  book<TH1F>("pt_jet1_rebinned", "p_{T}^{jet 1}", n_jet_pt_bins - 1, jet_pt_bins);
  book<TH1F>("pt_jet2_rebinned", "p_{T}^{jet 2}", n_jet_pt_bins - 1, jet_pt_bins);
  book<TH1F>("pt_jet3_rebinned", "p_{T}^{jet 3}", n_jet_pt_bins - 1, jet_pt_bins);

  book<TH1F>("eta_jet0", "#eta^{jet 0}", 50, -5.2, 5.2);
  book<TH1F>("eta_jet1", "#eta^{jet 1}", 50, -5.2, 5.2);
  book<TH1F>("eta_jet2", "#eta^{jet 2}", 50, -5.2, 5.2);
  book<TH1F>("eta_jet3", "#eta^{jet 3}", 50, -5.2, 5.2);


  book<TH1F>("pt_bjet0", "p_{T}^{bjet 0}", 50, 10, 1000);
  book<TH1F>("pt_bjet0_rebinned", "p_{T}^{bjet 0}", n_jet_pt_bins - 1, jet_pt_bins);
  book<TH1F>("eta_bjet0", "#eta^{bjet 0}",  50, -5.2, 5.2);

  book<TH1F>("eta_HOTVR0", "#eta^{HOTVR 0}", 50, -5.2, 5.2);
  book<TH1F>("pt_HOTVR0", "p_{T}^{HOTVR 0}", 100, 10, 2000);
  book<TH1F>("pt_HOTVR0_rebinned", "p_{T}^{HOTVR 0}", n_jet_pt_bins - 1, jet_pt_bins);
  book<TH1F>("eta_HOTVR1", "#eta^{HOTVR 1}", 50, -5.2, 5.2);
  book<TH1F>("pt_HOTVR1", "p_{T}^{HOTVR 1}", 100, 10, 2000);
  book<TH1F>("pt_HOTVR1_rebinned", "p_{T}^{HOTVR 1}", n_jet_pt_bins - 1, jet_pt_bins);
  book<TH1F>("eta_HOTVR2", "#eta^{HOTVR 2}", 50, -5.2, 5.2);
  book<TH1F>("pt_HOTVR2", "p_{T}^{HOTVR 2}", 100, 10, 2000);
  book<TH1F>("pt_HOTVR2_rebinned", "p_{T}^{HOTVR 2}", n_jet_pt_bins - 1, jet_pt_bins);
  book<TH1F>("eta_HOTVR3", "#eta^{HOTVR 2}", 50, -5.2, 5.2);
  book<TH1F>("pt_HOTVR3", "p_{T}^{HOTVR 2}", 100, 10, 2000);
  book<TH1F>("pt_HOTVR3_rebinned", "p_{T}^{HOTVR 3}", n_jet_pt_bins - 1, jet_pt_bins);

  book<TH1F>("R_fat_jet", "R_{fat jet} (only valid for HOTVR)", 30, 0, 6);

  book<TH2D>("pt_mu_pt_ak4jet0", ";p_{T}^{#mu}; p_{T}^{AK4 jet 0}", 60, 10, 1100, 70, 10, 2500);
  book<TH2D>("pt_ele_pt_ak4jet0", ";p_{T}^{ele}; p_{T}^{AK4 jet 0}", 60, 10, 1100, 70, 10, 2500);
  book<TH2D>("pt_mu_pt_HOTVR0", ";p_{T}^{#mu}; p_{T}^{HOTVR jet 0}", 60, 10, 1100, 70, 10, 2500);
  book<TH2D>("pt_ele_pt_HOTVR0", ";p_{T}^{ele}; p_{T}^{HOTVR jet 0}", 60, 10, 1100, 70, 10, 2500);

  book<TH1F>("tau32_HOTVR0", "HOTVR 0 #tau_{32}", 20, 0, 1);
  book<TH1F>("jetmass_HOTVR0", "m_{HOTVR 0}", 25, 0, 500);

  // leptons
  book<TH1F>("N_mu", "N^{#mu}", 10, 0, 10);
  book<TH1F>("pt_mu", "p_{T}^{#mu} [GeV/c]", 50, 0, 1000);
  book<TH1F>("pt_mu_low", "p_{T}^{#mu} [GeV/c]", 50, 0, 100);
  book<TH1F>("pt_mu_rebinned", "p_{T}^{#mu} [GeV/c]", n_jet_pt_bins - 1, jet_pt_bins);
  book<TH1F>("eta_mu", "#eta^{#mu}", 50, -5.2, 5.2);
  book<TH1F>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);

  book<TH1F>("N_ele", "N^{e}", 10, 0, 10);
  book<TH1F>("pt_ele", "p_{T}^{ele} [GeV/c]", 50, 0, 1000);
  book<TH1F>("pt_ele_rebinned", "p_{T}^{ele} [GeV/c]", n_jet_pt_bins - 1, jet_pt_bins);
  book<TH1F>("eta_ele", "#eta^{ele}", 50, -5.2, 5.2);

  const int n_ele_pt_bins = 12;
  const int n_ele_eta_bins = 11;
  double ele_pt_bins[n_ele_pt_bins] = {55, 75, 95, 115, 135, 155, 175, 200, 250, 300, 400, 500};
  double ele_eta_bins[n_ele_eta_bins] = {-2.5, -2.0, -1.566, -1.444, -0.8, 0.0, 0.8, 1.444, 1.566, 2.0, 2.5};
  book<TH2D>("ele_pt_eta", ";p_{T}^{ele}; #eta^{ele}", n_ele_pt_bins-1, ele_pt_bins, n_ele_eta_bins-1, ele_eta_bins);

  book<TH1F>("N_photon", "N^{#gamma}", 15, 0, 15);
  book<TH1F>("pt_photon", "p_{T}^{#gamma} [GeV/c]", 50, 10, 2000);
  book<TH1F>("eta_photon", "#eta^{#gamma}", 50, -5.2, 5.2);

  book<TH1F>("pt_photon_1", "p_{T}^{leading #gamma} [GeV/c]", 50, 10, 2000);
  book<TH1F>("pt_photon_2", "p_{T}^{second #gamma} [GeV/c]", 50, 10, 2000);

  // primary vertices
  book<TH1F>("N_pv", "N^{PV}", 50, 0, 50);

  // MET
  const int n_met_pt_bins = 15;
  double met_pt_bins[n_met_pt_bins] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 750, 1000, 1500};
  book<TH1F>("pt_MET", "missing E_{T} [GeV]", 100, 0, 2000);
  book<TH1F>("pt_MET_rebinned", "missing E_{T} [GeV]", n_met_pt_bins - 1, met_pt_bins);
  book<TH1F>("phi_MET", "missing E_{T} #phi", 40, -4, 4);

  const int n_met_pxy_bins = 11;
  double met_pxy_bins[n_met_pxy_bins] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 500, 1000};
  book<TH1F>("px_MET", "missing E_{T} p_{x}", 25, 0, 1000);
  book<TH1F>("py_MET", "missing E_{T} p_{y}", 25, -0, 1000);
  book<TH1F>("px_MET_rebinned", "missing E_{T} p_{x}", n_met_pxy_bins - 1, met_pxy_bins);
  book<TH1F>("py_MET_rebinned", "missing E_{T} p_{y}", n_met_pxy_bins - 1, met_pxy_bins);

  // raw mett (lecker)
  book<TH1F>("pt_rawMET", "raw missing E_{T} [GeV]", 100, 0, 2000);
  book<TH1F>("pt_rawMET_rebinned", "raw missing E_{T} [GeV]", n_met_pt_bins - 1, met_pt_bins);
  book<TH1F>("phi_rawMET", "raw missing E_{T} #phi", 40, -4, 4);

  // CHS MET (assuming I use puppi...) TODO

  // binnings for the folliwng

  const int nbinsNjet = 8;
  double binsNjet[nbinsNjet] {0, 4, 5, 6, 7, 8, 10, 20};

  //const int nbins = 24;
  //double bins[nbins] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 4000, 6000};
  const int nbins = 34;
  double bins[nbins] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3250, 4000, 6000};

  // HT AK4
  // // // TODO rename all these to _AK4 at some point
  // // // note that the btag yield corrections need this to be named this way for now!
  book<TH1F>("pt_HT", "H_{T} [GeV]", 40, 0, 4000);
  book<TH2D>("pt_HT_N_jet", ";H_{T}; N(jet)", 40, 0, 4000, 20, 0, 20);
  book<TH2D>("pt_HT_N_jet_rebinned", ";H_{T}; N(jet)", nbins-1, bins, nbinsNjet-1, binsNjet);
  book<TH1F>("pt_HT_rebinned", "H_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_AK4", "S_{T} (AK4) [GeV]", 60, 0, 6000);
  book<TH1F>("pt_ST_AK4_rebinned", "S_{T} (AK4) [GeV]", nbins-1, bins);

  // HT & ST HOTVR
  book<TH1F>("pt_HT_HOTVR", "H_{T} [GeV]", 60, 0, 6000);
  book<TH1F>("pt_HT_HOTVR_rebinned", "H_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_HTlep", "H_{T} + p_{T} (#ell) [GeV]", 60, 0, 6000);
  book<TH1F>("pt_ST", "S_{T} [GeV]", 60, 0, 6000);
  book<TH1F>("pt_HTlep_rebinned", "H_{T} + p_{T} (#ell) [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_rebinned", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_HT_HOTVR_rebinned_nohighjets", "H_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_HT_HOTVR_rebinned_only3jets", "H_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_HTlep_rebinned_nohighjets", "H_{T} + p_{T} (#ell) [GeV]", nbins-1, bins);
  book<TH1F>("pt_HTlep_rebinned_only3jets", "H_{T} + p_{T} (#ell) [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_rebinned_nohighjets", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_rebinned_only3jets", "S_{T} [GeV]", nbins-1, bins);

  // pt asymmetries
  book<TH1F>("pt_asym12", "#Delta p_{T} (jet 1, jet 2) [GeV]", 20, 0, 2000);
  book<TH1F>("pt_asym13", "#Delta p_{T} (jet 1, jet 3) [GeV]", 20, 0, 2000);
  book<TH1F>("pt_asym14", "#Delta p_{T} (jet 1, jet 4) [GeV]", 20, 0, 2000);
  book<TH1F>("pt_asym23", "#Delta p_{T} (jet 2, jet 3) [GeV]", 20, 0, 2000);
  book<TH1F>("pt_asym24", "#Delta p_{T} (jet 2, jet 4) [GeV]", 20, 0, 2000);
  book<TH1F>("pt_asym34", "#Delta p_{T} (jet 3, jet 4) [GeV]", 20, 0, 2000);
  book<TH1F>("pt_asym12_over_ST", "#Delta p_{T} (jet 1, jet 2) / S_{T} [GeV]", 15, 0, 1.5);


  // deltaR observables
  book<TH1F>("dR_lepton_closestJet", "#DeltaR (lepton, closestJet)", 30, 0, 6);
  book<TH1F>("dR_lepton_closestJet_fine", "#DeltaR (lepton, closestJet)", 40, 0, 2);
  book<TH1F>("dR_lepton_closestHOTVR", "#DeltaR (lepton, closest HOTVR)", 30, 0, 6);
  book<TH1F>("dR_lepton_closestHOTVR_fine", "#DeltaR (lepton, closest HOTVR)", 40, 0, 2);
  book<TH1F>("dR_fatjet0_fatjet1", "#DeltaR (hotvrjet 0, hotvrjet 1)", 20, 0, 6);
  book<TH1F>("dR_lepton_fatjet0", "#DeltaR (lepton, hotvrjet 0)", 20, 0, 6);
  book<TH1F>("dR_lepton_fatjet1", "#DeltaR (lepton, hotvrjet 1)", 20, 0, 6);
  book<TH1F>("dR_ttagjet_gluonjet", "#DeltaR (ttagjet, gluonjet)", 20, 0, 6);

  // deltaPhi observables
  book<TH1F>("dphi_fatjet0_fatjet1", "#Delta#phi (hotvrjet 0, hotvrjet 1)", 21, 0, 7);
  book<TH1F>("dphi_lepton_fatjet0", "#Delta#phi (lepton, hotvrjet 0)", 21, 0, 7);
  book<TH1F>("dphi_lepton_fatjet1", "#Delta#phi (lepton, hotvrjet 1)", 21, 0, 7);
  book<TH1F>("dphi_ttagjet_gluonjet", "#Delta#phi (ttagjet, gluonjet)", 21, 0, 7);

  // some invariant masses
  book<TH1F>("invmass_lep_b", "M_{\ell, b}", 25, 0, 500);
  book<TH1F>("invmass_lep_MET", "M_{\ell, MET}", 25, 0, 500);
  book<TH1F>("invmass_b_MET", "M_{b, MET}", 25, 0, 500);

  // reconstructred Tstar masses
  // HOTVR approach
  book<TH1F>("chi2_gHOTVR", "#chi^{2} gHOTVR", 25, 0, 100);
  book<TH1F>("M_Tstar_gHOTVR_had", "M_{T^{*}} gHOTVR had", 30, 0, 3000);
  book<TH1F>("M_Tstar_gHOTVR_lep", "M_{T^{*}} gHOTVR lep", 30, 0, 3000);
  book<TH1F>("M_Tstar_gHOTVR_uncorr", "M_{T^{*}} gHOTVR uncorr", 30, 0, 3000);
  book<TH1F>("M_Tstar_gHOTVR", "M_{T^{*}} gHOTVR", 30, 0, 3000);
  book<TH1F>("M_top_gHOTVR_had", "M_{t} gHOTVR had", 30, 0, 500);
  book<TH1F>("M_top_gHOTVR_lep", "M_{t} gHOTVR lep", 30, 0, 500);

  // AK4 jet approach
  book<TH1F>("chi2_gAK4", "#chi^{2} gAK4", 25, 0, 100);
  book<TH1F>("M_Tstar_gAK4_had", "M_{T^{*}} gAK4 had", 30, 0, 3000);
  book<TH1F>("M_Tstar_gAK4_lep", "M_{T^{*}} gAK4 lep", 30, 0, 3000);
  book<TH1F>("M_Tstar_gAK4_uncorr", "M_{T^{*}} gAK4 uncorr", 30, 0, 3000);
  book<TH1F>("M_Tstar_gAK4", "M_{T^{*}} gAK4", 30, 0, 3000);
  book<TH1F>("M_top_gAK4_had", "M_{t} gAK4 had", 30, 0, 500);
  book<TH1F>("M_top_gAK4_lep", "M_{t} gAK4 lep", 30, 0, 500);

  book<TH1F>("mjj", "m_{j j}", 30, 0, 3000);

  book<TH2D>("HEMcheck", ";#eta_{AK4}; #phi_{AK4}", 40, -4, 4, 40, -4, 4);


}


void TstarTstarHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  bool debug = false;

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  if(debug) cout << "Starting Tstar Hists." << endl;

  hist("weight")->Fill(weight);
  hist("weight_2")->Fill(weight);
  hist("weight_3")->Fill(weight);
  hist("weight_4")->Fill(weight);
  hist("weight_5")->Fill(weight);

  if(is_MC) {
    try {
      hist("prefiring_weight")->Fill(event.prefiringWeight, weight);
      hist("electron_id_weight")->Fill(event.get(h_weight_sfelec_idNominal), weight);
      hist("electron_reco_weight")->Fill(event.get(h_weight_sfelec_recoNominal), weight);
      hist("electron_trigger_weight")->Fill(event.get(h_weight_sfelec_triggerNominal), weight);
    } catch (...) {
      if(debug) std::cout << "electron weights not filled as not set yet" << std::endl;
    }
  }

  std::vector<Jet>* jets = event.jets;
  int Njets = jets->size();
  hist("N_jets")->Fill(Njets, weight);

  // b-tag
  BTag bJetID_loose = BTag(BTag::algo::DEEPJET, BTag::wp::WP_LOOSE);
  BTag bJetID_medium = BTag(BTag::algo::DEEPJET, BTag::wp::WP_MEDIUM);
  BTag bJetID_tight = BTag(BTag::algo::DEEPJET, BTag::wp::WP_TIGHT);

  int N_jets_btag_loose = 0;
  int N_jets_btag_medium = 0;
  int N_jets_btag_tight = 0;
  Jet bjet0;
  double highest_btag = -1;
  for(const auto & jet : *event.jets) {
    if(bJetID_loose(jet, event)) N_jets_btag_loose++;
    if(bJetID_medium(jet, event)) N_jets_btag_medium++;
    if(bJetID_tight(jet, event)) N_jets_btag_tight++;
    hist("DeepJetscore")->Fill(jet.btag_DeepJet(), weight);
    if (jet.btag_DeepJet() > highest_btag) {
      highest_btag = jet.btag_DeepJet();
      bjet0 = jet; 
    }
  }
  // fill bjet histograms
  if(Njets>=1){
    hist("pt_bjet0")->Fill(jets->at(0).pt(), weight);
    hist("pt_bjet0_rebinned")->Fill(jets->at(0).pt(), weight);
  }

  if(event.jets->size()>0) hist("DeepJetscore_firstJet")->Fill(event.jets->at(0).btag_DeepJet(), weight);

  hist("N_jets_btag_loose")->Fill(N_jets_btag_loose, weight);
  hist("N_jets_btag_medium")->Fill(N_jets_btag_medium, weight);
  hist("N_jets_btag_tight")->Fill(N_jets_btag_tight, weight);
  if(!event.isRealData)  hist("N_PU")->Fill(event.genInfo->pileup_TrueNumInteractions(), weight);

  if(Njets>=1){
    hist("pt_jet0")->Fill(jets->at(0).pt(), weight);
    hist("pt_jet0_rebinned")->Fill(jets->at(0).pt(), weight);
    hist("eta_jet0")->Fill(jets->at(0).eta(), weight);
  }
  if(Njets>=2){
    hist("eta_jet1")->Fill(jets->at(1).eta(), weight);
    hist("pt_jet1_rebinned")->Fill(jets->at(1).pt(), weight);
    hist("pt_jet1")->Fill(jets->at(1).pt(), weight);
  }
  if(Njets>=3){
    hist("eta_jet2")->Fill(jets->at(2).eta(), weight);
    hist("pt_jet2_rebinned")->Fill(jets->at(2).pt(), weight);
    hist("pt_jet2")->Fill(jets->at(2).pt(), weight);
  }
  if(Njets>=4){
    hist("eta_jet3")->Fill(jets->at(3).eta(), weight);
    hist("pt_jet3_rebinned")->Fill(jets->at(3).pt(), weight);
    hist("pt_jet3")->Fill(jets->at(3).pt(), weight);
  }

  if(event.topjets->size()>0) {
    hist("tau32_HOTVR0")->Fill(event.topjets->at(0).tau3_groomed()/event.topjets->at(0).tau2_groomed(), weight);
    hist("jetmass_HOTVR0")->Fill(inv_mass_2(event.topjets->at(0).v4()), weight);
  }

  if(debug) cout << "Finished filling jet observables." << endl;

  int Nmuons = event.muons->size();
  hist("N_mu")->Fill(Nmuons, weight);
  for (const Muon & thismu : *event.muons){
      hist("pt_mu")->Fill(thismu.pt(), weight);
      hist("pt_mu_low")->Fill(thismu.pt(), weight);
      hist("pt_mu_rebinned")->Fill(thismu.pt(), weight);
      hist("eta_mu")->Fill(thismu.eta(), weight);
      hist("reliso_mu")->Fill(thismu.relIso(), weight);
  }

  int Nele = event.electrons->size();
  hist("N_ele")->Fill(Nele, weight);
  for (const Electron & thisele : *event.electrons){
      hist("pt_ele")->Fill(thisele.pt(), weight);
      hist("pt_ele_rebinned")->Fill(thisele.pt(), weight);
      hist("eta_ele")->Fill(thisele.eta(), weight);
      ((TH2D*)hist("ele_pt_eta"))->Fill(thisele.pt(), thisele.eta(), weight);
  }
  if(debug) cout << "Finished filling lepton observables." << endl;

  int Ngamma = event.photons->size();
  hist("N_photon")->Fill(Ngamma, weight);
  int countgamma = 1;
  for (const Photon & thisgamma : *event.photons){
    hist("pt_photon")->Fill(thisgamma.pt(), weight);
    if(countgamma == 1){hist("pt_photon_1")->Fill(thisgamma.pt(), weight);}
    if(countgamma == 2){hist("pt_photon_2")->Fill(thisgamma.pt(), weight);}
    hist("eta_photon")->Fill(thisgamma.eta(), weight);
    countgamma++;
  }
  if(debug) cout << "Finished filling photon observables." << endl;

  int Npvs = event.pvs->size();
  hist("N_pv")->Fill(Npvs, weight);
  if(debug) cout << "Finished filling PV observables." << endl;

  hist("N_HOTVRjets")->Fill(event.topjets->size(), weight);

  int toptaggedjets = 0;
  for(auto & topjet : * event.topjets){
    double Rtopjet;
    if(topjet.pt() == 0) continue;
    Rtopjet = 600/topjet.pt();
    if(Rtopjet > 1.5) Rtopjet = 1.5;
    if(Rtopjet < 0.1) Rtopjet = 0.1;
    if(topjet.pt() > 0) hist("R_fat_jet")->Fill(Rtopjet, weight);
    if(topjetID(topjet, event)) toptaggedjets++;
  }
  hist("N_toptagged_HOTVRjets")->Fill(toptaggedjets, weight);
  if(debug) cout << "Finished filling ttag observables." << endl;

  if(event.topjets->size()>0){
    hist("eta_HOTVR0")->Fill(event.topjets->at(0).eta(), weight);
    hist("pt_HOTVR0")->Fill(event.topjets->at(0).pt(), weight);
    hist("pt_HOTVR0_rebinned")->Fill(event.topjets->at(0).pt(), weight);
  }
  if(event.topjets->size()>1){
    hist("eta_HOTVR1")->Fill(event.topjets->at(1).eta(), weight);
    hist("pt_HOTVR1")->Fill(event.topjets->at(1).pt(), weight);
    hist("pt_HOTVR1_rebinned")->Fill(event.topjets->at(1).pt(), weight);
  }
  if(event.topjets->size()>2){
    hist("eta_HOTVR2")->Fill(event.topjets->at(2).eta(), weight);
    hist("pt_HOTVR2")->Fill(event.topjets->at(2).pt(), weight);
    hist("pt_HOTVR2_rebinned")->Fill(event.topjets->at(2).pt(), weight);
  }
  if(event.topjets->size()>3){
    hist("eta_HOTVR3")->Fill(event.topjets->at(3).eta(), weight);
    hist("pt_HOTVR3")->Fill(event.topjets->at(3).pt(), weight);
    hist("pt_HOTVR3_rebinned")->Fill(event.topjets->at(3).pt(), weight);
  }
  if(debug) cout << "Finished filling HOTVR jet observables." << endl;

  if(Nele>0 && Njets>0) ((TH2D*)hist("pt_ele_pt_ak4jet0"))->Fill(event.electrons->at(0).pt(),jets->at(0).pt(), weight);
  if(Nmuons>0 && Njets>0) ((TH2D*)hist("pt_mu_pt_ak4jet0"))->Fill(event.muons->at(0).pt(),jets->at(0).pt(), weight);
  if(Nele>0 && event.topjets->size()>0) ((TH2D*)hist("pt_ele_pt_HOTVR0"))->Fill(event.electrons->at(0).pt(),event.topjets->at(0).pt(), weight);
  if(Nmuons>0 && event.topjets->size()>0) ((TH2D*)hist("pt_mu_pt_HOTVR0"))->Fill(event.muons->at(0).pt(),event.topjets->at(0).pt(), weight);
  if(debug) cout << "Finished filling pt observables." << endl;

  hist("pt_MET")->Fill(event.met->pt(), weight);
  hist("phi_MET")->Fill(event.met->phi(), weight);
  hist("px_MET")->Fill(event.met->v4().px(), weight);
  hist("py_MET")->Fill(event.met->v4().py(), weight);
  hist("px_MET_rebinned")->Fill(event.met->v4().px(), weight);
  hist("py_MET_rebinned")->Fill(event.met->v4().py(), weight);
  hist("pt_MET_rebinned")->Fill(event.met->pt(), weight);
  if(debug) cout << "Finished filling MET observables." << endl;

  hist("pt_rawMET")->Fill(event.met->uncorr_pt(), weight);
  hist("pt_rawMET_rebinned")->Fill(event.met->uncorr_pt(), weight);
  hist("phi_rawMET")->Fill(event.met->uncorr_phi(), weight);

  // ST & HT block

  // AK-4 based HT
  double ht_ak4 = 0;
  for(const auto & jet : *event.jets) {
    ht_ak4 += jet.pt();
  }

  if(ht_ak4 > 4000) {
    hist("pt_HT")->Fill(3999.9, weight);
    hist("pt_HT_rebinned")->Fill(3999.9, weight);
    ((TH2D*)hist("pt_HT_N_jet"))->Fill(3999.9, event.jets->size(), weight);
    ((TH2D*)hist("pt_HT_N_jet_rebinned"))->Fill(3999.9, event.jets->size(), weight);
  }
  else {
    hist("pt_HT")->Fill(ht_ak4, weight);
    hist("pt_HT_rebinned")->Fill(ht_ak4, weight);
    ((TH2D*)hist("pt_HT_N_jet"))->Fill(ht_ak4, event.jets->size(), weight);
    ((TH2D*)hist("pt_HT_N_jet_rebinned"))->Fill(ht_ak4, event.jets->size(), weight);
  }

  double ht = 0.;
  double ht_nohighjets = 0.;
  double ht_only3jets = 0.;
  int jeti = 0;
  for(const auto & topjet : *event.topjets) {
    ht += topjet.pt();
    if (topjet.pt() < 1500) ht_nohighjets += topjet.pt();
    if (jeti < 3) ht_only3jets += topjet.pt();
    jeti++;
  }
  if(ht > 6000) {
    hist("pt_HT_HOTVR")->Fill(5999.9, weight);
    hist("pt_HT_HOTVR_rebinned")->Fill(5999.9, weight);
  }
  else {
    hist("pt_HT_HOTVR")->Fill(ht, weight);
    hist("pt_HT_HOTVR_rebinned")->Fill(ht, weight);
  }

  if(ht_nohighjets > 6000) {
    hist("pt_HT_HOTVR_rebinned_nohighjets")->Fill(5999.9, weight);
  } else {
    hist("pt_HT_HOTVR_rebinned_nohighjets")->Fill(ht_nohighjets, weight);
  }

  if(ht_only3jets > 6000) {
    hist("pt_HT_HOTVR_rebinned_only3jets")->Fill(5999.9, weight);
  } else {
    hist("pt_HT_HOTVR_rebinned_only3jets")->Fill(ht_only3jets, weight);
  }

  double st = ht;
  double st_nohighjets = ht_nohighjets;
  double st_only3jets = ht_only3jets;
  for(const auto & lepton : *event.electrons) {
    st += lepton.pt();
    st_nohighjets += lepton.pt();
    st_only3jets += lepton.pt();
  }
  for(const auto & lepton : *event.muons) {
    st += lepton.pt();
    st_nohighjets += lepton.pt();
    st_only3jets += lepton.pt();
  }
  if(st > 6000) {
    hist("pt_HTlep_rebinned")->Fill(5999.9, weight);
    hist("pt_HTlep")->Fill(5999.9, weight);
  }
  else{
    hist("pt_HTlep_rebinned")->Fill(st, weight);
    hist("pt_HTlep")->Fill(st, weight);
  }

  if(st_nohighjets > 6000) {
    hist("pt_HTlep_rebinned_nohighjets")->Fill(5999.9, weight);
  } else {
    hist("pt_HTlep_rebinned_nohighjets")->Fill(st_nohighjets, weight);
  }

  if(st_only3jets > 6000) {
    hist("pt_HTlep_rebinned_only3jets")->Fill(5999.9, weight);
  } else {
    hist("pt_HTlep_rebinned_only3jets")->Fill(st_only3jets, weight);
  }

  // getting the AK4-based ST
  try {
    double stak4 = event.get(h_ST_AK4);
    if(stak4 > 6000) hist("pt_ST_AK4")->Fill(3999.9, weight);
    else hist("pt_ST_AK4")->Fill(stak4, weight);
    if(stak4 > 6000) hist("pt_ST_AK4_rebinned")->Fill(5999.9, weight);
    else hist("pt_ST_AK4_rebinned")->Fill(stak4, weight);
  } catch(...) {}

  try {
    double stHOTVR = event.get(h_ST_HOTVR);

    if (!(stHOTVR == st + event.met->pt())) std::cout << "ATTENTION: STs do not agree!!!" << std::endl;

    if(stHOTVR > 6000) hist("pt_ST")->Fill(3999.9, weight);
    else hist("pt_ST")->Fill(stHOTVR, weight);
    if(stHOTVR > 6000) hist("pt_ST_rebinned")->Fill(5999.9, weight);
    else hist("pt_ST_rebinned")->Fill(stHOTVR, weight);

  } catch(...) {}

  if(st_nohighjets + event.met->pt() > 6000) {
    hist("pt_ST_rebinned_nohighjets")->Fill(5999.9, weight);
  } else {
    hist("pt_ST_rebinned_nohighjets")->Fill(st_nohighjets + event.met->pt(), weight);
  }

  if(st_only3jets + event.met->pt() > 6000) {
    hist("pt_ST_rebinned_only3jets")->Fill(5999.9, weight);
  } else {
    hist("pt_ST_rebinned_only3jets")->Fill(st_only3jets + event.met->pt(), weight);
  }

  // p_T asymmetry
  if(event.jets->size() > 3) {
    hist("pt_asym12")->Fill(event.jets->at(0).pt()-event.jets->at(1).pt(), weight);
    hist("pt_asym13")->Fill(event.jets->at(0).pt()-event.jets->at(2).pt(), weight);
    hist("pt_asym14")->Fill(event.jets->at(0).pt()-event.jets->at(3).pt(), weight);
    hist("pt_asym23")->Fill(event.jets->at(1).pt()-event.jets->at(2).pt(), weight);
    hist("pt_asym24")->Fill(event.jets->at(1).pt()-event.jets->at(3).pt(), weight);
    hist("pt_asym34")->Fill(event.jets->at(2).pt()-event.jets->at(3).pt(), weight);
    hist("pt_asym12_over_ST")->Fill((event.jets->at(0).pt()-event.jets->at(1).pt())/st, weight);
  }

  for (const auto & jet : *event.jets){
    ((TH2D*)hist("HEMcheck"))->Fill(jet.eta(),jet.phi(), weight);
  }

  // dR stuff
  FlavorParticle primary_lepton;
  try {
    primary_lepton = event.get(h_primlep);
  } catch(...) {return;}
  double min_deltaR = 999;
  for(auto &jet : *event.jets){
    double cur_deltaR = deltaR(jet, primary_lepton);
    if(cur_deltaR < min_deltaR) min_deltaR = cur_deltaR;
  }
  hist("dR_lepton_closestJet")->Fill(min_deltaR, weight);
  hist("dR_lepton_closestJet_fine")->Fill(min_deltaR, weight);

  double min_deltaR_hotvr = 999;
  for(auto &jet : *event.topjets){
    double cur_deltaR_hotvr = deltaR(jet, primary_lepton);
    if(cur_deltaR_hotvr < min_deltaR_hotvr) min_deltaR_hotvr = cur_deltaR_hotvr;
  }
  hist("dR_lepton_closestHOTVR")->Fill(min_deltaR_hotvr, weight);
  hist("dR_lepton_closestHOTVR_fine")->Fill(min_deltaR_hotvr, weight);

  if(event.topjets->size()>1){
    hist("dR_fatjet0_fatjet1")->Fill(deltaR(event.topjets->at(0), event.topjets->at(1)), weight);
    hist("dphi_fatjet0_fatjet1")->Fill(abs(event.topjets->at(0).phi() - event.topjets->at(1).phi()), weight);
  }
  if(event.topjets->size()>0){
    hist("dR_lepton_fatjet0")->Fill(deltaR(primary_lepton, event.topjets->at(0)), weight);
    hist("dphi_lepton_fatjet0")->Fill(abs(primary_lepton.phi() - event.topjets->at(0).phi()), weight);
  }
  if(event.topjets->size()>1){
    hist("dR_lepton_fatjet1")->Fill(deltaR(primary_lepton, event.topjets->at(1)), weight);
    hist("dphi_lepton_fatjet1")->Fill(abs(primary_lepton.phi() - event.topjets->at(1).phi()), weight);
  }

  bool found_ttagjet = false;
  bool found_gluonjet = false;

  //find ttagged jet
  std::vector<TopJet> ttagjets;
  for(auto & topjet : * event.topjets){
    if(topjetID(topjet, event)){
      ttagjets.push_back(topjet);
      found_ttagjet = true;
      break;
    }
  }
  //match gluon jets
  TopJet gluonjet;
  TopJet ttagjet_best;
  double currentbestdpt = 999;
  double topmass = 173.0;
  for (const auto & ttagjet : ttagjets){
    for(const auto & topjet : * event.topjets){
      if(topjet.pt() == ttagjet.pt()) continue; // dont take same jet twice
      if(topjet.pt() == 0) continue;
      double topjetradius = 600/(topjet.pt());
      if(topjetradius < 0.1) topjetradius = 0.1;
      else if(topjetradius > 1.5) topjetradius = 1.5;
      if(deltaR(primary_lepton, topjet) < topjetradius) continue;
      if(abs(topjet.pt() - (ttagjet.pt() + topmass)) < currentbestdpt){
        currentbestdpt = abs(topjet.pt() - (ttagjet.pt() + topmass));
        gluonjet = topjet;
        ttagjet_best = ttagjet;
        found_gluonjet = true;
      }
    }
  }
  if(found_ttagjet && found_gluonjet){
    hist("dR_ttagjet_gluonjet")->Fill(deltaR(ttagjet_best, gluonjet), weight);
    hist("dphi_ttagjet_gluonjet")->Fill(abs(ttagjet_best.phi() - gluonjet.phi()), weight);
  }

  for(const auto & jet : *event.jets) {
    if(jet.btag_DeepCSV() > 0.2219) {
      hist("invmass_lep_b")->Fill(inv_mass_2(primary_lepton.v4() + jet.v4()), weight);
      hist("invmass_b_MET")->Fill(inv_mass_2(jet.v4() + event.met->v4()), weight);
      break;
    }
  }
  hist("invmass_lep_MET")->Fill(inv_mass_2(primary_lepton.v4() + event.met->v4()), weight);

  // reco plots
  //gHOTVR
  try {
    ReconstructionTstarHypothesis hyp_gHOTVR = event.get(h_tstartstar_hyp_gHOTVR);
    hist("M_Tstar_gHOTVR_uncorr")->Fill(hyp_gHOTVR.tstarlep_v4().pt(), weight/2);
    if(hyp_gHOTVR.tstarhad_v4().pt() != 0) {
      ReconstructionHypothesis ttbarhyp_gHOTVR = hyp_gHOTVR.ttbar_hyp();
      hist("chi2_gHOTVR")->Fill(hyp_gHOTVR.chi2(), weight);
      hist("M_top_gHOTVR_had")->Fill(inv_mass_2(ttbarhyp_gHOTVR.tophad_v4()), weight);
      hist("M_top_gHOTVR_lep")->Fill(inv_mass_2(ttbarhyp_gHOTVR.toplep_v4()), weight);
      hist("M_Tstar_gHOTVR")->Fill(inv_mass_2(hyp_gHOTVR.tstarlep_v4()), weight/2);
      hist("M_Tstar_gHOTVR_lep")->Fill(inv_mass_2(hyp_gHOTVR.tstarlep_v4()), weight);
      hist("M_Tstar_gHOTVR")->Fill(inv_mass_2(hyp_gHOTVR.tstarhad_v4()), weight/2);
      hist("M_Tstar_gHOTVR_had")->Fill(inv_mass_2(hyp_gHOTVR.tstarhad_v4()), weight);
    }
  } catch(...) {}

  //gAK4
  try {
    ReconstructionTstarHypothesis hyp_gAK4 = event.get(h_tstartstar_hyp_gAK4);
    hist("M_Tstar_gAK4_uncorr")->Fill(hyp_gAK4.tstarlep_v4().pt(), weight/2);
    if(hyp_gAK4.tstarhad_v4().pt() != 0) {
      ReconstructionHypothesis ttbarhyp_gAK4 = hyp_gAK4.ttbar_hyp();
      hist("chi2_gAK4")->Fill(hyp_gAK4.chi2(), weight);
      hist("M_top_gAK4_had")->Fill(inv_mass_2(ttbarhyp_gAK4.tophad_v4()), weight);
      hist("M_top_gAK4_lep")->Fill(inv_mass_2(ttbarhyp_gAK4.toplep_v4()), weight);
      hist("M_Tstar_gAK4")->Fill(inv_mass_2(hyp_gAK4.tstarlep_v4()), weight/2);
      hist("M_Tstar_gAK4_lep")->Fill(inv_mass_2(hyp_gAK4.tstarlep_v4()), weight);
      hist("M_Tstar_gAK4")->Fill(inv_mass_2(hyp_gAK4.tstarhad_v4()), weight/2);
      hist("M_Tstar_gAK4_had")->Fill(inv_mass_2(hyp_gAK4.tstarhad_v4()), weight);
    }
  } catch(...) {}

  if(event.topjets->size() > 1) {
    hist("mjj")->Fill(inv_mass_2(event.topjets->at(0).v4() + event.topjets->at(1).v4()));
  }


  if(debug) cout << "Finished Tstar Hists!" << endl;

}

TstarTstarHists::~TstarTstarHists(){}
