#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

namespace {

  // invariant mass of a lorentzVector, but safe for timelike / spacelike vectors
  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

}

TstarTstarRecoTstarHists::TstarTstarRecoTstarHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  // ###### Handles ######
  h_tstartstar_hyp = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_Hyp");
  h_tstartstar_hyp_vector = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>("TstarTstar_Hyp_Vector");

  // ###### Histograms ######

  book<TH1F>("N_hyps","N_{hypotheses}", 25, 0, 50);
  book<TH1F>("N_hyps_2","N_{hypotheses}", 25, 0, 250);

  // #### M_Tstar
  book<TH1F>("M_Tstar", "m_{T*} [GeV]", 30, 0, 3000);
  book<TH1F>("newST", "Sum of all Fat jets plus lepton", 40, 0, 4000);
  // tgtg
  book<TH1F>("M_Tstar_lep", "m_{T*_{lep}} [GeV]", 50, 0, 3000);
  book<TH1F>("M_Tstar_had", "m_{T*_{had}} [GeV]", 50, 0, 3000);
  book<TH1F>("deltaM_Tstar", "Mass difference of T^{*}_{tgtg,lep} and T^{*}_{tgtg,had} [GeV]", 50, -2, 2);
  book<TH1F>("Chi2_bestchi2", "#chi^{2}_{best hyp}", 20, 0, 500);

  // #### Delta_R
  book<TH1F>("DeltaR_toplep_gluon", "#DeltaR (t_{lep}, g)", 15, 0, 6);
  book<TH1F>("DeltaR_tophad_gluon", "#DeltaR (t_{had}, g)", 15, 0, 6);

  // #### ttbar
  book<TH1F>("Mtoplep", "M_{toplep} [GeV]", 50, 0, 500);
  book<TH1F>("Mtophad", "M_{tophad} [GeV]", 50, 0, 500);
  book<TH1F>("t_lep_used_jets", "N_{Jets} (leptonic top)", 10, 0, 10);
  book<TH1F>("t_had_used_jets", "N_{Jets} (hadronic top)", 11, -1, 10);

  // #### Pt
  // tgtg
  book<TH1F>("Pt_gluon_Tstar_lep", "p_{T} gluon, T^{*}_{tgtg,lep} [GeV]", 50, 0, 1000);
  book<TH1F>("Pt_gluon_Tstar_had", "p_{T} gluon, T^{*}_{tgtg,had} [GeV]", 50, 0, 1000);
  book<TH1F>("Pt_ratio_Tstar_lep_to_Tstar_had", "(p_{T} T^{*}_{tgtg,lep})/(p_{T} T^{*}_{tgtg,had})", 50, 0, 1e1);

  // #### DNN Hists ####

  book<TH1F>("DNN_hadtopjet_pt", "p_{T, jet_{had top}} [GeV]", 50, 0, 1000);
  book<TH1F>("DNN_hadtopjet_eta", "#eta_{jet_{had top}} [GeV]", 30, -3, 3);
  book<TH1F>("DNN_hadtopjet_phi", "#phi_{jet_{had top}} [GeV]", 25, -3.2, 3.2);

  book<TH1F>("DNN_leptopjet_pt", "p_{T, jet_{lep top}} [GeV]", 50, 0, 1000);
  book<TH1F>("DNN_leptopjet_eta", "#eta_{jet_{lep top}} [GeV]", 30, -3, 3);
  book<TH1F>("DNN_leptopjet_phi", "#phi_{jet_{lep top}} [GeV]", 25, -3.2, 3.2);

  book<TH1F>("DNN_lepton_pt", "p_{T, lepton} [GeV]", 50, 0, 1000);
  book<TH1F>("DNN_lepton_eta", "#eta_{lepton} [GeV]", 30, -3, 3);
  book<TH1F>("DNN_lepton_phi", "#phi_{lepton} [GeV]", 25, -3.2, 3.2);

  book<TH1F>("DNN_MET_pt", "p_{T, MET} [GeV]", 50, 0, 1000);
  book<TH1F>("DNN_MET_eta", "#eta_{MET} [GeV]", 30, -3, 3);
  book<TH1F>("DNN_MET_phi", "#phi_{MET} [GeV]", 25, -3.2, 3.2);

  book<TH1F>("DNN_gluon1_pt", "p_{T, g_1} [GeV]", 50, 0, 1000);
  book<TH1F>("DNN_gluon1_eta", "#eta_{g_1} [GeV]", 30, -3, 3);
  book<TH1F>("DNN_gluon1_phi", "#phi_{g_1} [GeV]", 25, -3.2, 3.2);

  book<TH1F>("DNN_gluon2_pt", "p_{T, g_2} [GeV]", 50, 0, 1000);
  book<TH1F>("DNN_gluon2_eta", "#eta_{g_2} [GeV]", 30, -3, 3);
  book<TH1F>("DNN_gluon2_phi", "#phi_{g_2} [GeV]", 25, -3.2, 3.2);


}





void TstarTstarRecoTstarHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  bool debug = false;

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  // ###### Handles ######
  std::vector<ReconstructionTstarHypothesis> hyps = event.get(h_tstartstar_hyp_vector);
  ReconstructionHypothesis hyp = event.get(h_tstartstar_hyp).ttbar_hyp();
  ReconstructionTstarHypothesis tstar_hyp_best = event.get(h_tstartstar_hyp);

  if(debug){cout << "Starting filling histograms..." << endl;}

  // ###### Histgram Filling ######

  hist("N_hyps")->Fill(hyps.size(), weight);
  hist("N_hyps_2")->Fill(hyps.size(), weight);

  // #### M_Tstar
  // tgtg
  double mTstarlep = inv_mass(tstar_hyp_best.tstarlep_v4());
  double mTstarhad = inv_mass(tstar_hyp_best.tstarhad_v4());
  double mTstar = (mTstarlep+mTstarhad)/2;
  hist("M_Tstar")->Fill(mTstar, weight);

  double newST = 0;
  for(const auto & jet : *event.topjets){
    newST += jet.pt();
  }
  newST += hyp.lepton().v4().pt();
  hist("newST")->Fill(newST, weight);

  hist("M_Tstar_lep")->Fill(mTstarlep, weight);
  hist("M_Tstar_had")->Fill(mTstarhad, weight);
  hist("deltaM_Tstar")->Fill((mTstarlep-mTstarhad)*2/(mTstarlep+mTstarhad), weight);
  hist("Chi2_bestchi2")->Fill(tstar_hyp_best.chi2(), weight);

  if(debug){cout << "Finished writing Tstar masses." << endl;}

  // tgtg
  double DeltaR_toplep_gluon =  deltaR(tstar_hyp_best.gluon1_v4(),tstar_hyp_best.ttbar_hyp().toplep_v4());
  double DeltaR_tophad_gluon =  deltaR(tstar_hyp_best.gluon2_v4(),tstar_hyp_best.ttbar_hyp().tophad_v4());
  hist("DeltaR_toplep_gluon")->Fill(DeltaR_toplep_gluon, weight);
  hist("DeltaR_tophad_gluon")->Fill(DeltaR_tophad_gluon, weight);

  if(debug){cout << "Finished writing DeltaRs." << endl;}

  // #### ttbar
  hist("Mtoplep")->Fill(inv_mass(hyp.toplep_v4()), weight);
  hist("Mtophad")->Fill(inv_mass(hyp.tophad_v4()), weight);
  hist("t_lep_used_jets")->Fill(hyp.toplep_jets().size(), weight);
  if(!hyp.tophad_topjet_ptr()){hist("t_had_used_jets")->Fill(hyp.tophad_jets().size(), weight);}
  else {hist("t_had_used_jets")->Fill(-0.5, weight);}

  if(debug){cout << "Finished writing ttbar plots." << endl;}

  // #### Pt
  hist("Pt_gluon_Tstar_lep")->Fill(tstar_hyp_best.gluon1_v4().pt(), weight);
  hist("Pt_gluon_Tstar_had")->Fill(tstar_hyp_best.gluon2_v4().pt(), weight);
  hist("Pt_ratio_Tstar_lep_to_Tstar_had")->Fill(tstar_hyp_best.tstarlep_v4().pt()/tstar_hyp_best.tstarhad_v4().pt(), weight);

  // ##########################
  // Filling output for DNN
  // TODO put this in extra file that you pass an hypothesis
  if(debug) cout << "Start filling of DNN stuff" << endl;
  ReconstructionTstarHypothesis best_hyp = event.get(h_tstartstar_hyp);
  ReconstructionHypothesis best_hyp_ttbar = best_hyp.ttbar_hyp();

  hist("DNN_hadtopjet_pt")->Fill(best_hyp_ttbar.tophad_v4().pt(), weight);
  hist("DNN_hadtopjet_eta")->Fill(best_hyp_ttbar.tophad_v4().eta(), weight);
  hist("DNN_hadtopjet_phi")->Fill(best_hyp_ttbar.tophad_v4().phi(), weight);
  if(debug) cout << "Done with ttagjet." << endl;

  hist("DNN_lepton_pt")->Fill( best_hyp_ttbar.lepton().pt(), weight);
  hist("DNN_lepton_eta")->Fill( best_hyp_ttbar.lepton().eta(), weight);
  hist("DNN_lepton_phi")->Fill( best_hyp_ttbar.lepton().phi(), weight);
  if(debug) cout << "Done with lepton." << endl;

  hist("DNN_MET_pt")->Fill( best_hyp_ttbar.neutrino_v4().pt(), weight);
  hist("DNN_MET_eta")->Fill(best_hyp_ttbar.neutrino_v4().eta(), weight);
  hist("DNN_MET_phi")->Fill( best_hyp_ttbar.neutrino_v4().phi(), weight);
  if(debug) cout << "Done with neutrino." << endl;

  hist("DNN_leptopjet_pt")->Fill( best_hyp_ttbar.blep_v4().pt(), weight);
  hist("DNN_leptopjet_eta")->Fill( best_hyp_ttbar.blep_v4().eta(), weight);
  hist("DNN_leptopjet_phi")->Fill( best_hyp_ttbar.blep_v4().phi(), weight);
  if(debug) cout << "Done with leptopjet." << endl;

  hist("DNN_gluon1_pt")->Fill(best_hyp.gluon1_v4().pt(), weight);
  hist("DNN_gluon1_eta")->Fill( best_hyp.gluon1_v4().eta(), weight);
  hist("DNN_gluon1_phi")->Fill( best_hyp.gluon1_v4().phi(), weight);
  if(debug) cout << "Done with gluon1." << endl;

  hist("DNN_gluon2_pt")->Fill(best_hyp.gluon2_v4().pt(), weight);
  hist("DNN_gluon2_eta")->Fill( best_hyp.gluon2_v4().eta(), weight);
  hist("DNN_gluon2_phi")->Fill( best_hyp.gluon2_v4().phi(), weight);
  if(debug) cout << "Done with gluon2." << endl;

  if(debug){cout << "Finished filling hists, return to main!" << endl;}


}


void TstarTstarRecoTstarHists::fill_ttbarhyps(const Event & event, const ReconstructionHypothesis & hyp){

  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  hist("Mtoplep_all")->Fill(inv_mass(hyp.toplep_v4()), weight);
  hist("Mtophad_all")->Fill(inv_mass(hyp.tophad_v4()), weight);

}

TstarTstarRecoTstarHists::~TstarTstarRecoTstarHists(){}
