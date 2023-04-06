#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/TstarTstar/include/ReconstructionTstarHypothesis.h"
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/core/include/LorentzVector.h>
#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>

using namespace std;
using namespace uhh2;

namespace {

  // invariant mass of a lorentzVector, but safe for timelike / spacelike vectors
  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }


  double chi2Function(ReconstructionTstarHypothesis hyp, bool is_HOTVR){
    double chi2 = 0;

    ReconstructionHypothesis ttbar_hyp = hyp.ttbar_hyp();

    double mtophad_ = 170.5;
    double sigmatophad_ = 28.02;
    double chi2_top_had;
    if(is_HOTVR) chi2_top_had = pow((inv_mass(ttbar_hyp.tophad_v4()) - (mtophad_))/(sigmatophad_) ,2);
    else chi2_top_had = pow((ttbar_hyp.tophad_topjet_ptr()->softdropmass() - (mtophad_))/(sigmatophad_) ,2);
    chi2 += chi2_top_had;

    double mtoplep_ = 173.8;
    double sigmatoplep_ = 25.56;
    double chi2_top_lep = pow((inv_mass(ttbar_hyp.toplep_v4()) - (mtoplep_))/(sigmatoplep_) ,2);
    chi2 += chi2_top_lep;

    double sigma_deltaM_Tstar = 0.1749;
    double mean_deltaM_Tstar = 0.3513;
    double deltaM_norm = (inv_mass(hyp.tstarlep_v4()) - inv_mass(hyp.tstarhad_v4()))*2/(inv_mass(hyp.tstarlep_v4()) + inv_mass(hyp.tstarhad_v4()));
    double chi2_deltaM_Tstar = pow( (deltaM_norm - mean_deltaM_Tstar)/(sigma_deltaM_Tstar) ,2);
    chi2 += chi2_deltaM_Tstar;

    //cout << "Hyp detailed chi2:    top_had = " << chi2_top_had << "    top_lep = " << chi2_top_lep << "    deltaM_Tstar = " << chi2_deltaM_Tstar << endl;

    return chi2;
  }


}

bool TopTagMassWindow::operator()(const TopJet & topjet, const uhh2::Event &) const {
  //auto subjets = topjet.subjets();
  //if(subjets.size() < 3) return false;

  float mjet = topjet.softdropmass();

  if(mjet < m_mlower) return false;
  if(mjet > m_mupper) return false;

  return true;
}

TopTagMassWindow::TopTagMassWindow(double mlower, double mupper): m_mlower(mlower), m_mupper(mupper) {}

bool TopTagSubbtag::operator()(const TopJet & topjet, const uhh2::Event &) const {
  auto subjets = topjet.subjets();
  //if(subjets.size() < 3) return false;

  for(const auto & subjet : subjets){
    if(subjet.btag_DeepCSV() > m_btag){return true;}
  }

  return false;
}

TopTagSubbtag::TopTagSubbtag(double btag): m_btag(btag) {}


TstarTstar_tgtg_TopTag_Reconstruction::TstarTstar_tgtg_TopTag_Reconstruction(Context & ctx, const NeutrinoReconstructionMethod & neutrinofunction, TopJetId tjetid):
  m_neutrinofunction(neutrinofunction), topjetID_(tjetid) {

  // If a min_DR is -1, it is calculated live.
  // this is useful for HOTVR jets with variable radius

  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_tstartstar_hyp_vector = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>("TstarTstar_Hyp_Vector");
  h_flag_toptagevent = ctx.get_handle<int>("flag_toptagevent");
}

// #########################################
// ## Begin of Hypothesis Creation Module ##
// #########################################


bool TstarTstar_tgtg_TopTag_Reconstruction::process(uhh2::Event& event){

  bool debug = false;

  bool is_toptagevt = event.get(h_flag_toptagevent);

  if(debug){cout << "Hello World from Tstartstar_Full_Toptag_Reconstruction!" << endl;}

  // Build all posible hypothesieseses
  // 1) Hadronic Top candidates: ttagged AK8 jets
  // 2) Leptonic Top jet candidates: btagged AK4 jets
  // 3) Leptonic Top other stuff: MET->Neutrinofunction and Lepton
  // 4) Gluons: AK8 jets that are not the ttagged jet and not overlapping with AK4 from lep. top

  // Memory Structure
  // (to keep relatively compatible with present code:
  // Vector of ReconstructionTstarHypothesis
  // // ttbar information in a ReconstructionHypothesis.
  // // Gluon information
  // // Place for a chi2

  std::vector<ReconstructionTstarHypothesis> TstarTstar_hyp_vector;

  if(debug){cout << "Starting to find TopTagged AK8 Jets for the Hadronic Top." << endl;}

  assert(event.jets && event.topjets);
  assert(event.met);

  if(event.topjets->size() < 3) return false; // need three Fat jets for Reconstruction
  if(event.jets->size() < 1) return false; // need 1 ak4 jet for leptonic b

  const Particle& lepton = event.get(h_primlep); // Primary Lepton has to be set
  std::vector<LorentzVector> neutrinos = m_neutrinofunction(lepton.v4(), event.met->v4());

  for(uint i = 0; i < event.topjets->size(); i++){
    TopJet tj = event.topjets->at(i);

    if(is_toptagevt && !topjetID_(tj, event)) continue; // loop over all top tagged jets if one is present

    if(debug){cout << "Finding Gluon candidates..." << endl;}
    std::vector<LorentzVector> gluonCands;
    for(uint j = 0; j < event.topjets->size(); j++){
       TopJet other_tj = event.topjets->at(j);

       if(i == j) continue;
       else if(gluonCands.size() == 2) continue; // only take leading two gluon candidates. shortens calculation time and improves result.
       gluonCands.push_back(event.topjets->at(j).v4());
     }

     if(debug) cout << "Check gluons... " << endl;
     if(gluonCands.size() < 2) {
       if(debug) cout << "Not enough gluons" << endl;
       continue;
     }

    for(uint g1 = 0; g1 < gluonCands.size(); g1++){
   	 for(uint g2 = 0; g2 < gluonCands.size(); g2++){
   	   if(g1 == g2) continue;

    if(debug){cout << "Finding bjet candidate for leptonic top." << endl;}
    std::vector<const Jet*> tlep_jets;
    tlep_jets.reserve(event.jets->size());
    double R_toptaggedjet = (tj.pt() > 0) ? 600/(tj.pt()) : 1.5;
    if(R_toptaggedjet > 1.5) R_toptaggedjet = 1.5;
    if(R_toptaggedjet < 0.1) R_toptaggedjet = 0.1;
    LorentzVector gluon_cand_1 = gluonCands.at(g1);
    double R_g1 = (gluon_cand_1.pt() > 0) ? 600/(gluon_cand_1.pt()) : 1.5;
    if(R_g1 > 1.5) R_g1 = 1.5;
    if(R_g1 < 0.1) R_g1 = 0.1;
    LorentzVector gluon_cand_2 = gluonCands.at(g2);
    double R_g2 = (gluon_cand_2.pt() > 0) ? 600/(gluon_cand_2.pt()) : 1.5;
    if(R_g2 > 1.5) R_g2 = 1.5;
    if(R_g2 < 0.1) R_g2 = 0.1;
    for(const auto & jet : *event.jets){
      if((deltaR(tj, jet) > R_toptaggedjet+0.4) && (deltaR(gluon_cand_1, jet) > R_g1+0.4) && (deltaR(gluon_cand_2, jet) > R_g2+0.4)){tlep_jets.push_back(&jet);}
    }

    for(const auto & tlep_jet : tlep_jets){
      if(tlep_jet->btag_DeepCSV() <= 0.2219) continue;
      // TODO Save Btag in ReconstructionTstarHypothesis
      if(debug){cout << "Finding Neutrino candidates..." << endl;}
      for(const auto& neutrino_p4 : neutrinos){ // loop over Neutrinos

	if(debug){cout << "Starting Hypothesis construction..." << endl;}
	ReconstructionHypothesis ttbar_hyp;

	LorentzVector toplep_v4(lepton.v4() + neutrino_p4 + tlep_jet->v4());
	ttbar_hyp.set_lepton(lepton);
	ttbar_hyp.set_neutrino_v4(neutrino_p4);
	ttbar_hyp.add_toplep_jet(*tlep_jet);
	ttbar_hyp.set_blep_v4(tlep_jet->v4());

	LorentzVector tophad_v4(tj.v4());
	ttbar_hyp.add_tophad_jet(tj);
	ttbar_hyp.set_tophad_topjet_ptr(&event.topjets->at(i));

	ttbar_hyp.set_toplep_v4(toplep_v4);
	ttbar_hyp.set_tophad_v4(tophad_v4);

	    LorentzVector gluon1_v4 = gluonCands.at(g1);
	    LorentzVector gluon2_v4 = gluonCands.at(g2);

	    ReconstructionTstarHypothesis hyp;

	    hyp.set_ttbar(ttbar_hyp);

	    LorentzVector tstarlep_v4(toplep_v4 + gluon1_v4);
	    LorentzVector tstarhad_v4(tophad_v4 + gluon2_v4);

	    hyp.set_tstarlep_v4(tstarlep_v4);
	    hyp.set_tstarhad_v4(tstarhad_v4);

	    hyp.set_tstar1gluon_v4(tstarlep_v4);
	    hyp.set_tstar2gluon_v4(tstarhad_v4);

	    hyp.set_gluon1_v4(gluon1_v4);
	    hyp.set_gluon2_v4(gluon2_v4);

	    if(debug){cout << "Writing Hypothesis..." << endl;}
	    TstarTstar_hyp_vector.push_back(hyp);

	  }
	}
      }
    }
  }


  if(debug){cout << "Done, return to main" << endl;}
  event.set(h_tstartstar_hyp_vector, TstarTstar_hyp_vector);

  if(TstarTstar_hyp_vector.size() == 0){return false;}

  return true;
}

// ############################ new stuff begins here ########################

TstarTstar_tgtg_AK4_Reconstruction::TstarTstar_tgtg_AK4_Reconstruction(Context & ctx, const NeutrinoReconstructionMethod & neutrinofunction, TopJetId tjetid):
  m_neutrinofunction(neutrinofunction), topjetID_(tjetid) {

  // If a min_DR is -1, it is calculated live.
  // this is useful for HOTVR jets with variable radius

  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_tstartstar_hyp_vector = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>("TstarTstar_Hyp_Vector");
  h_flag_toptagevent = ctx.get_handle<int>("flag_toptagevent");
}

// #########################################
// ## Begin of Hypothesis Creation Module ##
// #########################################


bool TstarTstar_tgtg_AK4_Reconstruction::process(uhh2::Event& event){

  bool debug = false;

  bool is_toptagevt = event.get(h_flag_toptagevent);

  if(debug){cout << "Hello World from Tstartstar_AK4_Reconstruction!" << endl;}

  // Build all posible hypothesieseses
  // 1) Hadronic Top candidates: ttagged AK8 jets
  // 2) Leptonic Top jet candidates: btagged AK4 jets
  // 3) Leptonic Top other stuff: MET->Neutrinofunction and Lepton
  // 4) Gluons: AK8 jets that are not the ttagged jet and not overlapping with AK4 from lep. top

  // Memory Structure
  // (to keep relatively compatible with present code:
  // Vector of ReconstructionTstarHypothesis
  // // ttbar information in a ReconstructionHypothesis.
  // // Gluon information
  // // Place for a chi2

  std::vector<ReconstructionTstarHypothesis> TstarTstar_hyp_vector;

  if(debug){cout << "Starting to find TopTagged AK8 Jets for the Hadronic Top." << endl;}

  assert(event.jets && event.topjets);
  assert(event.met);

  const Particle& lepton = event.get(h_primlep); // Primary Lepton has to be set
  std::vector<LorentzVector> neutrinos = m_neutrinofunction(lepton.v4(), event.met->v4());

  for(uint i = 0; i < event.topjets->size(); i++){
    TopJet tj = event.topjets->at(i);

    if(is_toptagevt && !topjetID_(tj, event)) continue; // loop over all top tagged jets if one is present

    double R_toptaggedjet = (tj.pt() > 0) ? 600/(tj.pt()) : 1.5;
    if(R_toptaggedjet > 1.5) R_toptaggedjet = 1.5;
    if(R_toptaggedjet < 0.1) R_toptaggedjet = 0.1;

    if(debug){cout << "Finding bjet candidate for leptonic top." << endl;}
    for(uint k = 0; k < event.jets->size(); k++){
      Jet bj = event.jets->at(k);
      if(deltaR(tj, bj) <= R_toptaggedjet+0.4) continue;
      if(bj.btag_DeepCSV() <= 0.2219) continue;

      std::vector<LorentzVector> gluonCands;
      for(uint j = 0; j < event.jets->size(); j++){
        if (k == j) continue;
	Jet jet = event.jets->at(j);
	if(deltaR(tj, jet) > R_toptaggedjet+0.4) gluonCands.push_back(event.jets->at(j).v4());
      }
      if(gluonCands.size() < 2) continue;

      for(uint g1 = 0; g1 < gluonCands.size(); g1++){
	for(uint g2 = 0; g2 < gluonCands.size(); g2++){
	  if(g1 == g2) continue;

	  if(debug){cout << "Finding Neutrino candidates..." << endl;}
	  for(const auto& neutrino_p4 : neutrinos){ // loop over Neutrinos

	    if(debug){cout << "Starting Hypothesis construction..." << endl;}
	    ReconstructionHypothesis ttbar_hyp;
	    LorentzVector toplep_v4(lepton.v4() + neutrino_p4 + bj.v4());
	    ttbar_hyp.set_lepton(lepton);
	    ttbar_hyp.set_neutrino_v4(neutrino_p4);
	    ttbar_hyp.add_toplep_jet(bj);
	    ttbar_hyp.set_blep_v4(bj.v4());

	    LorentzVector tophad_v4(tj.v4());
	    ttbar_hyp.add_tophad_jet(tj);
	    ttbar_hyp.set_tophad_topjet_ptr(&event.topjets->at(i));

	    ttbar_hyp.set_toplep_v4(toplep_v4);
	    ttbar_hyp.set_tophad_v4(tophad_v4);

	    LorentzVector gluon1_v4 = gluonCands.at(g1);
	    LorentzVector gluon2_v4 = gluonCands.at(g2);

	    ReconstructionTstarHypothesis hyp;

	    hyp.set_ttbar(ttbar_hyp);

	    LorentzVector tstarlep_v4(toplep_v4 + gluon1_v4);
	    LorentzVector tstarhad_v4(tophad_v4 + gluon2_v4);

	    hyp.set_tstarlep_v4(tstarlep_v4);
	    hyp.set_tstarhad_v4(tstarhad_v4);

	    hyp.set_tstar1gluon_v4(tstarlep_v4);
	    hyp.set_tstar2gluon_v4(tstarhad_v4);

	    hyp.set_gluon1_v4(gluon1_v4);
	    hyp.set_gluon2_v4(gluon2_v4);

	    if(debug){cout << "Writing Hypothesis..." << endl;}
	    TstarTstar_hyp_vector.push_back(hyp);
	  }
	}
      }
    }
  }

  if(debug){cout << "Done, return to main" << endl;}
  event.set(h_tstartstar_hyp_vector, TstarTstar_hyp_vector);

  if(TstarTstar_hyp_vector.size() == 0){return false;}

  return true;
}

// ######################## new stuff ends here ############

TstarTstar_Discrimination::TstarTstar_Discrimination(Context & ctx){
  h_tstartstar_hyp_vector = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>("TstarTstar_Hyp_Vector");
  h_tstartstar_hyp = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_Hyp");

}


bool TstarTstar_Discrimination::process(uhh2::Event& event){

  bool debug = false;

  // Read in all Hypothesiseses
  // save chi2 in each Hyp.
  // Chose hyp with lowest chi2
  // save in h_tstartstar_hyp

  std::vector<ReconstructionTstarHypothesis> hyp_vector = event.get(h_tstartstar_hyp_vector);

  double chi2_min=1e9;
  ReconstructionTstarHypothesis hyp_best;
  for(auto & hyp : hyp_vector){
    double chi2 = chi2Function(hyp, true);
    hyp.set_chi2(chi2);

    if(chi2 < chi2_min){
      chi2_min = chi2;
      hyp_best = hyp;
    }
  }

  if(debug)cout << "lowewst chi2 = " << chi2_min << endl;

  event.set(h_tstartstar_hyp_vector, hyp_vector);
  event.set(h_tstartstar_hyp, hyp_best);

  return true;
}
