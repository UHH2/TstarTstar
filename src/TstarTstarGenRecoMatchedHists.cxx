#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/core/include/Event.h"
#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>
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


TstarTstarGenRecoMatchedHists::TstarTstarGenRecoMatchedHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

  // book all histograms here
  book<TH1F>("dR_ele", "dR_{e}(GEN,RECO)", 40, 0, 1.2);
  book<TH1F>("Pt_ratio_ele", "p_{T}^{e, RECO}/p_{T}^{e, GEN}", 40, 0, 2);
  book<TH1F>("Pt_ele", "p_{T}^{e, RECO matched}", 50, 0, 1100);
  book<TH1F>("Eta_ele", "#eta^{e, RECO matched}", 50, -5.2, 5.2);

  book<TH1F>("dR_mu", "dR_{#mu}(GEN,RECO)", 40, 0, 1.2);
  book<TH1F>("Pt_ratio_mu", "p_{T}^{#mu, RECO}/p_{T}^{#mu, GEN}", 40, 0, 2);
  book<TH1F>("Pt_mu", "p_{T}^{#mu, RECO matched}", 50, 0, 1100);
  book<TH1F>("Eta_mu", "#eta^{#mu, RECO matched}", 50, -5.2, 5.2);

  book<TH1F>("dR_photon", "dR_{#gamma}(GEN,RECO)", 40, 0, 1.2);
  book<TH1F>("Pt_ratio_photon", "p_{T}^{#gamma, RECO}/p_{T}^{#gamma, GEN}", 40, 0, 2);
  book<TH1F>("Pt_photon", "p_{T}^{#gamma, RECO matched}", 50, 0, 2500);
  book<TH1F>("Eta_photon", "#eta^{#gamma, RECO matched}", 50, -5.2, 5.2);

  book<TH1F>("dPhi_neutrino", "d#phi_{#nu}(GEN,RECO)", 40, 0, 4.0);
  book<TH1F>("Pt_ratio_neutrino", "p_{T}^{#nu, RECO}/p_{T}^{#nu, GEN}", 40, 0, 2);

  book<TH1F>("Pt_neutrino", "p_{T}^{#nu, RECO}", 100, 0, 2500);

  book<TH1F>("AK4_dR_b", "dR_{b}(GEN,AK4 RECO)", 20, 0, 1.2);
  book<TH1F>("AK4_Pt_ratio_b", "p_{T}^{AK4 RECO b}/p_{T}^{b, GEN}", 20, 0, 2);
  book<TH1F>("AK4_Pt_b", "p_{T}^{AK4 RECO}_{b, RECO matched}", 100, 0, 2500);
  book<TH1F>("AK4_Eta_b", "#eta^{AK4 RECO}_{b, RECO matched}", 50, -5.2, 5.2);

  book<TH1F>("AK4_dR_gluon", "dR_{gluon}(GEN,AK4 RECO)", 20, 0, 1.2);
  book<TH1F>("AK4_id_dRmin_gluon", "ID of the close to gluon AK4 jets", 30, 0, 29);
  book<TH1F>("AK4_Pt_ratio_gluon", "p_{T}^{AK4 RECO gluon}/p_{T}^{gluon, GEN}", 20, 0, 2);
  book<TH1F>("AK4_Pt_gluon", "p_{T}^{AK4 gluon, RECO matched}", 100, 0, 2500);
  book<TH1F>("AK4_Eta_gluon", "#eta^{AK4 gluon, RECO matched}", 50, -5.2, 5.2);
  book<TH1F>("AK4_dR_top", "dR_{top}(GEN,AK4 RECO)", 20, 0, 1.2);
  book<TH1F>("AK4_Pt_ratio_top", "p_{T}^{AK4 RECO top}/p_{T}^{top, GEN}", 20, 0, 2);
  book<TH1F>("AK4_Pt_top", "p_{T}^{AK4 RECO}_{top, RECO matched}", 100, 0, 3300);
  book<TH1F>("AK4_Eta_top", "#eta^{AK4 RECO}_{top, RECO matched}", 50, -5.2, 5.2);
  book<TH1F>("N_b_AK4matched", "N_{AK4 jets} matched to GEN b quarks", 10, 0, 10);
  book<TH1F>("N_gluon_AK4matched", "N_{AK4 jets} matched to GEN gluons", 10, 0, 10);
  book<TH1F>("N_top_AK4matched", "N_{AK4 jets} matched to GEN tops", 10, 0, 10);

  book<TH1F>("AK8_dR_b", "dR_{b}(GEN,AK8 RECO)", 20, 0, 1.2);
  book<TH1F>("AK8_Pt_ratio_b", "p_{T}^{AK8 RECO b}/p_{T}^{b, GEN}", 20, 0, 2);
  book<TH1F>("AK8_Pt_b", "p_{T}^{AK8 RECO}_{b, RECO matched}", 100, 0, 2500);
  book<TH1F>("AK8_Eta_b", "#eta^{AK8 RECO}_{b, RECO matched}", 50, -5.2, 5.2);

  book<TH1F>("AK8_dR_gluon", "dR_{gluon}(GEN,AK8 RECO)", 20, 0, 1.2);
  book<TH1F>("AK8_id_dRmin_gluon", "ID of the close to gluon AK8 jets", 30, 0, 29);
  book<TH1F>("AK8_Pt_ratio_gluon", "p_{T}^{AK8 RECO gluon}/p_{T}^{gluon, GEN}", 20, 0, 2);

  book<TH1F>("AK8_Pt_gluon", "p_{T}^{AK8 RECO}_{gluon, RECO matched}", 100, 0, 2500);
  book<TH1F>("AK8_Eta_gluon", "#eta^{AK8 RECO}_{gluon, RECO matched}", 50, -5.2, 5.2);
  book<TH1F>("AK8_dR_top", "dR_{top}(GEN,AK8 RECO)", 20, 0, 1.2);
  book<TH1F>("AK8_Pt_ratio_top", "p_{T}^{AK8 RECO top}/p_{T}^{top, GEN}", 20, 0, 2);
  book<TH1F>("AK8_Pt_top", "p_{T}^{AK8 RECO}_{top, RECO matched}", 100, 0, 3300);
  book<TH1F>("AK8_Eta_top", "#eta^{AK8 RECO}_{top, RECO matched}", 50, -5.2, 5.2);

  book<TH1F>("N_b_AK8matched", "N_{AK8 jets} matched to GEN b quarks", 10, 0, 10);
  book<TH1F>("N_gluon_AK8matched", "N_{AK8 jets} matched to GEN gluons", 10, 0, 10);
  book<TH1F>("N_top_AK8matched", "N_{AK8 jets} matched to GEN tops", 10, 0, 10);

  book<TH1F>("topmatches", "Number of GEN quarks in leading top jet", 4, 0, 4);

  is_mc = ctx.get("dataset_type") == "MC";
}


void TstarTstarGenRecoMatchedHists::fill(const Event & event){
  //  cout<<"TstarTstarGenRecoMatchedHists::fill start"<<endl;
  if(!is_mc) return;
  assert(event.genparticles);
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  // For pdgID mapping see https://twiki.cern.ch/twiki/bin/view/Main/PdgId
  // Find gen particles from expected topology Tstar+Tstar -> (t+gluon)+(t+gluon) or Tstar+Tstar -> (t+photon)+(t+photon)
  GenParticle tstar, antitstar, top, antitop, gluon1, gluon2, photon1, photon2;
  bool found_tstar = false, found_antitstar = false, found_top = false, found_antitop = false;
  bool found_gluon1 = false, found_gluon2 = false, found_photon1 = false, found_photon2 = false;
  int n_gluons = 0, n_photons = 0;
  //  vector<int> tstar_lep_jets, tstar_had_jets;
  for(const GenParticle & gp : *event.genparticles){
    if(gp.pdgId() == 6){
      top = gp;
      found_top = true;
    }
    else if(gp.pdgId() == -6){
      antitop = gp;
      found_antitop = true;
    }
    else if(gp.pdgId() == 9000005){
      tstar = gp;
      found_tstar = true;
    }
    else if(gp.pdgId() == -9000005){
      antitstar = gp;
      found_antitstar = true;
    }
    else if(gp.pdgId() == 21 && gp.status()==23){//only gluons from Tstar decay
      n_gluons++;
      if(!found_gluon1){
	gluon1 = gp;
	found_gluon1 = true;
      }
      else{
	gluon2 = gp;
	found_gluon2 = true;
      }
    }
    else if(gp.pdgId() == 22 && gp.status()==23){//only photons from Tstar decay
      n_photons++;
      if(!found_photon1){
	photon1 = gp;
	found_photon1 = true;
      }
      else{
	photon2 = gp;
	found_photon2 = true;
      }
    }


  }
  //  cout<<"TstarTstarGenRecoMatchedHists::fill step 2"<<endl;
  if(!found_tstar || !found_antitstar) return;
  if(!found_top || !found_antitop) return;
  if(!found_gluon1 && !found_gluon2 && !found_photon1 && !found_photon2) return;

  //Find gen particles from t->W+b
  auto b1 = top.daughter(event.genparticles,1);
  auto W1 = top.daughter(event.genparticles,2);
  if(fabs(W1->pdgId()) == 5 && fabs(b1->pdgId()) == 24){
    b1 = top.daughter(event.genparticles,2);
    W1 = top.daughter(event.genparticles,1);
  }

  auto b2 = antitop.daughter(event.genparticles,1);
  auto W2 = antitop.daughter(event.genparticles,2);
  if(fabs(W2->pdgId()) == 5 && fabs(b2->pdgId()) == 24){
    b2 = antitop.daughter(event.genparticles,2);
    W2 = antitop.daughter(event.genparticles,1);
  }
  if(fabs(W1->pdgId())!=fabs(W2->pdgId()) || fabs(W1->pdgId())!=24) return;//make sure Ws are Ws
  if(fabs(b1->pdgId())!=fabs(b2->pdgId()) || fabs(b1->pdgId())!=5) return;//make sure bs are bs
  //  cout<<"W1 id = "<<fabs(W1->pdgId())<<" W2 id = "<<fabs(W2->pdgId())<<endl;
  //Find gen particles from W->l+neutrino (leptonic) or W->q/g+q/g (hadronic decay)
  auto W1d1 = W1->daughter(event.genparticles,1);
  auto W1d2 = W1->daughter(event.genparticles,2);
  auto W2d1 = W2->daughter(event.genparticles,1);
  auto W2d2 = W2->daughter(event.genparticles,2);
  // cout<<"Look at Ws daughters "<<endl;
  // cout<<"W1 id1 = "<<fabs(W1d1->pdgId())<<" W2 id1 = "<<fabs(W2d1->pdgId())<<endl;
  // cout<<"W1 id2 = "<<fabs(W1d2->pdgId())<<" W2 id2 = "<<fabs(W2d2->pdgId())<<endl;

  bool W1_islep = true, W2_islep = true;//== isTop_lep = false and isAntiTop_lep=false;
  if(fabs(W1d1->pdgId()) < 7 && fabs(W1d2->pdgId()) < 7) W1_islep = false;
  if(fabs(W2d1->pdgId()) < 7 && fabs(W2d2->pdgId()) < 7) W2_islep = false;
  if(W1_islep && W2_islep) return;//non of W decays hadronicly

  auto lepton = W1d1;
  auto neutrino = W1d2;
  //  GenParticle lepton, neutrino;
  if(W1_islep){
    if((abs(W1d1->pdgId()) == 11) || (abs(W1d1->pdgId()) == 13)){
      lepton = W1d1; neutrino = W1d2;
    }
    if((abs(W1d2->pdgId()) == 11) || (abs(W1d2->pdgId()) == 13)){
      lepton = W1d2; neutrino = W1d1;
    }
  }
  if(W2_islep){
    if((abs(W2d1->pdgId()) == 11) || (abs(W2d1->pdgId()) == 13)){
      lepton = W2d1; neutrino = W2d2;
    }
    if((abs(W2d2->pdgId()) == 11) || (abs(W2d2->pdgId()) == 13)){
      lepton = W2d2; neutrino = W2d1;
    }
  }
  //  cout<<"Lepton and neutrino are assigned"<<endl;


  // Preparation for finding best hypothesis (if already reconstructed)


  //Now let's find match from RECO objects
  // lepton
  if(fabs(lepton->pdgId()) == 11){
    //    cout<<"Our lepton is electron"<<endl;
    for(const auto & ele : *event.electrons){
      //      cout<<"Look we loop over RECO electrons"<<endl;
      if(deltaR(*lepton,ele) <= 0.2){
	hist("dR_ele")->Fill(deltaR(*lepton,ele), weight);
	hist("Pt_ratio_ele")->Fill(ele.pt()/lepton->pt(), weight);
	hist("Pt_ele")->Fill(ele.pt(), weight);
	hist("Eta_ele")->Fill(ele.eta(), weight);
      }
    }
  }
  else if(fabs(lepton->pdgId()) == 13){
    //    cout<<"Our lepton is muon"<<endl;
    for(const auto & mu : *event.muons){
      if(deltaR(*lepton,mu) <= 0.2){
	hist("dR_mu")->Fill(deltaR(*lepton,mu), weight);
	hist("Pt_ratio_mu")->Fill(mu.pt()/lepton->pt(), weight);
	hist("Pt_mu")->Fill(mu.pt(), weight);
	hist("Eta_mu")->Fill(mu.eta(), weight);
      }
    }
  }
  //neutrino
  hist("Pt_neutrino")->Fill(neutrino->pt(),weight);
  hist("Pt_ratio_neutrino")->Fill(event.met->pt()/neutrino->pt(),weight);
  hist("dPhi_neutrino")->Fill(uhh2::deltaPhi(*event.met,*neutrino),weight);

  //photons matching
  if(found_photon1){
    for(const auto & gamma : *event.photons){
      if(deltaR(photon1,gamma) <= 0.2){
	hist("dR_photon")->Fill(deltaR(photon1,gamma), weight);
	hist("Pt_ratio_photon")->Fill(gamma.pt()/photon1.pt(), weight);
	hist("Pt_photon")->Fill(gamma.pt(), weight);
	hist("Eta_photon")->Fill(gamma.eta(), weight);
      }
    }
  }
  if(found_photon2){
    for(const auto & gamma : *event.photons){
      if(deltaR(photon2,gamma) <= 0.2){
	hist("dR_photon")->Fill(deltaR(photon2,gamma), weight);
	hist("Pt_ratio_photon")->Fill(gamma.pt()/photon2.pt(), weight);
	hist("Pt_photon")->Fill(gamma.pt(), weight);
	hist("Eta_photon")->Fill(gamma.eta(), weight);
      }
    }
  }

  //  cout<<"Now let's go to jets"<<endl;
  //b-jets, W-hadronic and gluons
  // Consider AK4 jets first
  double dR_max = 0.2;
  int bAK4match = 0, gluonAK4match = 0, topAK4match=0;
  int jetid = 0;
  for(const auto & jet : *event.jets){
    if(deltaR(*b1,jet) <= dR_max){
      hist("AK4_dR_b")->Fill(deltaR(*b1,jet), weight);
      hist("AK4_Pt_ratio_b")->Fill(jet.pt()/b1->pt(), weight);
      hist("AK4_Pt_b")->Fill(jet.pt(), weight);
      hist("AK4_Eta_b")->Fill(jet.eta(), weight);
      bAK4match++;
    }
    if(deltaR(*b2,jet) <= dR_max){
      hist("AK4_dR_b")->Fill(deltaR(*b2,jet), weight);
      hist("AK4_Pt_ratio_b")->Fill(jet.pt()/b2->pt(), weight);
      hist("AK4_Pt_b")->Fill(jet.pt(), weight);
      hist("AK4_Eta_b")->Fill(jet.eta(), weight);
      bAK4match++;
    }
    if(found_gluon1 && deltaR(gluon1,jet) <= dR_max){
      hist("AK4_dR_gluon")->Fill(deltaR(gluon1,jet), weight);
      hist("AK4_Pt_ratio_gluon")->Fill(jet.pt()/gluon1.pt(), weight);
      hist("AK4_Pt_gluon")->Fill(jet.pt(), weight);
      hist("AK4_Eta_gluon")->Fill(jet.eta(), weight);
      hist("AK4_id_dRmin_gluon")->Fill(jetid, weight);
      gluonAK4match++;
    }
    if(found_gluon2 && deltaR(gluon2,jet) <= dR_max){
      hist("AK4_dR_gluon")->Fill(deltaR(gluon2,jet), weight);
      hist("AK4_Pt_ratio_gluon")->Fill(jet.pt()/gluon2.pt(), weight);
      hist("AK4_Pt_gluon")->Fill(jet.pt(), weight);
      hist("AK4_Eta_gluon")->Fill(jet.eta(), weight);
      hist("AK4_id_dRmin_gluon")->Fill(jetid, weight);
      gluonAK4match++;
    }
    if(!W1_islep && deltaR(top,jet) <= dR_max){//t->Wb,W->hadronic
      hist("AK4_dR_top")->Fill(deltaR(top,jet), weight);
      hist("AK4_Pt_ratio_top")->Fill(jet.pt()/top.pt(), weight);
      hist("AK4_Pt_top")->Fill(jet.pt(), weight);
      hist("AK4_Eta_top")->Fill(jet.eta(), weight);
      topAK4match++;
    }
    if(!W2_islep && deltaR(antitop,jet) <= dR_max){//tbar->Wb,W->hadronic
      hist("AK4_dR_top")->Fill(deltaR(antitop,jet), weight);
      hist("AK4_Pt_ratio_top")->Fill(jet.pt()/antitop.pt(), weight);
      hist("AK4_Pt_top")->Fill(jet.pt(), weight);
      hist("AK4_Eta_top")->Fill(jet.eta(), weight);
      topAK4match++;
    }
    jetid++;
  }
  hist("N_b_AK4matched")->Fill(bAK4match, weight);
  hist("N_gluon_AK4matched")->Fill(gluonAK4match, weight);
  hist("N_top_AK4matched")->Fill(topAK4match, weight);

  // Now consider AK8 jets
  dR_max = 0.4;
  int bAK8match = 0, gluonAK8match = 0, topAK8match = 0;
  int topjetid = 0;
  for(const auto & jet : *event.topjets){
    if(deltaR(*b1,jet) <= dR_max){
      hist("AK8_dR_b")->Fill(deltaR(*b1,jet), weight);
      hist("AK8_Pt_ratio_b")->Fill(jet.pt()/b1->pt(), weight);
      hist("AK8_Pt_b")->Fill(jet.pt(), weight);
      hist("AK8_Eta_b")->Fill(jet.eta(), weight);
      bAK8match++;
    }
    if(deltaR(*b2,jet) <= dR_max){
      hist("AK8_dR_b")->Fill(deltaR(*b2,jet), weight);
      hist("AK8_Pt_ratio_b")->Fill(jet.pt()/b2->pt(), weight);
      hist("AK8_Pt_b")->Fill(jet.pt(), weight);
      hist("AK8_Eta_b")->Fill(jet.eta(), weight);
      bAK8match++;
    }
    if(found_gluon1 && deltaR(gluon1,jet) <= dR_max){
      hist("AK8_dR_gluon")->Fill(deltaR(gluon1,jet), weight);
      hist("AK8_Pt_ratio_gluon")->Fill(jet.pt()/gluon1.pt(), weight);
      hist("AK8_Pt_gluon")->Fill(jet.pt(), weight);
      hist("AK8_Eta_gluon")->Fill(jet.eta(), weight);
      hist("AK8_id_dRmin_gluon")->Fill(topjetid, weight);
      gluonAK8match++;
    }
    if(found_gluon2 && deltaR(gluon2,jet) <= dR_max){
      hist("AK8_dR_gluon")->Fill(deltaR(gluon2,jet), weight);
      hist("AK8_Pt_ratio_gluon")->Fill(jet.pt()/gluon2.pt(), weight);
      hist("AK8_Pt_gluon")->Fill(jet.pt(), weight);
      hist("AK8_Eta_gluon")->Fill(jet.eta(), weight);
      hist("AK8_id_dRmin_gluon")->Fill(topjetid, weight);
      gluonAK8match++;
    }
    if(!W1_islep && deltaR(top,jet) <= dR_max){
      hist("AK8_dR_top")->Fill(deltaR(top,jet), weight);
      hist("AK8_Pt_ratio_top")->Fill(jet.pt()/top.pt(), weight);
      hist("AK8_Pt_top")->Fill(jet.pt(), weight);
      hist("AK8_Eta_top")->Fill(jet.eta(), weight);
      topAK8match++;
    }
    if(!W2_islep && deltaR(antitop,jet) <= dR_max){
      hist("AK8_dR_top")->Fill(deltaR(antitop,jet), weight);
      hist("AK8_Pt_ratio_top")->Fill(jet.pt()/antitop.pt(), weight);
      hist("AK8_Pt_top")->Fill(jet.pt(), weight);
      hist("AK8_Eta_top")->Fill(jet.eta(), weight);
      topAK8match++;
    }
    topjetid++;
  }
  hist("N_b_AK8matched")->Fill(bAK8match, weight);
  hist("N_gluon_AK8matched")->Fill(gluonAK8match, weight);
  hist("N_top_AK8matched")->Fill(topAK8match, weight);

  // get ttbargen
  TTbarGen ttbargen = event.get(h_ttbargen);
  if(!ttbargen.IsSemiLeptonicDecay()) return;

  // find jet matching hadtop
  bool topjetfound = false;
  TopJet topjet;
  double dR_best = 999;
  for(const auto & jet : *event.topjets){
    double dR_cur = deltaR(jet, ttbargen.TopHad());
    if(dR_cur < dR_best) {
      topjetfound = true;
      dR_best = dR_cur;
      topjet = jet;
    }
  }
  if(topjetfound) {
    // calculate jet radius
    double R_topjet = (topjet.pt() > 0) ? 600/(topjet.pt()) : 1.5;
    if(R_topjet > 1.5) R_topjet = 1.5;
    if(R_topjet < 0.1) R_topjet = 0.1;
    // check overlaps
    int mergedCount = 0;
    if(deltaR(topjet, ttbargen.BHad()) < R_topjet) mergedCount++;
    if(deltaR(topjet, ttbargen.Q1()) < R_topjet) mergedCount++;
    if(deltaR(topjet, ttbargen.Q2()) < R_topjet) mergedCount++;
    hist("topmatches")->Fill(mergedCount, weight);
  }
  else hist("topmatches")->Fill(0., weight);

  //  cout<<"TstarTstarGenRecoMatchedHists: All hists are filled "<<endl;

 //  book<TH1F>("Lep_dR_gamma", "dR_{#gamma}(GEN,RECO)", 20, 0, 1.2);
 //  book<TH1F>("Lep_Pt_ratio_gamma", "p_{T}^{#gamma}/p_{T}^{#gamma}", 20, 0, 2);

}

TstarTstarGenRecoMatchedHists::~TstarTstarGenRecoMatchedHists(){}
