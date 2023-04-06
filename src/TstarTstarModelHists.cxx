#include "UHH2/TstarTstar/include/TstarTstarModelHists.h"
#include "UHH2/core/include/Event.h"
#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

namespace {
  // invariant mass of a lorentzVector, but safe for timelike / spacelike vectors
  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }
}


TstarTstarModelHists::TstarTstarModelHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here

  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

  book<TH1F>("weight", "weight", 40, 0, 2);

  // #################
  // ## Basic hists ##
  // #################

  book<TH1F>("N_Tstar", "N_{T*} gen", 5, 0, 5);
  book<TH1F>("Tstar_spin", "Spin^{T*} gen", 30, -3, 3);
  book<TH1F>("Tstar_mass", "M^{T*} gen", 250, 0, 5000);
  book<TH1F>("Tstar_pt", "p_{T}^{T*} gen", 50, 0, 2000);
  book<TH1F>("Tstar_eta", "#eta^{T*} gen", 50, -5.2, 5.2);
  book<TH1F>("Tstar_phi", "#phi^{T*} gen", 30, -3.14, 3.14);

  book<TH1F>("top_spin", "Spin^{t} gen", 30, -3, 3);
  book<TH1F>("top_mass", "M^{t} gen", 100, 0, 2000);
  book<TH1F>("top_pt", "p_{T}^{t} gen", 50, 0, 2000);
  book<TH1F>("top_pt_initial", "p_{T}^{t (initial)} gen", 50, 0, 2000);
  book<TH1F>("top_pt_fromTstar", "p_{T}^{t (from T*)} gen", 50, 0, 2000);
  book<TH1F>("top_eta", "#eta^{t} gen", 52, -5.2, 5.2);
  book<TH1F>("top_phi", "#phi^{t} gen", 30, -3.14, 3.14);

  book<TH1F>("N_gluon", "N_{g} gen", 5, 0, 5);
  book<TH1F>("gluon_pt", "p_{T}^{g} gen", 50, 0, 2000);
  book<TH1F>("gluon_eta", "#eta^{g} gen", 52, -5.2, 5.2);
  book<TH1F>("gluon_phi", "#phi^{g} gen", 30, -3.14, 3.14);

  book<TH1F>("electron_pt", "p_{T}^{e} gen", 50, 0, 2000);
  book<TH1F>("electron_eta", "#eta^{e} gen", 52, -5.2, 5.2);
  book<TH1F>("electron_phi", "#phi^{e} gen", 30, -3.14, 3.14);

  book<TH1F>("muon_pt", "p_{T}^{#mu} gen", 50, 0, 2000);
  book<TH1F>("muon_eta", "#eta^{#mu} gen", 52, -5.2, 5.2);
  book<TH1F>("muon_phi", "#phi^{#mu} gen", 30, -3.14, 3.14);

  book<TH1F>("b_pt", "p_{T}^{b} gen", 50, 0, 2000);
  book<TH1F>("b_eta", "#eta^{b} gen", 52, -5.2, 5.2);
  book<TH1F>("b_phi", "#phi^{b} gen", 30, -3.14, 3.14);

  book<TH1F>("q_pt", "p_{T}^{q} gen", 50, 0, 2000);
  book<TH1F>("q_eta", "#eta^{q} gen", 52, -5.2, 5.2);
  book<TH1F>("q_phi", "#phi^{q} gen", 30, -3.14, 3.14);

  // #################
  // ## Other hists ##
  // #################

  book<TH1F>("dR_tops", "#DeltaR_{T* T* or t T*} gen", 30, 0, 6);
  book<TH1F>("deta_tops", "#Delta#eta_{T* T* or t T*} gen", 30, 0, 6);
  book<TH1F>("dphi_tops", "#Delta#phi_{T* T* or t T*} gen", 20, 0, 4);

  book<TH1F>("invmass_top", "m_{t t} gen", 100, 0, 2000);
  book<TH1F>("invmass_gluon", "m_{g g} gen", 100, 0, 2000);

  book<TH1F>("dR_gluons", "#DeltaR_{g g} gen", 30, 0, 6);
  book<TH1F>("deta_gluons", "#Delta#eta_{g g} gen", 30, 0, 6);
  book<TH1F>("dphi_gluons", "#Delta#phi_{g g} gen", 20, 0, 4);

  book<TH1F>("dR_gluon_top", "#DeltaR_{g t} gen", 30, 0, 6);
  book<TH1F>("deta_gluon_top", "#Delta#eta_{g t} gen", 30, 0, 6);
  book<TH1F>("dphi_gluon_top", "#Delta#phi_{g t} gen", 20, 0, 4);

  book<TH1F>("dR_qs", "#DeltaR_{q q} gen", 30, 0, 6);
  book<TH1F>("deta_qs", "#Delta#eta_{q q} gen", 30, 0, 6);
  book<TH1F>("dphi_qs", "#Delta#phi_{q q} gen", 20, 0, 4);

  book<TH1F>("recoTstar_mass", "m_{t g} gen", 100, 0, 2000);
  book<TH1F>("recoTstar_pt", "p_{T}^{t g} gen", 50, 0, 2000);

  is_mc = ctx.get("dataset_type") == "MC";

}


void TstarTstarModelHists::fill(const Event & event){
  if(!is_mc) return;

  assert(event.genparticles);
  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  hist("weight")->Fill(weight);

  // finding Tstars and determining which model is used
  std::vector<GenParticle> Tstars;
  bool isSpin12 = false;
  for(const GenParticle & gp : *event.genparticles){
    if((gp.pdgId() == 25001) || (gp.pdgId() == -25001)) {
      isSpin12 = true;
      Tstars.push_back(gp);
    }
    if((gp.pdgId() == 600) || (gp.pdgId() == -600)) {
      isSpin12 = true;
      Tstars.push_back(gp);
    }
    else if((gp.pdgId() == 9000005 && (gp.status()==23 || gp.status()==22)) || (gp.pdgId() == -9000005 && (gp.status()==23 || gp.status()==22))) {
      if(isSpin12) std::cout << "Error, models mixed?" << std::endl;
      Tstars.push_back(gp);
    }
  }

  // N_Tstar to determine whether pair production or single production
  int N_Tstar = Tstars.size();

  // filling gluons, depending on used model
  std::vector<GenParticle> gluons;
  if(isSpin12) {
    for(const GenParticle & gp : *event.genparticles){
      if((gp.pdgId() == 21) && (gp.status() > 0)) gluons.push_back(gp);
    }
  }
  else {
    for(const GenParticle & gp : *event.genparticles){
      if((gp.pdgId() == 21) && (gp.status() == 23)) gluons.push_back(gp);
    }
  }

  if(N_Tstar == 2){
    hist("dR_tops")->Fill(deltaR(Tstars.at(0), Tstars.at(1)), weight);
    hist("deta_tops")->Fill(abs(Tstars.at(0).eta() - Tstars.at(1).eta()), weight);
    hist("dphi_tops")->Fill(deltaPhi(Tstars.at(0), Tstars.at(1)), weight);
  }

  hist("N_Tstar")->Fill(N_Tstar, weight);
  for (const auto & Tstar : Tstars){
    hist("Tstar_spin")->Fill(Tstar.spin(), weight/N_Tstar);
    hist("Tstar_mass")->Fill(inv_mass(Tstar.v4()), weight/N_Tstar);
    hist("Tstar_pt")->Fill(Tstar.pt(), weight/N_Tstar);
    hist("Tstar_eta")->Fill(Tstar.eta(), weight/N_Tstar);
    hist("Tstar_phi")->Fill(Tstar.phi(), weight/N_Tstar);
  }

  // get ttbargen object
  TTbarGen ttbargen = event.get(h_ttbargen);
  if(!ttbargen.IsSemiLeptonicDecay()){return;}

  // top quarks
  hist("top_spin")->Fill(ttbargen.Top().spin(), weight/2);
  hist("top_mass")->Fill(inv_mass(ttbargen.Top().v4()), weight/2);
  hist("top_pt")->Fill(ttbargen.Top().pt(), weight/2);
  hist("top_eta")->Fill(ttbargen.Top().eta(), weight/2);
  hist("top_phi")->Fill(ttbargen.Top().phi(), weight/2);

  hist("top_spin")->Fill(ttbargen.Antitop().spin(), weight/2);
  hist("top_mass")->Fill(inv_mass(ttbargen.Antitop().v4()), weight/2);
  hist("top_pt")->Fill(ttbargen.Antitop().pt(), weight/2);
  hist("top_eta")->Fill(ttbargen.Antitop().eta(), weight/2);
  hist("top_phi")->Fill(ttbargen.Antitop().phi(), weight/2);

  std::vector<GenParticle> tops;
  tops.push_back(ttbargen.Top());
  tops.push_back(ttbargen.Antitop());
  int N_top = tops.size();
  int N_gluon = gluons.size();
  hist("N_gluon")->Fill(N_gluon, weight);
  int count = 0; // crosscheck, should in the end be the same as N_Tstar
  for (const auto & gluon : gluons){
    hist("gluon_pt")->Fill(gluon.pt(), weight/N_gluon);
    hist("gluon_eta")->Fill(gluon.eta(), weight/N_gluon);
    hist("gluon_phi")->Fill(gluon.phi(), weight/N_gluon);
    for (const auto & top : tops){
      if(top.mother1() == gluon.mother1()) {
        hist("dR_gluon_top")->Fill(deltaR(gluon, top), weight/N_Tstar);
        hist("deta_gluon_top")->Fill(abs(gluon.eta() - top.eta()), weight/N_Tstar);
        hist("dphi_gluon_top")->Fill(deltaPhi(gluon, top), weight/N_Tstar);
        hist("recoTstar_pt")->Fill((gluon.v4()+top.v4()).pt(), weight/N_Tstar);
        hist("recoTstar_mass")->Fill(inv_mass(gluon.v4()+top.v4()), weight/N_Tstar);
        hist("top_pt_fromTstar")->Fill(top.pt(), weight/N_top);
        count++;
      }
      else if(N_Tstar == 1) hist("top_pt_initial")->Fill(top.pt(), weight/N_top); // TODO this is an ugly hack. Fix it Finn!
    }
  }
  if(count > N_Tstar) std::cout << "Found too many top-gluon pairs." << std::endl;

  if(N_gluon == 2){
    hist("dR_gluons")->Fill(deltaR(gluons.at(0), gluons.at(1)), weight);
    hist("deta_gluons")->Fill(abs(gluons.at(0).eta() - gluons.at(1).eta()), weight);
    hist("dphi_gluons")->Fill(deltaPhi(gluons.at(0), gluons.at(1)), weight);
  }

  // invariant masses
  // currently only if exactly two tops / gluons are present (which should be the case in the hard process)
  if(N_top == 2) hist("invmass_top")->Fill(inv_mass(tops.at(0).v4() + tops.at(1).v4()), weight);
  if(N_gluon == 2) hist("invmass_gluon")->Fill(inv_mass(gluons.at(0).v4() + gluons.at(1).v4()), weight);

  if(abs(ttbargen.ChargedLepton().pdgId()) == 11) {
    hist("electron_pt")->Fill(ttbargen.ChargedLepton().pt(), weight);
    hist("electron_eta")->Fill(ttbargen.ChargedLepton().eta(), weight);
    hist("electron_phi")->Fill(ttbargen.ChargedLepton().phi(), weight);
  }
  else {
    hist("muon_pt")->Fill(ttbargen.ChargedLepton().pt(), weight);
    hist("muon_eta")->Fill(ttbargen.ChargedLepton().eta(), weight);
    hist("muon_phi")->Fill(ttbargen.ChargedLepton().phi(), weight);
  }

  std::vector<GenParticle> bs;
  bs.push_back(ttbargen.bTop());
  bs.push_back(ttbargen.bAntitop());
  int N_b = bs.size();
  for (const auto & b : bs){
    hist("b_pt")->Fill(b.pt(), weight/N_b);
    hist("b_eta")->Fill(b.eta(), weight/N_b);
    hist("b_phi")->Fill(b.phi(), weight/N_b);
  }

  std::vector<GenParticle> qs;
  qs.push_back(ttbargen.Q1());
  qs.push_back(ttbargen.Q2());
  int N_q = qs.size();
  for (const auto & q : qs){
    hist("q_pt")->Fill(q.pt(), weight/N_q);
    hist("q_eta")->Fill(q.eta(), weight/N_q);
    hist("q_phi")->Fill(q.phi(), weight/N_q);
  }

  if(N_q == 2){
    hist("dR_qs")->Fill(deltaR(qs.at(0), qs.at(1)), weight);
    hist("deta_qs")->Fill(abs(qs.at(0).eta() - qs.at(1).eta()), weight);
    hist("dphi_qs")->Fill(deltaPhi(qs.at(0), qs.at(1)), weight);
  }

}

TstarTstarModelHists::~TstarTstarModelHists(){}
