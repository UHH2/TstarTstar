#include "UHH2/TstarTstar/include/TstarTstarAllGenHists.h"
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
#include "Math/Boost.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzRotation.h"

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

namespace {
  // invariant mass of a lorentzVector, but safe for timelike / spacelike vectors
  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }
}


TstarTstarAllGenHists::TstarTstarAllGenHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here

  book<TH1F>("matches_found", "matches_found", 10, 0, 10);
  book<TH1F>("lab_pt_mother", "pt_mother", 20, 0, 2000);
  book<TH1F>("rest_pt_mother", "pt_mother", 20, 0, 2000);
  book<TH1F>("CS_pt_mother", "pt_mother", 20, 0, 2000);

  // lab frame
  book<TH1F>("lab_dR_top_matched_gluon", "lab_dR_top_matched_gluon", 20, 0, 10);
  book<TH1F>("lab_dphi_top_matched_gluon", "lab_dphi_top_matched_gluon", 24, 0, 8);
  book<TH1F>("lab_deta_top_matched_gluon", "lab_deta_top_matched_gluon", 24, 0, 8);

  book<TH1F>("lab_dR_top_matched_gluon_lowpt", "lab_dR_top_matched_gluon_lowpt", 20, 0, 10);
  book<TH1F>("lab_dphi_top_matched_gluon_lowpt", "lab_dphi_top_matched_gluon_lowpt", 24, 0, 8);
  book<TH1F>("lab_deta_top_matched_gluon_lowpt", "lab_deta_top_matched_gluon_lowpt", 24, 0, 8);

  book<TH1F>("lab_dR_top_matched_gluon_medpt", "lab_dR_top_matched_gluon_medpt", 20, 0, 10);
  book<TH1F>("lab_dphi_top_matched_gluon_medpt", "lab_dphi_top_matched_gluon_medpt", 24, 0, 8);
  book<TH1F>("lab_deta_top_matched_gluon_medpt", "lab_deta_top_matched_gluon_medpt", 24, -0, 8);

  book<TH1F>("lab_dR_top_matched_gluon_highpt", "lab_dR_top_matched_gluon_highpt", 20, 0, 10);
  book<TH1F>("lab_dphi_top_matched_gluon_highpt", "lab_dphi_top_matched_gluon_highpt", 24, 0, 8);
  book<TH1F>("lab_deta_top_matched_gluon_highpt", "lab_deta_top_matched_gluon_highpt", 24, 0, 8);

  book<TH1F>("lab_costheta_mother_gluon", "lab_costheta_mother_gluon", 24, -1, 1);
  book<TH1F>("lab_costheta_lowpt_mother_gluon", "lab_costheta_lowpt_mother_gluon", 24, -1, 1);
  book<TH1F>("lab_costheta_medpt_mother_gluon", "lab_costheta_medpt_mother_gluon", 24, -1, 1);
  book<TH1F>("lab_costheta_highpt_mother_gluon", "lab_costheta_highpt_mother_gluon", 24, -1, 1);

  book<TH1F>("lab_costheta_top_gluon", "lab_costheta_top_gluon", 24, -1, 1);
  book<TH1F>("lab_costheta_lowpt_top_gluon", "lab_costheta_lowpt_top_gluon", 24, -1, 1);
  book<TH1F>("lab_costheta_medpt_top_gluon", "lab_costheta_medpt_top_gluon", 24, -1, 1);
  book<TH1F>("lab_costheta_highpt_top_gluon", "lab_costheta_highpt_top_gluon", 24, -1, 1);


  // rest frame
  book<TH1F>("rest_dR_top_matched_gluon", "rest_dR_top_matched_gluon", 20, 0, 10);
  book<TH1F>("rest_dphi_top_matched_gluon", "rest_dphi_top_matched_gluon", 24, 0, 8);
  book<TH1F>("rest_deta_top_matched_gluon", "rest_deta_top_matched_gluon", 24, 0, 8);

  book<TH1F>("rest_dR_top_matched_gluon_lowpt", "rest_dR_top_matched_gluon_lowpt", 20, 0, 10);
  book<TH1F>("rest_dphi_top_matched_gluon_lowpt", "rest_dphi_top_matched_gluon_lowpt", 24, 0, 8);
  book<TH1F>("rest_deta_top_matched_gluon_lowpt", "rest_deta_top_matched_gluon_lowpt", 24, 0, 8);

  book<TH1F>("rest_dR_top_matched_gluon_medpt", "rest_dR_top_matched_gluon_medpt", 20, 0, 10);
  book<TH1F>("rest_dphi_top_matched_gluon_medpt", "rest_dphi_top_matched_gluon_medpt", 24, 0, 8);
  book<TH1F>("rest_deta_top_matched_gluon_medpt", "rest_deta_top_matched_gluon_medpt", 24, -0, 8);

  book<TH1F>("rest_dR_top_matched_gluon_highpt", "rest_dR_top_matched_gluon_highpt", 20, 0, 10);
  book<TH1F>("rest_dphi_top_matched_gluon_highpt", "rest_dphi_top_matched_gluon_highpt", 24, 0, 8);
  book<TH1F>("rest_deta_top_matched_gluon_highpt", "rest_deta_top_matched_gluon_highpt", 24, 0, 8);

  book<TH1F>("rest_costheta_mother_gluon", "rest_costheta_mother_gluon", 24, -1, 1);
  book<TH1F>("rest_costheta_lowpt_mother_gluon", "rest_costheta_lowpt_mother_gluon", 24, -1, 1);
  book<TH1F>("rest_costheta_medpt_mother_gluon", "rest_costheta_medpt_mother_gluon", 24, -1, 1);
  book<TH1F>("rest_costheta_highpt_mother_gluon", "rest_costheta_highpt_mother_gluon", 24, -1, 1);

  book<TH1F>("rest_costheta_top_gluon", "rest_costheta_top_gluon", 24, -1, 1);
  book<TH1F>("rest_costheta_lowpt_top_gluon", "rest_costheta_lowpt_top_gluon", 24, -1, 1);
  book<TH1F>("rest_costheta_medpt_top_gluon", "rest_costheta_medpt_top_gluon", 24, -1, 1);
  book<TH1F>("rest_costheta_highpt_top_gluon", "rest_costheta_highpt_top_gluon", 24, -1, 1);

  // CS frame
  book<TH1F>("CS_dR_top_matched_gluon", "CS_dR_top_matched_gluon", 20, 0, 10);
  book<TH1F>("CS_dphi_top_matched_gluon", "CS_dphi_top_matched_gluon", 24, 0, 8);
  book<TH1F>("CS_deta_top_matched_gluon", "CS_deta_top_matched_gluon", 24, 0, 8);

  book<TH1F>("CS_dR_top_matched_gluon_lowpt", "CS_dR_top_matched_gluon_lowpt", 20, 0, 10);
  book<TH1F>("CS_dphi_top_matched_gluon_lowpt", "CS_dphi_top_matched_gluon_lowpt", 24, 0, 8);
  book<TH1F>("CS_deta_top_matched_gluon_lowpt", "CS_deta_top_matched_gluon_lowpt", 24, 0, 8);

  book<TH1F>("CS_dR_top_matched_gluon_medpt", "CS_dR_top_matched_gluon_medpt", 20, 0, 10);
  book<TH1F>("CS_dphi_top_matched_gluon_medpt", "CS_dphi_top_matched_gluon_medpt", 24, 0, 8);
  book<TH1F>("CS_deta_top_matched_gluon_medpt", "CS_deta_top_matched_gluon_medpt", 24, -0, 8);

  book<TH1F>("CS_dR_top_matched_gluon_highpt", "CS_dR_top_matched_gluon_highpt", 20, 0, 10);
  book<TH1F>("CS_dphi_top_matched_gluon_highpt", "CS_dphi_top_matched_gluon_highpt", 24, 0, 8);
  book<TH1F>("CS_deta_top_matched_gluon_highpt", "CS_deta_top_matched_gluon_highpt", 24, 0, 8);

  book<TH1F>("CS_costheta_mother_gluon", "CS_costheta_mother_gluon", 24, -1, 1);
  book<TH1F>("CS_costheta_lowpt_mother_gluon", "CS_costheta_lowpt_mother_gluon", 24, -1, 1);
  book<TH1F>("CS_costheta_medpt_mother_gluon", "CS_costheta_medpt_mother_gluon", 24, -1, 1);
  book<TH1F>("CS_costheta_highpt_mother_gluon", "CS_costheta_highpt_mother_gluon", 24, -1, 1);

  book<TH1F>("CS_costheta_top_gluon", "CS_costheta_top_gluon", 24, -1, 1);
  book<TH1F>("CS_costheta_lowpt_top_gluon", "CS_costheta_lowpt_top_gluon", 24, -1, 1);
  book<TH1F>("CS_costheta_medpt_top_gluon", "CS_costheta_medpt_top_gluon", 24, -1, 1);
  book<TH1F>("CS_costheta_highpt_top_gluon", "CS_costheta_highpt_top_gluon", 24, -1, 1);

  is_mc = ctx.get("dataset_type") == "MC";

}


void TstarTstarAllGenHists::fill(const Event & event){
  if(!is_mc) return;
  assert(event.genparticles);
  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  // Find signal, top, and antitop
  GenParticle top, antitop;
  std::vector<GenParticle> gluons;
  // Find tops
  bool found_top = false, found_antitop = false;
  for(const GenParticle & gp : *event.genparticles){
    if(gp.pdgId() == 6){
      top = gp;
      if(found_top){
        std::cout << "Found multiple tops" << endl;
        return;
      }
      found_top = true;
    }
    else if(gp.pdgId() == -6){
      antitop = gp;
      if(found_antitop){
        std::cout << "Found multiple antitops" << endl;
        return;
      }
      found_antitop = true;
    }
    else if(gp.pdgId() == 21 && gp.status()==23){//only gluons from Tstar decay
      gluons.push_back(gp);
    }
  }

  if(!found_top || !found_antitop) return;

  //////////////////////////
  // The interesting part //
  //////////////////////////

  // matching
  std::vector<std::vector<GenParticle>> sets;
  for (const auto gluon : gluons) {
    std::vector<GenParticle> set;
    if (gluon.mother1() == top.mother1()) {
      set.push_back(event.genparticles->at(gluon.mother1()));
      set.push_back(top);
      set.push_back(gluon);
    }
    else if (gluon.mother1() == antitop.mother1()) {
      set.push_back(event.genparticles->at(gluon.mother1()));
      set.push_back(antitop);
      set.push_back(gluon);
    }
    if(set.size() == 3) sets.push_back(set);
    else if (set.size() != 0) { std::cout << "Set not size 3" << endl; return; }
  }

  int set_count = sets.size();
  hist("matches_found")->Fill(set_count, weight);

  for (const auto & set : sets) {
    std::vector<TString> systems = {"lab", "rest", "CS"};
    for (const auto & system : systems){
      LorentzVector top_v4, gluon_v4, mother_v4;
      if(system == "lab") {
        top_v4 = set.at(1).v4();
        gluon_v4 = set.at(2).v4();
        mother_v4 = top_v4 + gluon_v4;
      }
      else if(system == "rest") {
        top_v4 = set.at(1).v4();
        gluon_v4 = set.at(2).v4();
        mother_v4 = top_v4 + gluon_v4;

        auto boostvec_top = top_v4.BoostToCM(gluon_v4);
        auto boostvec_gluon = top_v4.BoostToCM(top_v4);
        ROOT::Math::Boost boost_top = ROOT::Math::Boost(boostvec_top);
        ROOT::Math::Boost boost_gluon = ROOT::Math::Boost(boostvec_gluon);
        top_v4 = boost_top(top_v4);
        gluon_v4 = boost_gluon(gluon_v4);
      }
      else if(system == "CS") {
        mother_v4 = set.at(1).v4() + set.at(2).v4();
        double mother_pt = mother_v4.pt();
        double qsquared = mother_v4.Dot(mother_v4);
        double q = std::sqrt(qsquared);
        double xt = std::sqrt(qsquared + mother_pt*mother_pt);
        ROOT::Math::LorentzRotation matrix = ROOT::Math::LorentzRotation(
          -mother_pt/q,       xt/q,                                 0,      0,
          0,                  0,                                    1,      0,
          -mother_v4.E()/q,   (mother_pt*mother_v4.px())/(q*xt),    0,      mother_v4.px()/xt,
          mother_v4.px()/q,   -(mother_pt*mother_v4.px())/(q*xt),   0,      -mother_v4.E()/xt
        );
        gluon_v4 = matrix(set.at(2).v4());
        top_v4 = matrix(set.at(1).v4());
        //mother_v4 = top_v4 + gluon_v4;
      }

      // calculating values
      double dR = deltaR(top_v4, gluon_v4);
      double deta = abs(top_v4.eta() - gluon_v4.eta());
      double dphi = abs(top_v4.phi() - gluon_v4.phi());
      double costheta_mother_gluon = ROOT::Math::VectorUtil::CosTheta(mother_v4,gluon_v4);
      double costheta_top_gluon = ROOT::Math::VectorUtil::CosTheta(top_v4,gluon_v4);
      double mother_pt = mother_v4.pt();

      // plotting
      // fill total Hists
      hist(system+"_pt_mother")->Fill(mother_pt, weight/set_count);
      hist(system+"_dR_top_matched_gluon")->Fill(dR, weight/set_count);
      hist(system+"_deta_top_matched_gluon")->Fill(deta, weight/set_count);
      hist(system+"_dphi_top_matched_gluon")->Fill(dphi, weight/set_count);
      hist(system+"_costheta_top_gluon")->Fill(costheta_top_gluon, weight/set_count);
      hist(system+"_costheta_mother_gluon")->Fill(costheta_mother_gluon, weight/set_count);

      // differentiate in mother pt
      if(mother_pt < 400){
        hist(system+"_dR_top_matched_gluon_lowpt")->Fill(dR, weight/set_count);
        hist(system+"_deta_top_matched_gluon_lowpt")->Fill(deta, weight/set_count);
        hist(system+"_dphi_top_matched_gluon_lowpt")->Fill(dphi, weight/set_count);
        hist(system+"_costheta_lowpt_top_gluon")->Fill(costheta_top_gluon, weight/set_count);
        hist(system+"_costheta_lowpt_mother_gluon")->Fill(costheta_mother_gluon, weight/set_count);
      }
      else if(mother_pt < 800){
        hist(system+"_dR_top_matched_gluon_medpt")->Fill(dR, weight/set_count);
        hist(system+"_deta_top_matched_gluon_medpt")->Fill(deta, weight/set_count);
        hist(system+"_dphi_top_matched_gluon_medpt")->Fill(dphi, weight/set_count);
        hist(system+"_costheta_lowpt_top_gluon")->Fill(costheta_top_gluon, weight/set_count);
        hist(system+"_costheta_lowpt_mother_gluon")->Fill(costheta_mother_gluon, weight/set_count);
      }
      else if(mother_pt >= 800){
        hist(system+"_dR_top_matched_gluon_highpt")->Fill(dR, weight/set_count);
        hist(system+"_deta_top_matched_gluon_highpt")->Fill(deta, weight/set_count);
        hist(system+"_dphi_top_matched_gluon_highpt")->Fill(dphi, weight/set_count);
        hist(system+"_costheta_lowpt_top_gluon")->Fill(costheta_top_gluon, weight/set_count);
        hist(system+"_costheta_lowpt_mother_gluon")->Fill(costheta_mother_gluon, weight/set_count);
      }
    }
  }
}

TstarTstarAllGenHists::~TstarTstarAllGenHists(){}
