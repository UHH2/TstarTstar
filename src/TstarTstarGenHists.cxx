#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
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


TstarTstarGenHists::TstarTstarGenHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here

  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

  book<TH1F>("weight", "weight", 40, 0, 2);

  book<TH1F>("x_g1", "x_g1", 20, 0, 1);
  book<TH1F>("x_g2", "x_g2", 20, 0, 1);

  book<TH1F>("Tstar_spin", "Tstar_Spin", 30, -3, 3);

  book<TH1F>("dR_top_gluon", "dR_top_gluon", 24, 0, 6);

  book<TH1F>("gluon_Efrac_AK4", "gluon_Efrac_AK4", 20, 0, 2);
  book<TH1F>("gluon_Efrac_HOTVR", "gluon_Efrac_HOTVR", 20, 0, 2);

  book<TH1F>("dR_hadtop_jet1", "dR_{t_{had} jet 1}", 30, 0, 3);
  book<TH1F>("dR_hadtop_jet2", "dR_{t_{had} jet 2}", 30, 0, 3);
  book<TH1F>("dR_hadtop_jet3", "dR_{t_{had} jet 3}", 30, 0, 3);
  book<TH1F>("dR_bquark_jets", "dR_{b jets}", 30, 0, 3);
  book<TH1F>("dR_gluon_jets", "dR_{gluon jets}", 30, 0, 3);

  book<TH1F>("M_tstar_gen", "M_{Tstar} gen", 100, 0, 2000);
  book<TH1F>("Pt_tstar_gen", "Pt_{Tstar} gen", 100, 0, 3000);
  book<TH1F>("eta_tstar_gen", "Pt_{Tstar} gen", 40, -4, 4);
  book<TH1F>("M_tstartstar_gen", "M_{TstarTstar} gen", 200, 0, 10000);

  const int nbinsX = 31;
  double binsX[nbinsX] = {0, 60, 120, 180, 240, 300, 360, 420, 480, 540, 600, 660, 720, 780, 840, 900, 960, 1020, 1080, 1140, 1200, 1260, 1320, 1380, 1440, 1500, 1620, 1740, 1860, 2100, 3000};
  const int nbinsY = 17;
  double binsY[nbinsY] = {0, 2000, 3000, 3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000, 6500, 8000, 10000};
  book<TH2F>("Pt_M_tstar_2D_gen", "2D pt M gen", nbinsX, binsX, nbinsY, binsY );
  book<TH2F>("Pt_eta_tstar_2D_gen", "2D pt eta gen", 100, 0, 3000, 40, -4, 4);

  book<TH1F>("M_ttbar_gen", "M_{ttbar} gen", 100, 0, 5000);
  book<TH1F>("M_top_gen", "M_{top} gen", 50, 0, 500);
  book<TH1F>("Pt_top_gen", "p^{t}_{T} [GeV/c]", 50, 0, 2500);
  book<TH1F>("Pt_top_gen_oldbins", "p^{t}_{T} [GeV/c]", 100, 0, 2500);
  book<TH1F>("dR_tstartstar_gen", "dR(Tstar,Tstar)", 30, 0, 6);
  book<TH1F>("dR_ttbar_gen", "dR(top,top)", 30, 0, 6);


  book<TH1F>("Pt_gluon1_gen", "Pt_{gluon1} gen", 100, 0, 3000);
  book<TH1F>("Pt_gluon2_gen", "Pt_{gluon2} gen", 100, 0, 3000);
  book<TH1F>("Eta_gluon1_gen", "#eta_{gluon1} gen", 100, -5.2, 5.2);
  book<TH1F>("Eta_gluon2_gen", "#eta_{gluon2} gen", 100, -5.2, 5.2);
  book<TH1F>("Phi_gluon1_gen", "#phi_{gluon1} gen", 100, -3.14, 3.14);
  book<TH1F>("Phi_gluon2_gen", "#phi_{gluon2} gen", 100, -3.14, 3.14);

  book<TH1F>("Pt_photon1_gen", "Pt_{#gamma1} gen", 100, 0, 3000);
  book<TH1F>("Pt_photon2_gen", "Pt_{#gamma2} gen", 100, 0, 3000);
  book<TH1F>("Eta_photon1_gen", "#eta_{#gamma1} gen", 100, -5.2, 5.2);
  book<TH1F>("Eta_photon2_gen", "#eta_{#gamma2} gen", 100, -5.2, 5.2);
  book<TH1F>("Phi_photon1_gen", "#phi_{#gamma1} gen", 100, -3.14, 3.14);
  book<TH1F>("Phi_photon2_gen", "#phi_{#gamma2} gen", 100, -3.14, 3.14);

  book<TH1F>("N_photon_gen", "N_{photons} from TstarTstar", 10, 0, 10);
  book<TH1F>("N_gluon_gen", "N_{gluons} from TstarTstar", 10, 0, 10);
  book<TH1F>("dR_tophad_gluon_gen", "dR (tophad, gluon)", 30, 0, 6);
  book<TH1F>("dR_toplep_gluon_gen", "dR (toplep, gluon)", 30, 0, 6);
  book<TH1F>("dR_min_gluon_top_gen", "dR^{min}(gluon,top)", 30, 0, 6);
  book<TH1F>("dR_min_photon_top_gen", "dR^{min}(#gamma,top)", 30, 0, 6);
  book<TH1F>("dR_photon_gluon_gen", "dR(#gamma,gluon) (only for Tgamma+Tgluon)", 30, 0, 6);

  book<TH1F>("dR_hadtop_b_q1", "dR(b, q1)", 30, 0, 6);
  book<TH1F>("dR_hadtop_b_q2", "dR(b, q2)", 30, 0, 6);
  book<TH1F>("dR_hadtop_q1_q2", "dR(q1, q2)", 30, 0, 6);
  book<TH1F>("max_dR_hadtop", "max dR(t_{had} decay products)", 30, 0, 6);
  book<TH2D>("2D_max_dR_hadtop", "pT t_{had} vs max(dR) of t_{had} decays", 100, 0, 3000, 100, 0, 6);

  book<TH1F>("isolated_partons", "number of isolated partons", 7, -1, 6);

  is_mc = ctx.get("dataset_type") == "MC";

}


void TstarTstarGenHists::fill(const Event & event){
  if(!is_mc) return;
  assert(event.genparticles);
  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  hist("weight")->Fill(weight);

  TTbarGen ttbargen = event.get(h_ttbargen);
  if(!ttbargen.IsSemiLeptonicDecay()){return;}

  // Find signal, top, and antitop
  GenParticle tstar, antitstar, top, antitop, gluon1, gluon2, photon1, photon2;
  bool found_tstar = false, found_antitstar = false, found_top = false, found_antitop = false;
  bool found_gluon1 = false, found_gluon2 = false, found_photon1 = false, found_photon2 = false;
  int n_gluons = 0, n_photons = 0;
  for(const GenParticle & gp : *event.genparticles){
    //    cout<<"gp.pdgId() = "<<gp.pdgId()<<" status() = "<<gp.status()<<endl;
    //    if(gp.status()==22) cout<<"gp.pdgId() = "<<gp.pdgId()<<" status() = "<<gp.status()<<endl;
    //    if(gp.pdgId()==21) cout<<"gp.pdgId() = "<<gp.pdgId()<<" status() = "<<gp.status()<<" pt = "<<gp.pt()<<endl;
    if(gp.pdgId() == 6){
      top = gp;
      found_top = true;
    }
    else if(gp.pdgId() == -6){
      antitop = gp;
      found_antitop = true;
    }
    else if(gp.pdgId() == 600 && (gp.status()==23 || gp.status()==22)){
      tstar = gp;
      if(found_tstar) std::cout << "Error: found two T*" << std::endl;
      found_tstar = true;
    }
    else if(gp.pdgId() == -600 && (gp.status()==23 || gp.status()==22)){
      antitstar = gp;
      if(found_antitstar) std::cout << "Error: found two anti T*" << std::endl;
      found_antitstar = true;
    }
    else if(abs(gp.pdgId()) == 21 || abs(gp.pdgId()) == 9){//only gluons from Tstar decay
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
    else if(gp.pdgId() == 22){//only photons from Tstar decay
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

  if(!found_tstar || !found_antitstar) return;
  if(!found_top || !found_antitop) return;

  /**
    if(tstar.daughter(event.genparticles, 1)) std::cout << "Tstar daughter 1 pdgID: " << tstar.daughter(event.genparticles, 1)->pdgId() << std::endl;
    else std::cout << "Tstar daughter 1 not found." << std::endl;
    if(tstar.daughter(event.genparticles, 2)) std::cout << "Tstar daughter 2 pdgID: " << tstar.daughter(event.genparticles, 2)->pdgId() << std::endl;
    else std::cout << "Tstar daughter 2 not found." << std::endl;
    if(antitstar.daughter(event.genparticles, 1)) std::cout << "Anti-Tstar daughter 1 pdgID: " << antitstar.daughter(event.genparticles, 1)->pdgId() << std::endl;
    else std::cout << "Anti-Tstar daughter 1 not found." << std::endl;
    if(antitstar.daughter(event.genparticles, 2)) std::cout << "Anti-Tstar daughter 2 pdgID: " << antitstar.daughter(event.genparticles, 2)->pdgId() << std::endl;
    else std::cout << "Anti-Tstar daughter 2 not found." << std::endl;
  **/

  int isolated_partons = 4;
  std::vector<GenParticle> partons = {top, antitop, gluon1, gluon2};
  for(uint i = 0; i < 4; i++){
    for(uint j = 0; j < 4; j++){
      if(j == i) continue;
      if(partons.at(i).pt() == 0){isolated_partons--;}
      else if(deltaR(partons.at(i), partons.at(j)) < 600/(partons.at(i).pt())){
        isolated_partons--;
        break;
      }
    }
  }
  hist("isolated_partons")->Fill(isolated_partons, weight);

  double dR_tstarstar = deltaR(tstar,antitstar);
  double dR_ttbar = deltaR(top,antitop);
  hist("dR_tstartstar_gen")->Fill(dR_tstarstar, weight);
  hist("dR_ttbar_gen")->Fill(dR_ttbar, weight);
  // book<TH1F>("dR_tstartstar_gen", "dR(Tstar,Tstar)", 30, 0, 6);
  // book<TH1F>("dR_ttbar_gen", "dR(top,top)", 30, 0, 6);

  if(found_tstar){
    float m_tstar = inv_mass(tstar.v4());
    hist("M_tstar_gen")->Fill(m_tstar, weight);
    hist("Pt_tstar_gen")->Fill(tstar.pt(), weight);
    hist("eta_tstar_gen")->Fill(tstar.eta(), weight);
    ((TH2F*)hist("Pt_eta_tstar_2D_gen"))->Fill(tstar.pt(), tstar.eta(), weight);
  }
  if(found_tstar && found_antitstar){
    float m_tstartstar = inv_mass(tstar.v4()+antitstar.v4());
    hist("M_tstartstar_gen")->Fill(m_tstartstar, weight);
    ((TH2F*)hist("Pt_M_tstar_2D_gen"))->Fill(m_tstartstar, tstar.pt(), weight);
  }
  float m_ttbar = inv_mass(top.v4() + antitop.v4());
  hist("M_ttbar_gen")->Fill(m_ttbar, weight);
  hist("Pt_top_gen")->Fill(top.pt(), weight/2);
  hist("Pt_top_gen")->Fill(antitop.pt(), weight/2);
  hist("Pt_top_gen_oldbins")->Fill(top.pt(), weight/2);
  hist("Pt_top_gen_oldbins")->Fill(antitop.pt(), weight/2);
  hist("M_top_gen")->Fill(inv_mass(top.v4()), weight);
  hist("M_top_gen")->Fill(inv_mass(antitop.v4()), weight);

  if(found_gluon1){
    hist("Pt_gluon1_gen")->Fill(gluon1.pt(),weight);
    hist("Eta_gluon1_gen")->Fill(gluon1.eta(),weight);
    hist("Phi_gluon1_gen")->Fill(gluon1.phi(),weight);
  }
  if(found_gluon2){
    hist("Pt_gluon2_gen")->Fill(gluon2.pt(),weight);
    hist("Eta_gluon2_gen")->Fill(gluon2.eta(),weight);
    hist("Phi_gluon2_gen")->Fill(gluon2.phi(),weight);
  }
  if(found_photon1){
    hist("Pt_photon1_gen")->Fill(photon1.pt(),weight);
    hist("Eta_photon1_gen")->Fill(photon1.eta(),weight);
    hist("Phi_photon1_gen")->Fill(photon1.phi(),weight);
  }
  if(found_photon2){
    hist("Pt_photon2_gen")->Fill(photon2.pt(),weight);
    hist("Eta_photon2_gen")->Fill(photon2.eta(),weight);
    hist("Phi_photon2_gen")->Fill(photon2.phi(),weight);
  }
  hist("N_gluon_gen")->Fill(n_gluons,weight);
  hist("N_photon_gen")->Fill(n_photons,weight);

  if(found_gluon1 && found_gluon2){
    double dR11 = deltaR(gluon1, top);
    double dR12 = deltaR(gluon1, antitop);
    double dR21 = deltaR(gluon2, top);
    double dR22 = deltaR(gluon2, antitop);
    double dR_min1 = dR11;
    if(dR12<dR11) dR_min1 = dR12;
    double dR_min2 = dR21;
    if(dR22<dR21) dR_min2 = dR22;
    //    cout<<"dR11, dR12, dR21, dR22 = "<<dR11<<", "<<dR12<<", "<<dR21<<", "<<dR22<<endl;
    //    cout<<"dR_min1, dR_min2 = "<<dR_min1<<", "<<dR_min2<<endl;
    hist("dR_min_gluon_top_gen")->Fill(dR_min1,weight);
    hist("dR_min_gluon_top_gen")->Fill(dR_min2,weight);
  }
  if(found_photon1 && found_photon2){
    double dR11 = deltaR(photon1, top);
    double dR12 = deltaR(photon1, antitop);
    double dR21 = deltaR(photon2, top);
    double dR22 = deltaR(photon2, antitop);
    double dR_min1 = dR11;
    if(dR12<dR11) dR_min1 = dR12;
    double dR_min2 = dR21;
    if(dR22<dR21) dR_min2 = dR22;
    hist("dR_min_photon_top_gen")->Fill(dR_min1,weight);
    hist("dR_min_photon_top_gen")->Fill(dR_min2,weight);
  }

  if(((found_gluon1 && !found_gluon2) || (found_gluon2 && !found_gluon2)) && ((found_photon1 && !found_photon2) || (!found_photon1 && found_photon2))){//for Tgluon+Tgamma
    GenParticle gluon, photon;
    if(found_gluon1 && !found_gluon2) gluon = gluon1;
    if(!found_gluon1 && found_gluon2) gluon = gluon2;
    if(found_photon1 && !found_photon2) photon = photon1;
    if(!found_photon1 && found_photon2) photon = photon2;
    double dR11 = deltaR(gluon, top);
    double dR12 = deltaR(gluon, antitop);
    double dR21 = deltaR(photon, top);
    double dR22 = deltaR(photon, antitop);
    double dR_min1 = dR11;
    if(dR12<dR11) dR_min1 = dR12;
    hist("dR_min_gluon_top_gen")->Fill(dR_min1,weight);
    double dR_min2 = dR21;
    if(dR22<dR21) dR_min2 = dR22;
    hist("dR_min_photon_top_gen")->Fill(dR_min2,weight);

    //Plot this only for the interesting Tgluon+Tgamma case
    /**
    double dRgluongamma = deltaR(gluon, photon);
    hist("dR_photon_gluon_gen")->Fill(dRgluongamma, weight);
    **/
  }

  GenParticle toplep = ttbargen.TopLep();
  GenParticle tophad = ttbargen.TopHad();

  GenParticle gluonlep, gluonhad;
  if(toplep.mother1() ==  gluon1.mother1()){
    gluonlep = gluon1;
    gluonhad = gluon2;
  }
  else{
    gluonlep = gluon2;
    gluonhad = gluon1;
  }

  double dR_toplep_gluon = deltaR(toplep, gluonlep);
  double dR_tophad_gluon = deltaR(tophad, gluonhad);
  hist("dR_toplep_gluon_gen")->Fill(dR_toplep_gluon,weight);
  hist("dR_tophad_gluon_gen")->Fill(dR_tophad_gluon,weight);

  double dR_hadtop_b_q1 = deltaR(ttbargen.BHad(), ttbargen.Q1());
  double dR_hadtop_b_q2 = deltaR(ttbargen.BHad(), ttbargen.Q2());
  double dR_hadtop_q1_q2 = deltaR(ttbargen.Q1(), ttbargen.Q2());
  hist("dR_hadtop_b_q1")->Fill(dR_hadtop_b_q1, weight);
  hist("dR_hadtop_b_q2")->Fill(dR_hadtop_b_q2, weight);
  hist("dR_hadtop_q1_q2")->Fill(dR_hadtop_q1_q2, weight);
  hist("max_dR_hadtop")->Fill(max({dR_hadtop_b_q1, dR_hadtop_b_q2, dR_hadtop_q1_q2}), weight);
  ((TH2D*)hist("2D_max_dR_hadtop"))->Fill(tophad.pt(), max({dR_hadtop_b_q1, dR_hadtop_b_q2, dR_hadtop_q1_q2}), weight);

  GenInfo* genInfo = event.genInfo;

  hist("x_g1")->Fill(genInfo->pdf_x1(), weight);
  hist("x_g2")->Fill(genInfo->pdf_x2(), weight);

  std::vector<GenParticle> gluons;
  std::vector<GenParticle> tops;
  int ngluons = 0;
  for(const GenParticle & gp : *event.genparticles){
    if(gp.pdgId() == 600 && (gp.status()==23 || gp.status()==22)){
      hist("Tstar_spin")->Fill(gp.spin(), weight);
    }
    else if(gp.pdgId() == -600 && (gp.status()==23 || gp.status()==22)){
      hist("Tstar_spin")->Fill(gp.spin(), weight);
    }
      else if(gp.pdgId() == 21 && (gp.status()==23 || gp.status()==22)){
      gluons.push_back(gp);
      ngluons++;
    }
    else if((gp.pdgId() == 6) || (gp.pdgId() == -6))
      tops.push_back(gp);
  }

  int combsFound = 0;
  for(const auto & top : tops){
    for(const auto & gluon : gluons){
      if(top.mother1() == gluon.mother1()){
      	combsFound++;
      	hist("dR_top_gluon")->Fill(deltaR(top, gluon), weight);
      }
    }
  }

  if(combsFound != 2){
    //cout << "Error: not exactly two top gluon sibling pairs found. We have " << combsFound << "!" << endl;
  }

  // Check energy in jets corresponding to gluons
  for(const auto & gluon : gluons){
    bool matchedAK4 = false;
    double bestmatchAK4 = 999;
    Jet matchedAK4jet;
    for(const auto & jet : *event.jets){
      double dR = deltaR(gluon, jet);
      if(dR < 0.4 && dR < bestmatchAK4){
      	matchedAK4 = true;
      	bestmatchAK4 = dR;
      	matchedAK4jet = jet;
      }
    } // found (best) matched AK4 jet

    bool matchedAK8 = false;
    double bestmatchAK8 = 999;
    Jet matchedAK8jet;
    for(const auto & jet : *event.topjets){
      double R_jet = (jet.pt() > 0) ? 600/(jet.pt()) : 1.5;
      if(R_jet > 1.5) R_jet = 1.5;
      if(R_jet < 0.1) R_jet = 0.1;
      double dR = deltaR(gluon, jet);
      if(dR < R_jet && dR < bestmatchAK8){
      	matchedAK8 = true;
      	bestmatchAK8 = dR;
      	matchedAK8jet = jet;
      }
    } // found (best) matched AK4 jet

    double gluon_E = gluon.energy();
    if(matchedAK4 && matchedAK8){
      hist("gluon_Efrac_AK4")->Fill(matchedAK4jet.energy()/gluon_E, weight);
      hist("gluon_Efrac_HOTVR")->Fill(matchedAK8jet.energy()/gluon_E, weight);
    }

    hist("dR_gluon_jets")->Fill(deltaR(matchedAK4jet, matchedAK8jet), weight);

  }

  // Find hadronic top jet
  double dR_hadtopjet = 999;
  TopJet hadtopjet;
  for(const auto & jet : *event.topjets){
    double dR = deltaR(jet, ttbargen.TopHad());
    if(dR < dR_hadtopjet){
      dR_hadtopjet = dR;
      hadtopjet = jet;
    }
  }

  // Find leptonic b jet
  double dR_lepbjet = 999;
  Jet lepbjet;
  for(const auto & jet : *event.jets){
    //if(deltaR(jet, hadtopjet) < 1.2) continue;
    double dR = deltaR(jet, ttbargen.BLep());
    if(dR < dR_lepbjet){
      dR_lepbjet = dR;
      lepbjet = jet;
    }
  }

  // Find hadronic slim jet
  double dR_hadtopjet_slim = 999;
  Jet hadtopjet_slim1;
  Jet hadtopjet_slim2;
  Jet hadtopjet_slim3;
  int count_found = 0;
  for(const auto & jet : *event.jets){
    double dR = deltaR(jet, ttbargen.TopHad());
    if(dR < dR_hadtopjet_slim){
      count_found++;
      dR_hadtopjet_slim = dR;
      hadtopjet_slim3 = hadtopjet_slim2;
      hadtopjet_slim2 = hadtopjet_slim1;
      hadtopjet_slim1 = jet;
    }
  }
  if(count_found > 0){
    hist("dR_hadtop_jet1")->Fill(deltaR(hadtopjet, hadtopjet_slim1), weight);
    if(count_found > 1){
      hist("dR_hadtop_jet2")->Fill(deltaR(hadtopjet, hadtopjet_slim2), weight);
      if(count_found > 2){
	hist("dR_hadtop_jet3")->Fill(deltaR(hadtopjet, hadtopjet_slim3), weight);
      }
    }
  }

  // Find leptonic b jet
  double dR_lepbjet_fat = 999;
  TopJet lepbjet_fat;
  for(const auto & jet : *event.topjets){
    double dR = deltaR(jet, ttbargen.BLep());
    if(dR < dR_lepbjet_fat){
      dR_lepbjet_fat = dR;
      lepbjet_fat = jet;
    }
  }
  hist("dR_bquark_jets")->Fill(deltaR(lepbjet, lepbjet_fat), weight);


}


TstarTstarGenHists::~TstarTstarGenHists(){}



// #########################







TstarTstarMergedHists::TstarTstarMergedHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  book<TH1F>("Merge_leading_gluon", "Merge_g1 gen", 4, 0, 4);
  book<TH1F>("Merge_second_gluon", "Merge_g2 gen", 4, 0, 4);
  book<TH1F>("Merge_TstarLep", "Merge_TstarLep gen", 4, 0, 4);
  book<TH1F>("Merge_TstarHad", "Merge_TstarHad gen", 4, 0, 4);
  book<TH1F>("Merge_tLep", "Merge_tLep gen", 4, 0, 4);
  book<TH1F>("Merge_tHad", "Merge_tHad gen", 4, 0, 4);

  bool is_tgtg = false;
  if(ctx.get("channel") == "tgtg") is_tgtg = true;

  if(is_tgtg) h_recohyp_tstartstar_best_ = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_tgtg_best");

  is_mc = ctx.get("dataset_type") == "MC";

}

void TstarTstarMergedHists::fill(const Event & event){

  // HELP CLASS FOR THE MOMENT ONLY VALID FOR TGTG!!!!
  //if(is_tgtgamma){return false;}

  if(!is_mc) return;
  assert(event.genparticles);
  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  // Find signal, top, and antitop
  GenParticle tstar, antitstar, top, antitop, gluon1, gluon2;
  bool found_tstar = false, found_antitstar = false, found_top = false, found_antitop = false;
  bool found_gluon1 = false, found_gluon2 = false;
  for(const GenParticle & gp : *event.genparticles){
    if(gp.pdgId() == 6){
      top = gp;
      found_top = true;
    }
    else if(gp.pdgId() == -6){
      antitop = gp;
      found_antitop = true;
    }
    else if(gp.pdgId() == 600 && (gp.status()==23 || gp.status()==22)){
      tstar = gp;
      found_tstar = true;
    }
    else if(gp.pdgId() == -600 && (gp.status()==23 || gp.status()==22)){
      antitstar = gp;
      found_antitstar = true;
    }
    else if(gp.pdgId() == 21 && gp.status()==23){//only gluons from Tstar decay
      if(!found_gluon1){
      	gluon1 = gp;
      	found_gluon1 = true;
      }
      else{
    	gluon2 = gp;
    	found_gluon2 = true;
      }
    }
  }

  if(!(found_tstar && found_antitstar && found_top && found_antitop && found_gluon1 && found_gluon2)) {return;}

  ReconstructionTstarHypothesis tstar_hyp_best = event.get(h_recohyp_tstartstar_best_);

  double fill = 0;
  if(deltaR(gluon1.v4(), tstar_hyp_best.gluon1_v4()) < 0.8){fill+=1;}
  if(deltaR(gluon2.v4(), tstar_hyp_best.gluon1_v4()) < 0.8){fill+=1;}
  hist("Merge_leading_gluon")->Fill(fill +0.5);

  fill = 0;
  if(deltaR(gluon1.v4(), tstar_hyp_best.gluon2_v4()) < 0.8){fill+=1;}
  if(deltaR(gluon2.v4(), tstar_hyp_best.gluon2_v4()) < 0.8){fill+=1;}
  hist("Merge_second_gluon")->Fill(fill+0.5);

  fill = 0;
  if(deltaR(tstar.v4(), tstar_hyp_best.tstarlep_v4()) < 0.8){fill+=1;}
  if(deltaR(antitstar.v4(), tstar_hyp_best.tstarlep_v4()) < 0.8){fill+=1;}
  hist("Merge_TstarLep")->Fill(fill+0.5);

  fill = 0;
  if(deltaR(tstar.v4(), tstar_hyp_best.tstarhad_v4()) < 0.8){fill+=1;}
  if(deltaR(antitstar.v4(), tstar_hyp_best.tstarhad_v4()) < 0.8){fill+=1;}
  hist("Merge_TstarHad")->Fill(fill+0.5);

  fill = 0;
  if(deltaR(top.v4(), tstar_hyp_best.ttbar_hyp().toplep_v4()) < 0.8){fill+=1;}
  if(deltaR(antitop.v4(), tstar_hyp_best.ttbar_hyp().toplep_v4()) < 0.8){fill+=1;}
  hist("Merge_tLep")->Fill(fill+0.5);

  fill = 0;
  if(deltaR(top.v4(), tstar_hyp_best.ttbar_hyp().tophad_v4()) < 0.8){fill+=1;}
  if(deltaR(antitop.v4(), tstar_hyp_best.ttbar_hyp().tophad_v4()) < 0.8){fill+=1;}
  hist("Merge_tHad")->Fill(fill+0.5);

}

//TstarTstarMergedHists::~TstarTstarMergedHists(){}
