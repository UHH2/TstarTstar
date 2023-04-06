#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/core/include/Event.h"

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <memory>
//using namespace uhh2examples;
using namespace uhh2;
using namespace std;

//From Z'->ttbar
bool uhh2::TwoDCut::passes(const uhh2::Event& event){

  assert(event.muons || event.electrons);
  //  if((event.muons->size()+event.electrons->size()) != 1) throw runtime_error("In TwoDCut: Event does not have exactly one muon and electron.");
  if((event.muons->size()+event.electrons->size()) <1) throw runtime_error("In TwoDCut: Event does not have at least one muon and electron.");

  bool pass = false;
  if(event.muons->size()){
    float drmin, ptrel;
    drmin = event.muons->at(0).get_tag(Muon::twodcut_dRmin);
    ptrel = event.muons->at(0).get_tag(Muon::twodcut_pTrel);
    if((drmin > min_deltaR_) || (ptrel > min_pTrel_)) pass = true;
  }
  if(event.electrons->size()){
    float drmin, ptrel;
    drmin = event.electrons->at(0).get_tag(Electron::twodcut_dRmin);
    ptrel = event.electrons->at(0).get_tag(Electron::twodcut_pTrel);
    if((drmin > min_deltaR_) || (ptrel > min_pTrel_)) pass = true;
  }

  return pass;
}

//From Z'->ttbar
TTbarSemiLepMatchableSelection::TTbarSemiLepMatchableSelection(){
  Wlep = GenParticle(); Whad = GenParticle();
  blep =  GenParticle(); bhad = GenParticle();
  thad =  GenParticle(); tlep =  GenParticle();
  lepton =  GenParticle(); neutrino =  GenParticle(); 
  Whadd1 =  GenParticle(); Whadd2 =  GenParticle();
}
bool TTbarSemiLepMatchableSelection::passes(const Event & event){
  // For pdgID mapping see https://twiki.cern.ch/twiki/bin/view/Main/PdgId
  if(event.isRealData) return false;
  assert(event.genparticles);

  //check, if one top decays had and one decays lep
  bool found_had = false, found_lep = false;

  //Loop over genparticles
  for(const auto & gp : *event.genparticles){
    //    cout<<" gp.pdgId() = "<<gp.pdgId()<<" status() = "<<gp.status()<<endl;
    //Get tops
    if(fabs(gp.pdgId()) == 6){

      //Get b and W
      auto b = gp.daughter(event.genparticles,1);
      auto W = gp.daughter(event.genparticles,2);
      if(fabs(W->pdgId()) == 5 && fabs(b->pdgId()) == 24){
        b = gp.daughter(event.genparticles,2);
        W = gp.daughter(event.genparticles,1);
      }

      if(abs(W->pdgId()) != 24) {
        for(unsigned int j = 0; j < event.genparticles->size(); ++j) {
          const GenParticle & genp = event.genparticles->at(j);
          auto m1 = genp.mother(event.genparticles, 1);
          auto m2 = genp.mother(event.genparticles, 2);
          bool has_top_mother = ((m1 && m1->index() == gp.index()) || (m2 && m2->index() == gp.index()));
          if(has_top_mother && (abs(genp.pdgId()) == 24)) {
            W = &genp;
            break;
          }
        }
      }
      if(abs(b->pdgId()) != 5 && abs(b->pdgId()) != 3 && abs(b->pdgId()) != 1) {//b,s,d
        for(unsigned int j = 0; j < event.genparticles->size(); ++j) {
          const GenParticle & genp = event.genparticles->at(j);
          auto m1 = genp.mother(event.genparticles, 1);
          auto m2 = genp.mother(event.genparticles, 2);
          bool has_top_mother = ((m1 && m1->index() == gp.index()) || (m2 && m2->index() == gp.index()));
          if(has_top_mother && (abs(genp.pdgId()) == 5 || abs(genp.pdgId()) == 3 || abs(genp.pdgId()) == 1)) {
            b = &genp;
            break;
          }
        }
      }
      //      cout<<" W->pdgId() = "<<W->pdgId()<<" b->pdgId() = "<<b->pdgId()<<endl;
      //      if(!((fabs(b->pdgId()) == 5 || fabs(b->pdgId()) == 3 || fabs(b->pdgId()) == 1) && fabs(W->pdgId()) == 24)) cout<<"!!! Top decay products (W,b) are not found"<<endl;
      if(!((fabs(b->pdgId()) == 5 || fabs(b->pdgId()) == 3 || fabs(b->pdgId()) == 1) && fabs(W->pdgId()) == 24)) return false;

      //To identify decay type, check ID of W daughters
      auto Wd1 = W->daughter(event.genparticles,1);
      auto Wd2 = W->daughter(event.genparticles,2);
      //      cout<<"  Wd1 ID = "<<fabs(Wd1->pdgId())<<" Wd2 ID = "<<fabs(Wd2->pdgId())<<endl;

      //try to match the b quarks
      bool matched_b_ak4 = false;

      // Consider AK4 jets first
      for(const auto & jet : *event.jets){
	//	cout<<"   Distance between b and AK4 jets:"<<deltaR(*b,jet)<<endl;
        if(deltaR(*b,jet) <= 0.4) matched_b_ak4 = true;
      }

      bool matched_b_ak8 = false;
      // Now consider AK8 jets
      int idx_matched_topjet = -1;
      int idx = 0;
      for(const auto & jet : *event.topjets){
	//	cout<<"    Distance between b and AK8 jets:"<<deltaR(*b,jet)<<endl;
        if(deltaR(*b,jet) <= 0.8){
          matched_b_ak8 = true;
          idx_matched_topjet = idx;
        }
        idx++;
      }
      //      if(!matched_b_ak4 && !matched_b_ak8) cout<<"!!! b jet is not matched to neither AK4 neither to AK8, eta of b gen = "<<b->eta()<<", pt of b gen = "<<b->pt()<<endl;
      if(!matched_b_ak4 && !matched_b_ak8) return false;

      //Check decaymodes of W
    
      //hadronic
      if(fabs(Wd1->pdgId()) < 7 && fabs(Wd2->pdgId()) < 7){
        if(found_had) return false;
        found_had = true;
	Whad = *W;
	bhad = *b;
	thad = gp;
	Whadd1 = *Wd1;
	Whadd2 = *Wd2;
        //check if both daughters can be matched by jets
        bool matched_d1_ak4 = false, matched_d2_ak4 = false;
        bool matched_d1_ak8 = false, matched_d2_ak8 = false;
        // Consider AK4 jets first
        for(const auto & jet : *event.jets){
          if(deltaR(*Wd1, jet) <= 0.4) matched_d1_ak4 = true;
          if(deltaR(*Wd2, jet) <= 0.4) matched_d2_ak4 = true;
        }

        // Now consider the one AK8 jet also used for the b-jet
        if(!matched_b_ak8){
          matched_d1_ak8 = false;
          matched_d2_ak8 = false;
        }
        else{
          if(deltaR(*Wd1, event.topjets->at(idx_matched_topjet)) <= 0.8) matched_d1_ak8 = true;
          if(deltaR(*Wd2, event.topjets->at(idx_matched_topjet)) <= 0.8) matched_d2_ak8 = true;
        }

        // if(!(matched_d1 && matched_d2)) return false;
	// cout<<" Wd1->eta(), Wd1->pt(): "<<Wd1->eta()<<", "<<Wd1->pt()<<" Wd2->eta(), Wd2->pt(): "<<Wd2->eta()<<", "<<Wd2->pt()<<" dR(Wd1,Wd2) = "<<deltaR(*Wd1,*Wd2)<<endl;
	// cout<<"Easy case: AK4 jets only, match to b, Wd1, Wd2: "<<matched_b_ak4<<", "<<matched_d1_ak4<<", "<<matched_d2_ak4<<endl;
	// cout<<"Easy case: AK8 jets only, match to b, Wd1, Wd2: "<<matched_b_ak8<<", "<<matched_d1_ak8<<", "<<matched_d2_ak8<<endl;

	// if(!(matched_b_ak4 && matched_d1_ak4 && matched_d2_ak4) && !(matched_b_ak8 && matched_d1_ak8 && matched_d2_ak8))
	//   cout<<"!!! Hadronic side is not matched, ak4:"<<matched_b_ak4<<" "<<matched_d1_ak4<<" "<<matched_d2_ak4<<" ak8:"<<matched_b_ak8<<" "<<matched_d1_ak8<<" "<<matched_d2_ak8<<endl;
        if(!(matched_b_ak4 && matched_d1_ak4 && matched_d2_ak4) && !(matched_b_ak8 && matched_d1_ak8 && matched_d2_ak8)) return false;
      }

      //leptonic
      else if((abs(Wd1->pdgId()) == 11 || abs(Wd1->pdgId()) == 13) || (abs(Wd2->pdgId()) == 11 || abs(Wd2->pdgId()) == 13)){
        if(found_lep) return false;

        // Escape cases where the W radiates an intermediate photon, that splits into llbar
        if(Wd1->pdgId() == -Wd2->pdgId()){
          // cout << "Entered the escape-part" << endl;
          // Find 2 genparts with 11,12,13,14 that follow each other in the list and don't have the same fabs
          int idx = 0;
          for(const auto & genp : *event.genparticles){
            if(found_lep) break;
            if(abs(genp.pdgId()) >= 11 && abs(genp.pdgId()) <= 14){
              bool is_charged = (abs(genp.pdgId()) == 11 || abs(genp.pdgId()) == 13);
              // cout << "Found a genpart at index " << idx << " with id " << genp.pdgId() << ", is_charged: " << is_charged << endl;

              // if the first one is charged, the second one has to have pdgId of +1 wrt. this genpart
              if(is_charged){
                // cout << "(charged) Going to check for next particle in list" << endl;
                if(abs(event.genparticles->at(idx+1).pdgId()) == abs(genp.pdgId()) + 1){
                  Wd1 = &genp;
                  Wd2 = &event.genparticles->at(idx+1);
                  found_lep = true;
                }
              }
              else{
                // cout << "(neutral) Going to check for next particle in list" << endl;
                if(abs(event.genparticles->at(idx+1).pdgId()) == abs(genp.pdgId()) - 1){
                  Wd2 = &genp;
                  Wd1 = &event.genparticles->at(idx+1);
                  found_lep = true;
                }
              }
            }
            idx++;
          }
	  //          if(!found_lep) cout<<"!!! Not found leptonic side";
          if(!found_lep) return false;
        }

        found_lep = true;

        //Find charged lepton
        auto lep = Wd1;
        auto nu = Wd2;
        if(fabs(Wd2->pdgId()) == 11 || fabs(Wd2->pdgId()) == 13){
          lep = Wd2;
          nu = Wd1;
        }
        if(!(abs(lep->pdgId()) == 11 && abs(nu->pdgId()) == 12) && !(abs(lep->pdgId()) == 13 && abs(nu->pdgId()) == 14)) throw runtime_error("In TTbarSemiLepMatchable: The leptonic W does not decay into a lepton and its neutrino.");
	Wlep = *W;
	blep = *b;
	tlep = gp;
	lepton = *lep;
	neutrino = *nu;

        //check, if lepton can be matched
        bool matched_lep = false;
        if(fabs(lep->pdgId()) == 11){
          for(const auto & ele : *event.electrons){
            if(deltaR(*lep,ele) <= 0.2) matched_lep = true;
          }
        }
        else if(fabs(lep->pdgId()) == 13){
          for(const auto & mu : *event.muons){
            if(deltaR(mu,*lep) <= 0.2) matched_lep = true;
          }
        }
        else throw runtime_error("In TTbarSemiLepMatchable: Lepton from W decay is neither e nor mu.");
        if(!matched_lep) return false;
      }
      //tau-decays
      else return false;
    }
  }

  //  if(!(found_had && found_lep)) cout<<"!!! found_had = "<<found_had<<" found_lep = "<<found_lep<<endl;
  if(!(found_had && found_lep)) return false;
  //  cout<<"### found_had = "<<found_had<<" found_lep = "<<found_lep<<endl;
  return true;
}

std::pair<bool,double> TTbarSemiLepMatchableSelection::check_reco(const ReconstructionHypothesis hyp){

  //Following check with matching to GEN particles as described in AN2015-107 (Z' with 2015 data)
  /**
  //Hadronic top V1
  bool tophad_match = false;
  double dR_Wd1_min = 1e6;
  double dR_Wd2_min = 1e6;
  double dR_bhad_min = 1e6;
  double dR_other_tophadjets = 0; // to prevent this from strongly favouring hypopthesises with way to many jets, calculate dR for each "false" jet and add that to result.
  if(!hyp.tophad_topjet_ptr()){//hadronic top reconstructed as set of AK4 jets
    for (uint i = 0; i < hyp.tophad_jets().size(); i++){
      bool merged = false;
      double dR_Wd1 = deltaR(hyp.tophad_jets().at(i).v4(), Whadd1.v4());
      double dR_Wd2 = deltaR(hyp.tophad_jets().at(i).v4(), Whadd2.v4());
      double dR_bhad = deltaR(hyp.tophad_jets().at(i).v4(), bhad.v4());
      double min_dR_tmp = 1e6;
      if(dR_Wd1_min>dR_Wd1) { min_dR_tmp = min(min_dR_tmp, dR_Wd1_min); dR_Wd1_min = dR_Wd1; merged = true;}
      if(dR_Wd2_min>dR_Wd2) { min_dR_tmp = min(min_dR_tmp, dR_Wd2_min); dR_Wd2_min = dR_Wd2; merged = true;}
      if(dR_bhad_min>dR_bhad) { min_dR_tmp = min(min_dR_tmp, dR_bhad_min); dR_bhad_min = dR_bhad; merged = true;}
      if(merged && (min_dR_tmp < 1e6)){dR_other_tophadjets += min_dR_tmp;} // This is still wrong. Now something can be merged for two jets, be replaced for one and be counted twice!
      else{dR_other_tophadjets += min({dR_Wd1, dR_Wd2, dR_bhad});}
    }
    if(dR_Wd1_min<0.4 && dR_Wd2_min<0.4 && dR_bhad_min<0.4) tophad_match = true;

  }
  else{//hadronic top reconstructed as AK8 jet
    dR_Wd1_min = deltaR(hyp.tophad_topjet_ptr()->v4(), Whadd1.v4());
    dR_Wd2_min = deltaR(hyp.tophad_topjet_ptr()->v4(), Whadd2.v4());
    dR_bhad_min = deltaR(hyp.tophad_topjet_ptr()->v4(), bhad.v4());
    if(dR_Wd1_min<0.8 && dR_Wd2_min<0.8 && dR_bhad_min<0.8) tophad_match = true;
  }
  **/

  //Hadronic top V2
  bool tophad_match = false;
  double dR_hadtop = 0.;
  if(!hyp.tophad_topjet_ptr()){ //hadronic top reconstructed as set of AK4 jets
    bool Wd1_merged = false;
    bool Wd2_merged = false;
    bool bhad_merged = false;
    for (uint i = 0; i < hyp.tophad_jets().size(); i++){
      double dR_Wd1 = deltaR(hyp.tophad_jets().at(i).v4(), Whadd1.v4());
      double dR_Wd2 = deltaR(hyp.tophad_jets().at(i).v4(), Whadd2.v4());
      double dR_bhad = deltaR(hyp.tophad_jets().at(i).v4(), bhad.v4());
      if(dR_Wd1 < 0.4){Wd1_merged = true;}
      if(dR_Wd2 < 0.4){Wd2_merged = true;}
      if(dR_bhad < 0.4){bhad_merged = true;}
 
      dR_hadtop += min({dR_Wd1, dR_Wd2, dR_bhad});
    }
    if(Wd1_merged && Wd2_merged && bhad_merged) tophad_match = true;

  }
  else{ //hadronic top reconstructed as AK8 jet
    double dR_Wd1_min = 1e6;
    double dR_Wd2_min = 1e6;
    double dR_bhad_min = 1e6;
    dR_Wd1_min = deltaR(hyp.tophad_topjet_ptr()->v4(), Whadd1.v4());
    dR_Wd2_min = deltaR(hyp.tophad_topjet_ptr()->v4(), Whadd2.v4());
    dR_bhad_min = deltaR(hyp.tophad_topjet_ptr()->v4(), bhad.v4());
    if(dR_Wd1_min<0.8 && dR_Wd2_min<0.8 && dR_bhad_min<0.8) tophad_match = true;
    dR_hadtop = dR_Wd1_min + dR_Wd2_min + dR_bhad_min;
  }

  //Leptonic top
  bool toplep_match = false;
  double dR_lep = deltaR(lepton.v4(),hyp.lepton().v4());
  double dR_neutrino = deltaR(neutrino.v4(),hyp.neutrino_v4());
  double dPhi_neutrino = deltaPhi(neutrino.v4(),hyp.neutrino_v4());  // for some reason we use dphi here, not dR. Why?
  double dR_blep_min = 1e6;
  double dR_other_toplepjets = 0;
  for (uint i = 0; i < hyp.toplep_jets().size(); i++){ // look into # of used jets here.
    double dR_blep = deltaR(hyp.toplep_jets().at(i).v4(), blep.v4());
    if(dR_blep_min>dR_blep){
      if(dR_blep_min < 1e6){dR_other_toplepjets += dR_blep_min;}
      dR_blep_min=dR_blep;
    }
    else{dR_other_toplepjets += dR_blep;}
  }
  if(dR_blep_min<0.4 && dR_lep<0.1 && dPhi_neutrino<0.3) toplep_match = true;
  
  double deltaM_lep = fabs(hyp.toplep_v4().M()-tlep.v4().M())/tlep.v4().M();
  double deltaM_had = fabs(hyp.tophad_v4().M()-thad.v4().M())/thad.v4().M();
  // cout<<"GEN: tlep.v4().M() = "<<tlep.v4().M()<<" thad.v4().M() = "<<thad.v4().M()<<endl;
  // cout<<"RECO: hyp.toplep_v4().M() = "<<hyp.toplep_v4().M()<<" hyp.tophad_v4().M() = "<<hyp.tophad_v4().M()<<endl;
  // cout<<"DELTA: "<<deltaM_lep<<" "<<deltaM_had<<endl;
  // cout<<" toplep_match, tophad_match :"<<toplep_match<<", "<<tophad_match<<endl;
  // cout<<" dR_neutrino = "<<dR_neutrino<<" dPhi_neutrino = "<<dPhi_neutrino<<" dRlep = "<<dR_lep<<" dR_blep_min = "<<dR_blep_min<<endl;
  // cout<<" dR_Wd1 = "<<dR_Wd1_min<<" dR_Wd2 = "<<dR_Wd2_min<<" dR_bhad = "<<dR_bhad_min<<endl;
  // vector<double> dR;//dRWd1had,dRWd2had,dRbhad,dRlep,dPhineutrino,dRblep
  // dR.push_back(dR_Wd1_min);
  // dR.push_back(dR_Wd2_min);
  // dR.push_back(dR_bhad_min);
  // dR.push_back(dR_lep);
  // dR.push_back(dPhi_neutrino);
  // dR.push_back(dRblep);

  double dR_sum = 0.;
  // hadronic top
  dR_sum+=dR_hadtop;   // other jets
  // leptonic top
  dR_sum+=dR_lep;                // lepton
  dR_sum+=dPhi_neutrino;         // neutrino
  dR_sum+=dR_blep_min;           // b quark
  dR_sum+=dR_other_toplepjets;   // other jets

  pair<bool,double> result;
  result.second = dR_sum;
  if(!toplep_match || !tophad_match){
    result.first = false;
  }
  else{
    result.first = true;
  }
  //  cout<<"### WE FOUND MATCH! ###"<<endl;
  //  return true;
  return result;
  // double dR_top_lep_reco_gen = deltaR( tlep.v4(), hyp.toplep_v4());
  // double dR_top_had_reco_gen = deltaR( thad.v4(), hyp.tophad_v4());
  // //  double dR_lep_reco_gen = deltaR(lepton.v4(),hyp.lepton().v4());
  // //  cout<<"Hi from TTbarSemiLepMatchableSelection::check_reco!"<<endl;
  // //  cout<<"dR_top_lep_reco_gen ="<<dR_top_lep_reco_gen<<" dR_top_had_reco_gen = "<<dR_top_had_reco_gen<<endl;
  // //  cout<<" dR_lep_reco_gen = "<<dR_lep_reco_gen<<endl;
  // if(dR_top_lep_reco_gen>0.4 || dR_top_had_reco_gen>0.4) return false;
 
  
  // if(deltaM_lep>1 || deltaM_had>1) return false;
  // return true;
}

////////////////////////////////////////////////////////


uhh2::METCut::METCut(float min_met, float max_met):
  min_met_(min_met), max_met_(max_met) {}

bool uhh2::METCut::passes(const uhh2::Event& event){

  assert(event.met);

  float MET = event.met->pt();
  return (MET > min_met_) && (MET < max_met_);
}

uhh2::STCut::STCut(float min_st, float max_st):
  min_st_(min_st), max_st_(max_st) {}

bool uhh2::STCut::passes(const uhh2::Event& event){

  assert(event.jets);

  float st_jets = 0;
  std::vector<Jet>* jets = event.jets;
  for(const auto & jet : *jets) st_jets += jet.pt();
  return (st_jets > min_st_) && (st_jets < max_st_);
}


uhh2::TopTagEventSelection::TopTagEventSelection(const TopJetId& tjetID):
  topjetID_(tjetID) {

}

bool uhh2::TopTagEventSelection::passes(const uhh2::Event& event){ 

  for(auto & topjet : * event.topjets){
    if(topjetID_(topjet, event)) return true; 
  }
  return false;
}
////////////////////////////////////////////////////////
