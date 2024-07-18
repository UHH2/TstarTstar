#include <UHH2/TstarTstar/include/PuppiCHSMatching.h>
#include <UHH2/core/include/LorentzVector.h>
#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>
#include "UHH2/core/include/Event.h"

using namespace std;
using namespace uhh2;

PuppiCHSMatching::PuppiCHSMatching(uhh2::Context& ctx, std::string chs_jet_handle = "jetsAk4CHS"){

    // this usage assumes that PUPPI jets have been set as the "main" jets of the analysis
    // thus we're just getting the CHS jets here as additional collections
    // in the end, CHS_matched is the resulting one that needs to be used
    h_CHSjets = ctx.get_handle< std::vector<Jet> >(chs_jet_handle);
    h_CHS_matched_ = ctx.declare_event_output<vector<Jet>>("CHS_matched");

  }

  bool PuppiCHSMatching::process(uhh2::Event& event){

    vector<Jet> CHSjets = event.get(h_CHSjets);
    std::vector<Jet> matched_jets;
    std::vector<Jet> matched_jets_PUPPI;

    for(const Jet & jet : *event.jets){ // PUPPI jets
      double deltaR_min = 99;

      for(const Jet & CHSjet : CHSjets){ // CHS jets
        double deltaR_CHS = deltaR(jet,CHSjet);
        if(deltaR_CHS<deltaR_min) deltaR_min = deltaR_CHS;
      } // end CHS loop

      if(deltaR_min>0.2) continue;

      for(const Jet & CHSjet : CHSjets){
        if(deltaR(jet,CHSjet)!=deltaR_min) continue;
        else{
          matched_jets.emplace_back(CHSjet);
          matched_jets_PUPPI.emplace_back(jet);
        }
      }

    } // end PUPPI loop
    std::swap(matched_jets_PUPPI, *event.jets);
    event.set(h_CHS_matched_, matched_jets);
    if(event.jets->size()==0) return false;
    return true;
  }