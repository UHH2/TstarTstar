#include "UHH2/TstarTstar/include/TstarTstarSFHists.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"


#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

TstarTstarSFHists::TstarTstarSFHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  h_weight_sfmu_id = ctx.get_handle<float>("weight_sfmu_id");
  h_weight_sfmu_id_down = ctx.get_handle<float>("weight_sfmu_id_down");
  h_weight_sfmu_id_up = ctx.get_handle<float>("weight_sfmu_id_up");

  h_weight_sfmu_isolation = ctx.get_handle<float>("weight_sfmu_isolation");
  h_weight_sfmu_isolation_down = ctx.get_handle<float>("weight_sfmu_isolation_down");
  h_weight_sfmu_isolation_up = ctx.get_handle<float>("weight_sfmu_isolation_up");

  h_weight_sfele_id = ctx.get_handle<float>("weight_sfelec_id");
  h_weight_sfele_id_down = ctx.get_handle<float>("weight_sfelec_id_down");
  h_weight_sfele_id_up = ctx.get_handle<float>("weight_sfelec_id_up");

  // book all histograms here
  book<TH1F>("weight", "weight", 55, -1, 10);
  book<TH1F>("weight_2", "weight_2", 55, -10, 100);
  book<TH1F>("weight_3", "weight_3", 55, -100, 1000);
  book<TH1F>("weight_4", "weight_4", 55, -1000, 10000);
  book<TH1F>("weight_5", "weight_5", 55, -10000, 100000);

  // a lot of SF hists to check if everything is working as it should
  book<TH1F>("ID_SF_muon", "SF muon ID", 50, 0, 2);
  book<TH1F>("ID_SF_muon_up", "SF muon ID up", 50, 0, 2);
  book<TH1F>("ID_SF_muon_down", "SF muon ID down", 50, 0, 2);

  book<TH1F>("ISO_SF_muon", "SF muon ISO", 50, 0, 2);
  book<TH1F>("ISO_SF_muon_up", "SF muon ISO up", 50, 0, 2);
  book<TH1F>("ISO_SF_muon_down", "SF muon ISO down", 50, 0, 2);

  book<TH1F>("ID_SF_ele", "SF ele ID", 50, 0, 2);
  book<TH1F>("ID_SF_ele_up", "SF ele ID up", 50, 0, 2);
  book<TH1F>("ID_SF_ele_down", "SF ele ID down", 50, 0, 2);

}


void TstarTstarSFHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  bool debug = false;

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  hist("weight")->Fill(weight);
  hist("weight_2")->Fill(weight);
  hist("weight_3")->Fill(weight);
  hist("weight_4")->Fill(weight);
  hist("weight_5")->Fill(weight);

  hist("ID_SF_muon")->Fill(event.get(h_weight_sfmu_id), weight);
  hist("ID_SF_muon_up")->Fill(event.get(h_weight_sfmu_id_up), weight);
  hist("ID_SF_muon_down")->Fill(event.get(h_weight_sfmu_id_down), weight);

  hist("ISO_SF_muon")->Fill(event.get(h_weight_sfmu_isolation), weight);
  hist("ISO_SF_muon_up")->Fill(event.get(h_weight_sfmu_isolation_up), weight);
  hist("ISO_SF_muon_down")->Fill(event.get(h_weight_sfmu_isolation_down), weight);

  hist("ID_SF_ele")->Fill(event.get(h_weight_sfele_id), weight);
  hist("ID_SF_ele_up")->Fill(event.get(h_weight_sfele_id_up), weight);
  hist("ID_SF_ele_down")->Fill(event.get(h_weight_sfele_id_down), weight);

  if(debug) cout << "Finished SF Hists!" << endl;

  }

TstarTstarSFHists::~TstarTstarSFHists(){}
