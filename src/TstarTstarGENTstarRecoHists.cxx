#include "UHH2/TstarTstar/include/TstarTstarGENTstarRecoHists.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/TstarTstar/include/ReconstructionTstarHypothesis.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

float inv_mass_77(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }


TstarTstarGENTstarRecoHists::TstarTstarGENTstarRecoHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  h_tstartstar_hyp_gHOTVR = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_GENHyp_gHOTVR");
  h_tstartstar_hyp_gAK4 = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_GENHyp_gAK4");

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

}

void TstarTstarGENTstarRecoHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  bool debug = false;

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  // reco plots
  //gHOTVR
  try {
    ReconstructionTstarHypothesis hyp_gHOTVR = event.get(h_tstartstar_hyp_gHOTVR);
    hist("M_Tstar_gHOTVR_uncorr")->Fill(inv_mass_77(hyp_gHOTVR.tstarlep_v4()), weight/2);
    if(hyp_gHOTVR.tstarhad_v4().pt() != 0) {
      ReconstructionHypothesis ttbarhyp_gHOTVR = hyp_gHOTVR.ttbar_hyp();
      hist("chi2_gHOTVR")->Fill(hyp_gHOTVR.chi2(), weight);
      hist("M_top_gHOTVR_had")->Fill(inv_mass_77(ttbarhyp_gHOTVR.tophad_v4()), weight);
      hist("M_top_gHOTVR_lep")->Fill(inv_mass_77(ttbarhyp_gHOTVR.toplep_v4()), weight);
      hist("M_Tstar_gHOTVR")->Fill(inv_mass_77(hyp_gHOTVR.tstarlep_v4()), weight/2);
      hist("M_Tstar_gHOTVR_lep")->Fill(inv_mass_77(hyp_gHOTVR.tstarlep_v4()), weight);
      hist("M_Tstar_gHOTVR")->Fill(inv_mass_77(hyp_gHOTVR.tstarhad_v4()), weight/2);
      hist("M_Tstar_gHOTVR_had")->Fill(inv_mass_77(hyp_gHOTVR.tstarhad_v4()), weight);
    }
  } catch(...) {}

  //gAK4
  try {
    ReconstructionTstarHypothesis hyp_gAK4 = event.get(h_tstartstar_hyp_gAK4);
    hist("M_Tstar_gAK4_uncorr")->Fill(inv_mass_77(hyp_gAK4.tstarlep_v4()), weight/2);
    if(hyp_gAK4.tstarhad_v4().pt() != 0) {
      ReconstructionHypothesis ttbarhyp_gAK4 = hyp_gAK4.ttbar_hyp();
      hist("chi2_gAK4")->Fill(hyp_gAK4.chi2(), weight);
      hist("M_top_gAK4_had")->Fill(inv_mass_77(ttbarhyp_gAK4.tophad_v4()), weight);
      hist("M_top_gAK4_lep")->Fill(inv_mass_77(ttbarhyp_gAK4.toplep_v4()), weight);
      hist("M_Tstar_gAK4")->Fill(inv_mass_77(hyp_gAK4.tstarlep_v4()), weight/2);
      hist("M_Tstar_gAK4_lep")->Fill(inv_mass_77(hyp_gAK4.tstarlep_v4()), weight);
      hist("M_Tstar_gAK4")->Fill(inv_mass_77(hyp_gAK4.tstarhad_v4()), weight/2);
      hist("M_Tstar_gAK4_had")->Fill(inv_mass_77(hyp_gAK4.tstarhad_v4()), weight);
    }
  } catch(...) {}

  if(debug) cout << "Finished Tstar Hists!" << endl;

}

TstarTstarGENTstarRecoHists::~TstarTstarGENTstarRecoHists(){}
