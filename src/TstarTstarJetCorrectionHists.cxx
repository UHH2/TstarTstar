#include "UHH2/TstarTstar/include/TstarTstarJetCorrectionHists.h"
#include "UHH2/core/include/Event.h"
#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

float inv_mass_3(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

using namespace std;
using namespace uhh2;

TstarTstarJetCorrectionHists::TstarTstarJetCorrectionHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

    book<TH1D>("ptrec_o_ptgen__AK4_leading", "leading AK4 p^{REC}_{T} / p^{GEN}_{T}", 100, 0, 2);
    book<TH1D>("ptrec_o_ptgen__AK4_second", "second AK4 p^{REC}_{T} / p^{GEN}_{T}", 100, 0, 2);
    book<TH1D>("ptrec_o_ptgen__AK4_all", "all AK4 p^{REC}_{T} / p^{GEN}_{T}", 100, 0, 2);

    book<TH1D>("ptrec_o_ptgen__HOTVR_leading", "leading HOTVR p^{REC}_{T} / p^{GEN}_{T}", 100, 0, 2);
    book<TH1D>("ptrec_o_ptgen__HOTVR_second", "second HOTVR p^{REC}_{T} / p^{GEN}_{T}", 100, 0, 2);
    book<TH1D>("ptrec_o_ptgen__HOTVR_all", "all HOTVR p^{REC}_{T} / p^{GEN}_{T}", 100, 0, 2);

    const int n_HOTVRpt_bins = 14;
    double HOTVRpt_bins[n_HOTVRpt_bins] = {0, 200, 225, 250, 275, 300, 350, 400, 500, 600, 800, 1000, 1500, 2000};
    book<TH2D>("2D_jetresponse_genjetpt_HOTVR", ";HOTVR jet response; HOTVR GEN pt", 100, 0, 2, n_HOTVRpt_bins-1, HOTVRpt_bins);

    const int n_AK4pt_bins = 14;
    double AK4pt_bins[n_AK4pt_bins] = {0, 30, 60, 90, 120, 150, 200, 250, 300, 400, 500, 750, 1000, 1500};
    book<TH2D>("2D_jetresponse_genjetpt_AK4", ";HOTVR jet response; HOTVR GEN pt", 100, 0, 2, n_AK4pt_bins-1, AK4pt_bins);
  
    book<TH1D>("lowresponse_jet_pt", "p_{T} [GeV]",  100, 0, 1000);
    book<TH1D>("lowresponse_jet_mass", "mass [GeV]",  100, 0, 1000);
    book<TH1D>("lowresponse_jet_eta", "#eta^{jet}", 50, -5.2, 5.2);
    book<TH1D>("lowresponse_jet_phi", "#phi^{jet}", 50, -3.2, 3.2);

    book<TH1D>("lowresponse_jet_ele_pt", "p_{T}^{ele} [GeV/c]", 50, 0, 1000);
    book<TH1D>("lowresponse_jet_ele_eta", "#eta^{ele}", 50, -5.2, 5.2);
    book<TH1D>("lowresponse_jet_ele_phi", "#phi^{ele}", 50, -3.2, 3.2);

    book<TH1D>("lowresponse_jet_mu_pt", "p_{T}^{#mu} [GeV/c]", 50, 0, 1000);
    book<TH1D>("lowresponse_jet_mu_eta", "#eta^{#mu}", 50, -5.2, 5.2);
    book<TH1D>("lowresponse_jet_mu_phi", "#phi^{#mu}", 50, -3.2, 3.2);

    book<TH1D>("lowresponse_dR_jet_lepton", "#DeltaR (jet, lepton)", 30, 0, 6);
    book<TH1D>("lowresponse_dR_jet_jet", "#DeltaR (jet, closest jet)", 30, 0, 6);

    book<TH1D>("highresponse_jet_pt", "p_{T} [GeV]",  100, 0, 1000);
    book<TH1D>("highresponse_jet_mass", "mass [GeV]",  100, 0, 1000);
    book<TH1D>("highresponse_jet_eta", "#eta^{jet}", 50, -5.2, 5.2);
    book<TH1D>("highresponse_jet_phi", "#phi^{jet}", 50, -3.2, 3.2);

    book<TH1D>("highresponse_jet_ele_pt", "p_{T}^{ele} [GeV/c]", 50, 0, 1000);
    book<TH1D>("highresponse_jet_ele_eta", "#eta^{ele}", 50, -5.2, 5.2);
    book<TH1D>("highresponse_jet_ele_phi", "#phi^{ele}", 50, -3.2, 3.2);

    book<TH1D>("highresponse_jet_mu_pt", "p_{T}^{#mu} [GeV/c]", 50, 0, 1000);
    book<TH1D>("highresponse_jet_mu_eta", "#eta^{#mu}", 50, -5.2, 5.2);
    book<TH1D>("highresponse_jet_mu_phi", "#phi^{#mu}", 50, -3.2, 3.2);

    book<TH1D>("highresponse_dR_jet_lepton", "#DeltaR (jet, lepton)", 30, 0, 6);
    book<TH1D>("highresponse_dR_jet_jet", "#DeltaR (jet, closest jet)", 30, 0, 6);

    matching_radius = 0.3;

}


void TstarTstarJetCorrectionHists::fill(const Event & event){

    double weight = event.weight;

    // AK4 part
    int AK4index = 0; // this will be used to determine which is the first and second

    std::vector<int> usedJets;                      // indices of all reco jets used for matching already

    for(const auto & AK4genjet : *event.genjets){ // loop over GENjets

        double mindR = 99999;                           // dR of the best match (starting at something huge)
        Jet matched_RECOjet;                            // the "current" best matched reco jet
        int current_matched_jet = -1;                   // index of the "currently" best matched jet
        int i = 0;                                      // loop index
        for(const auto & AK4recojet : *event.jets){

            double this_dR = deltaR(AK4recojet, AK4genjet);

            // only considering jets not used yet, and being closer than the previous best:
            if(this_dR < mindR && (std::find(usedJets.begin(), usedJets.end(), i)==usedJets.end())) {
                mindR = this_dR;
                matched_RECOjet = AK4recojet;
                current_matched_jet = i;
            }

            i++;

        }

        if (mindR < matching_radius) {

            if (current_matched_jet >= 0) usedJets.push_back(current_matched_jet); // log which jet is used

            // plotting
            if(AK4index == 0) hist("ptrec_o_ptgen__AK4_leading")->Fill( matched_RECOjet.pt() / AK4genjet.pt(), weight);
            if(AK4index == 1) hist("ptrec_o_ptgen__AK4_second")->Fill( matched_RECOjet.pt() / AK4genjet.pt(), weight);
            hist("ptrec_o_ptgen__AK4_all")->Fill( matched_RECOjet.pt() / AK4genjet.pt(), weight);
            ((TH2D*)hist("2D_jetresponse_genjetpt_AK4"))->Fill( matched_RECOjet.pt() / AK4genjet.pt(), AK4genjet.pt(), weight);
        
            if(matched_RECOjet.pt() / AK4genjet.pt() < 0.25) { // making some plots for low response jets

                hist("lowresponse_jet_pt")->Fill(matched_RECOjet.pt(), weight);
                hist("lowresponse_jet_mass")->Fill(inv_mass_3(matched_RECOjet.v4()), weight);
                hist("lowresponse_jet_eta")->Fill(matched_RECOjet.eta(), weight);
                hist("lowresponse_jet_phi")->Fill(matched_RECOjet.phi(), weight);

                // assuming this will be exactly one!
                Particle lepton;

                for (const Muon & thismu : *event.muons){
                    hist("lowresponse_jet_mu_pt")->Fill(thismu.pt(), weight);
                    hist("lowresponse_jet_mu_eta")->Fill(thismu.eta(), weight);
                    hist("lowresponse_jet_mu_phi")->Fill(thismu.phi(), weight);
                    lepton = thismu;
                }

                for (const Electron & thisele : *event.electrons){
                    hist("lowresponse_jet_ele_pt")->Fill(thisele.pt(), weight);
                    hist("lowresponse_jet_ele_eta")->Fill(thisele.eta(), weight);
                    hist("lowresponse_jet_ele_phi")->Fill(thisele.phi(), weight);
                    lepton = thisele;
                }

                hist("lowresponse_dR_jet_lepton")->Fill(deltaR(lepton, matched_RECOjet), weight);

                int j = 0;
                double mindR = 99999;
                Jet matched_matched_jet;
                for(const auto & AK4recojet : *event.jets){
                    double this_dR = deltaR(matched_RECOjet, AK4recojet);
                    if(this_dR < mindR && j != current_matched_jet) {
                        mindR = this_dR;
                        matched_matched_jet = AK4recojet;
                    }
                }
                hist("lowresponse_dR_jet_jet")->Fill(deltaR(matched_RECOjet, matched_matched_jet), weight);


            } else {

                hist("highresponse_jet_pt")->Fill(matched_RECOjet.pt(), weight);
                hist("highresponse_jet_mass")->Fill(inv_mass_3(matched_RECOjet.v4()), weight);
                hist("highresponse_jet_eta")->Fill(matched_RECOjet.eta(), weight);
                hist("highresponse_jet_phi")->Fill(matched_RECOjet.phi(), weight);

                // assuming this will be exactly one!
                Particle lepton;

                for (const Muon & thismu : *event.muons){
                    hist("highresponse_jet_mu_pt")->Fill(thismu.pt(), weight);
                    hist("highresponse_jet_mu_eta")->Fill(thismu.eta(), weight);
                    hist("highresponse_jet_mu_phi")->Fill(thismu.phi(), weight);
                    lepton = thismu;
                }

                for (const Electron & thisele : *event.electrons){
                    hist("highresponse_jet_ele_pt")->Fill(thisele.pt(), weight);
                    hist("highresponse_jet_ele_eta")->Fill(thisele.eta(), weight);
                    hist("highresponse_jet_ele_phi")->Fill(thisele.phi(), weight);
                    lepton = thisele;
                }

                hist("highresponse_dR_jet_lepton")->Fill(deltaR(lepton, matched_RECOjet), weight);

                int j = 0;
                double mindR = 99999;
                Jet matched_matched_jet;
                for(const auto & AK4recojet : *event.jets){
                    double this_dR = deltaR(matched_RECOjet, AK4recojet);
                    if(this_dR < mindR && j != current_matched_jet) {
                        mindR = this_dR;
                        matched_matched_jet = AK4recojet;
                    }
                }
                hist("highresponse_dR_jet_jet")->Fill(deltaR(matched_RECOjet, matched_matched_jet), weight);

            }

        }

        AK4index++;
    }

    // HOTVR part
    int HOTVRindex = 0; // this will be used to determine which is the first and second

    for(const auto & HOTVRgenjet : *event.gentopjets){

        // we have the jet, so we need to find the pt of the matching GEN jet

        double mindR = 99999;
        TopJet matched_RECOjet;
        int i = 0;
        int current_matched_jet = -1;
        std::vector<int> usedJets;
        for(const auto & HOTVRrecojet : *event.topjets){

            double this_dR = deltaR(HOTVRrecojet, HOTVRgenjet);
            if(this_dR < mindR && (std::find(usedJets.begin(), usedJets.end(), i)==usedJets.end())) {
                mindR = this_dR;
                matched_RECOjet = HOTVRrecojet;
                current_matched_jet = i;
            }

            if (current_matched_jet >= 0) usedJets.push_back(current_matched_jet);
            i++;

        }

        if (mindR < matching_radius) {
            if(HOTVRindex == 0) hist("ptrec_o_ptgen__HOTVR_leading")->Fill( matched_RECOjet.pt() / HOTVRgenjet.pt(), weight);
            if(HOTVRindex == 1) hist("ptrec_o_ptgen__HOTVR_second")->Fill( matched_RECOjet.pt() / HOTVRgenjet.pt(), weight);
            hist("ptrec_o_ptgen__HOTVR_all")->Fill( matched_RECOjet.pt() / HOTVRgenjet.pt(), weight);
            ((TH2D*)hist("2D_jetresponse_genjetpt_HOTVR"))->Fill( matched_RECOjet.pt() / HOTVRgenjet.pt(), HOTVRgenjet.pt(), weight);

        }
        HOTVRindex++;
    }

}


TstarTstarJetCorrectionHists::~TstarTstarJetCorrectionHists(){}
