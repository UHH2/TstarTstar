#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/TstarTstar/include/TstarTstarScaleFactors.h"

#include <stdexcept>
#include <set>
#include <regex>

using namespace std;
using namespace uhh2;

namespace {
  // invariant mass of a lorentzVector, but safe for timelike / spacelike vectors
  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }
}

// Generic Class for Applying SFs
void ScaleFactorsFromHistos::LoadHisto(TFile* file, std::string name, std::string hname) {
  histos[name].reset((TH1F*)file->Get(hname.c_str()));
  histos[name]->SetDirectory(0);
};

double ScaleFactorsFromHistos::Evaluator(std::string hname, double var) {
  // invalid cases
  if (var == uhh2::infinity) return 1.0;

  int firstBin = 1;
  int lastBin  = histos[hname]->GetNbinsX();
  double h_min = histos[hname]->GetBinCenter(firstBin)-0.5*histos[hname]->GetBinWidth(firstBin);
  double h_max = histos[hname]->GetBinCenter(lastBin)+0.5*histos[hname]->GetBinWidth(lastBin);
  double var_for_eval = var;
  var_for_eval = (var_for_eval > h_min) ? var_for_eval : h_min+0.001;
  var_for_eval = (var_for_eval < h_max) ? var_for_eval : h_max-0.001;
  return histos[hname]->GetBinContent(histos[hname]->FindBin(var_for_eval));
};

// Apply Theory weights
NLOCorrections::NLOCorrections(uhh2::Context& ctx) {

  std::string dataset_version = ctx.get("dataset_version");
  // Corrections for 2017 and 2018 are the same. 2016 is different
  is2016 = FindInString("UL16", dataset_version);

  //TODO it's arbitrary.
  is_Wjets  = FindInString("WJets",dataset_version);
  is_Znn    = FindInString("DY_inv",dataset_version);
  is_DY     = FindInString("DY",dataset_version) && !is_Znn;
  is_Zjets  = is_DY || is_Znn;

  std::string folder_ = ctx.get("NLOCorrections")+"/"; //TODO better name
  for (const std::string& proc: {"w","z"}) {
    TFile* file_ = new TFile((folder_+"merged_kfactors_"+proc+"jets.root").c_str());
    for (const std::string& corr: {"ewk","qcd","qcd_ewk"}) LoadHisto(file_, proc+"_"+corr, "kfactor_monojet_"+corr);
    file_->Close();
  }
  for (const std::string& proc: {"dy","znn"}) {
    TFile* file_ = new TFile((folder_+"kfac_"+proc+"_filter.root").c_str());
    LoadHisto(file_, proc+"_qcd_2017", "kfac_"+proc+"_filter");
    file_->Close();
  }
  TFile* file_ = new TFile((folder_+"2017_gen_v_pt_qcd_sf.root").c_str());
  LoadHisto(file_, "w_qcd_2017", "wjet_dress_inclusive");
  file_->Close();
  file_ = new TFile((folder_+"lindert_qcd_nnlo_sf.root").c_str());
  for (const std::string& proc: {"eej", "evj", "vvj"}) LoadHisto(file_, proc+"_qcd_nnlo", proc);
  file_->Close();

}

double NLOCorrections::GetPartonObjectPt(uhh2::Event& event, ParticleID objID) {
  for(const auto & gp : *event.genparticles) {if (gp.pdgId()==objID) return gp.pt(); }
  return uhh2::infinity;
};


bool NLOCorrections::process(uhh2::Event& event){
  // Sample dependant corrections
  if ((!is_Wjets && !is_Zjets) || event.isRealData) return true;
  // if ((is_Wjets || is_Zjets) && !is2016) return true;//check
  double objpt = uhh2::infinity, theory_weight = 1.0;
  std::string process = "";

  if (is_Zjets) objpt = GetPartonObjectPt(event,ParticleID::Z);
  if (is_Wjets) objpt = GetPartonObjectPt(event,ParticleID::W);

  if (is_Zjets) process = "z";
  if (is_Wjets) process = "w";

  if (do_QCD_EWK) theory_weight *= Evaluator(process+"_qcd_ewk",objpt);
  else {
    if (do_EWK) theory_weight *= Evaluator(process+"_ewk",objpt);
    if (do_QCD_NLO) {
      if (!is2016) {
        if (is_DY)  process = "dy";
        if (is_Znn) process = "znn";
      }
      theory_weight *= Evaluator(process+"_qcd"+(is2016?"":"_2017"),objpt);
    }
  }

  if (do_QCD_NNLO) {
    if (is_DY)    process = "eej";
    if (is_Znn)   process = "vvj";
    if (is_Wjets) process = "evj";
    theory_weight *= Evaluator(process+"_qcd_nnlo",objpt);
  }

  event.weight *= theory_weight;
  return true;
}

EtaPhiEventCleanerMC::EtaPhiEventCleanerMC(Context& ctx, float weight_factor_, float min_eta_, float max_eta_, float min_phi_, float max_phi_,
  string jetCollection_, bool doJets_, bool doElectrons_, bool doMuons_):
  weight_factor(weight_factor_),
  min_eta(min_eta_), max_eta(max_eta_),
  min_phi(min_phi_), max_phi(max_phi_), jetCollection(jetCollection_),
  doJets(doJets_), doElectrons(doElectrons_), doMuons(doMuons_){
    h_jets = ctx.get_handle<vector<Jet>>(jetCollection);
    h_topjets = ctx.get_handle<vector<TopJet>>(jetCollection);
  }

/* Remove events in data if at least one electron/muon/jet (depending on doJets,
doElectrons and doMuon) is found in area specified by the eta and phi
intervals. Runs only if the run rumber is in the interval [min_run, max_run]
(including bounds). Can pass -1 as max_run to ignore max_run. */
bool EtaPhiEventCleanerMC::process(Event& event){
  // only run on MC
  if (event.isRealData) {
    return true;
  }

  if (doJets){
    // Get the jets for the event, depening on whether they are jets or topJets.
    vector<Jet> jets;
    if (event.is_valid(h_topjets)) jets.assign((&event.get(h_topjets))->begin(), (&event.get(h_topjets))->end());
    else if (event.is_valid(h_jets)) jets.assign((&event.get(h_jets))->begin(), (&event.get(h_jets))->end());
    else throw logic_error("DetectorCleaning: No valid jet collection given.");

    for (auto & jet: jets){
      if (min_eta < jet.eta() && jet.eta() < max_eta){
        if (min_phi < jet.phi() && jet.phi() < max_phi){
          event.weight *= weight_factor;
          return true;
        }
      }
    }
  }

  if (doElectrons){
    for(const auto & ele : *event.electrons){
      if (min_eta < ele.eta() && ele.eta() < max_eta){
        if (min_phi < ele.phi() && ele.phi() < max_phi){
          event.weight *= weight_factor;
          return true;
        }
      }
    }
  }

  if (doMuons){
    for(const auto & muon : *event.muons){
      if (min_eta < muon.eta() && muon.eta() < max_eta){
        if (min_phi < muon.phi() && muon.phi() < max_phi){
          event.weight *= weight_factor;
          return true;
        }
      }
    }
  }
  return true;
}

HEMCleanerMCScale::HEMCleanerMCScale(Context& ctx, string jetCollection, bool doJets,
  bool doElectrons, bool doMuons):
  EtaPhiEventCleanerMC(ctx, weight_factor, min_eta_HEM,
    max_eta_HEM,min_phi_HEM, max_phi_HEM, jetCollection, doJets, doElectrons,
    doMuons){};


TstarTstarSpinScale::TstarTstarSpinScale(Context& ctx, TString path) {

  if( (ctx.get("SignalSpinScenario", "not set") == "_32") && (ctx.get("dataset_version").find("Tstar") != std::string::npos)) {

    doSpinReweighting = true;

    // detect mass string
    std::string masspoint = std::regex_replace(ctx.get("dataset_version"), std::regex("[^0-9]*([0-9]+).*"), std::string("$1"));
    std::cout << "Detected this to be a signal sample at mass point " << masspoint << " for spin reweighting." << std::endl;
    
    // if we are signal, read the file!
    TFile *f = new TFile(path + "/spin12tospin32_M-" + masspoint + ".root");
    scaleFactorHist = (TH1D*)f->Get("spinweights");

  } else {
    std::cout << "Not doing any spin reweighting." << std::endl;
  }

}

bool TstarTstarSpinScale::process(Event& event){

  if(!doSpinReweighting) return true; // not doing anything

  // Find tstars
  GenParticle tstar, antitstar;
  bool found_tstar = false, found_antitstar = false;
  for(const GenParticle & gp : *event.genparticles){
    if(gp.pdgId() == 600 && (gp.status()==23 || gp.status()==22)){
      tstar = gp;
      if(found_tstar) std::cout << "Error: found two T*" << std::endl;
      found_tstar = true;
    }
    else if(gp.pdgId() == -600 && (gp.status()==23 || gp.status()==22)){
      antitstar = gp;
      if(found_antitstar) std::cout << "Error: found two anti T*" << std::endl;
      found_antitstar = true;
    }
  }

  if(!found_tstar || !found_antitstar) return false;

  float m_tstartstar = inv_mass(tstar.v4()+antitstar.v4());

  double sf = scaleFactorHist->GetBinContent( scaleFactorHist->GetXaxis()->FindBin( m_tstartstar ));

  event.weight *= sf;

  return true;
}
