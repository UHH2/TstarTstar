#include "UHH2/TstarTstar/include/TstarTstarSignalRegionHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include "TFile.h"

using namespace std;
using namespace uhh2;

TstarTstarSignalRegionHists::TstarTstarSignalRegionHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  // ST handle
  h_ST_HOTVR = ctx.get_handle<double>("ST_HOTVR");
  is_MC = ctx.get("dataset_type") == "MC";

  if(ctx.get("debug", "<not set>") == "true") debug = true;

  std::vector<TString> samples_that_need_other_mcweight_handling = {"TstarTstar", "DYJets", "WJets", "QCD_HT"};
  needsOtherMCweightHandling = false;
  for (const auto sample : samples_that_need_other_mcweight_handling) {
    if (ctx.get("dataset_version").find(sample) != std::string::npos) needsOtherMCweightHandling = true;
  }
  if(needsOtherMCweightHandling && debug) std::cout << "We are signal so we need other idices" << std::endl;

  // Histogram containing decorrelation uncertainty
  TString sample_string = "total";      // won't be used, just to have it filled for mc only results :
  if(is_MC) { 
    if(ctx.get("dataset_version").find("TT") != std::string::npos) sample_string = "TTbar";
    else if(ctx.get("dataset_version").find("ST") != std::string::npos) sample_string = "ST";
  }
  TFile *decorrelationUncertaintyFile = new TFile("/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/decorrelationComparison_" + sample_string + "_smooth.root");
  decorrelationUncertainty = (TH1*)decorrelationUncertaintyFile->Get("decorrelation_uncertainty");

  // histograms for yield correction
  h_flag_muonevent = ctx.get_handle<bool>("is_muevt");
  TString year = ctx.get("year", "<not set>");
  if(year == "<not set>"){
    if(ctx.get("dataset_version").find("2016") != std::string::npos) year = "2016";
    else if(ctx.get("dataset_version").find("2017") != std::string::npos) year = "2017";
    else if(ctx.get("dataset_version").find("2018") != std::string::npos) year = "2018";
    else if(ctx.get("dataset_version").find("UL16preVFP") != std::string::npos) year = "UL16preVFP";
    else if(ctx.get("dataset_version").find("UL16postVFP") != std::string::npos) year = "UL16postVFP";
    else if(ctx.get("dataset_version").find("UL17") != std::string::npos) year = "UL17";
    else if(ctx.get("dataset_version").find("UL18") != std::string::npos) year = "UL18";
    else throw "No year found in dataset name!";
  }
  sample_string = ""; // we'll get all here, but only ST and TTbar are used
  if(ctx.get("dataset_version").find("TT") != std::string::npos) sample_string = "TTbar";
  else if(ctx.get("dataset_version").find("ST") != std::string::npos) sample_string = "ST";
  else if(ctx.get("dataset_version").find("WJets") != std::string::npos) sample_string = "WJets";
  else if(ctx.get("dataset_version").find("QCD") != std::string::npos) sample_string = "QCD";
  else if(ctx.get("dataset_version").find("Diboson") != std::string::npos) sample_string = "VV";
  else if(ctx.get("dataset_version").find("DY") != std::string::npos) sample_string = "DYJets";
  else if(ctx.get("dataset_version").find("TstarTstarToTgammaTgamma") != std::string::npos) sample_string = "TstarTstar";
  else if(ctx.get("dataset_version").find("TstarTstarToTgluonTgluon_Spin32") != std::string::npos) sample_string = "TstarTstar_Spin32";
  TFile *btagyield_eleFile = new TFile("/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/btagyield/btagYieldUncs_" + year + "_ele.root");
  btagyield_ele = (TH1*)btagyield_eleFile->Get(sample_string);
  TFile *btagyield_muFile = new TFile("/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/btagyield/btagYieldUncs_" + year + "_mu.root");
  btagyield_mu = (TH1*)btagyield_muFile->Get(sample_string);

  // weight handles

  // top pt reweighting
  h_weight_ttbar = ctx.get_handle<float>("weight_ttbar");

  // pileup reweighting
  h_weight_puNominal = ctx.get_handle<float>("weight_pu");
  h_weight_puUp = ctx.get_handle<float>("weight_pu_up");
  h_weight_puDown = ctx.get_handle<float>("weight_pu_down");

  // prefiring
  h_prefiringWeightNominal = ctx.get_handle<float>("prefiringWeight");
  h_prefiringWeightUp = ctx.get_handle<float>("prefiringWeightUp");
  h_prefiringWeightDown = ctx.get_handle<float>("prefiringWeightDown");

  // all the b-tagging stuff
  h_weight_btagdiscNominal = ctx.get_handle<float>("weight_btagdisc_central");

  h_weight_btagdisc_lfUp = ctx.get_handle<float>("weight_btagdisc_lf_up");
  h_weight_btagdisc_lfDown = ctx.get_handle<float>("weight_btagdisc_lf_down");
  h_weight_btagdisc_hfUp = ctx.get_handle<float>("weight_btagdisc_hf_up");
  h_weight_btagdisc_hfDown = ctx.get_handle<float>("weight_btagdisc_hf_down");
  h_weight_btagdisc_hfstats1Up = ctx.get_handle<float>("weight_btagdisc_hfstats1_up");
  h_weight_btagdisc_hfstats1Down = ctx.get_handle<float>("weight_btagdisc_hfstats1_down");
  h_weight_btagdisc_hfstats2Up = ctx.get_handle<float>("weight_btagdisc_hfstats2_up");
  h_weight_btagdisc_hfstats2Down = ctx.get_handle<float>("weight_btagdisc_hfstats2_down");
  h_weight_btagdisc_lfstats1Up = ctx.get_handle<float>("weight_btagdisc_lfstats1_up");
  h_weight_btagdisc_lfstats1Down = ctx.get_handle<float>("weight_btagdisc_lfstats1_down");
  h_weight_btagdisc_lfstats2Up = ctx.get_handle<float>("weight_btagdisc_lfstats2_up");
  h_weight_btagdisc_lfstats2Down = ctx.get_handle<float>("weight_btagdisc_lfstats2_down");
  h_weight_btagdisc_cferr1Up = ctx.get_handle<float>("weight_btagdisc_cferr1_up");
  h_weight_btagdisc_cferr1Down = ctx.get_handle<float>("weight_btagdisc_cferr1_down");
  h_weight_btagdisc_cferr2Up = ctx.get_handle<float>("weight_btagdisc_cferr2_up");
  h_weight_btagdisc_cferr2Down = ctx.get_handle<float>("weight_btagdisc_cferr2_down");

  // lepton stuff
  h_weight_sfelec_idNominal = ctx.get_handle<float>("weight_sfelec_id");
  h_weight_sfelec_idUp = ctx.get_handle<float>("weight_sfelec_id_up");
  h_weight_sfelec_idDown = ctx.get_handle<float>("weight_sfelec_id_down");

  h_weight_sfelec_triggerNominal = ctx.get_handle<float>("weight_sfelec_trigger");
  h_weight_sfelec_triggerUp = ctx.get_handle<float>("weight_sfelec_trigger_up");
  h_weight_sfelec_triggerDown = ctx.get_handle<float>("weight_sfelec_trigger_down");

  h_weight_sfelec_recoNominal = ctx.get_handle<float>("weight_sfelec_reco");
  h_weight_sfelec_recoUp = ctx.get_handle<float>("weight_sfelec_reco_up");
  h_weight_sfelec_recoDown = ctx.get_handle<float>("weight_sfelec_reco_down");

  h_weight_sfmu_idNominal = ctx.get_handle<float>("weight_sfmu_id");
  h_weight_sfmu_idUp = ctx.get_handle<float>("weight_sfmu_id_up");
  h_weight_sfmu_idDown = ctx.get_handle<float>("weight_sfmu_id_down");

  h_weight_sfmu_isoNominal = ctx.get_handle<float>("weight_sfmu_iso");
  h_weight_sfmu_isoUp = ctx.get_handle<float>("weight_sfmu_iso_up");
  h_weight_sfmu_isoDown = ctx.get_handle<float>("weight_sfmu_iso_down");

  h_weight_sfmu_triggerNominal = ctx.get_handle<float>("weight_sfmu_trigger");
  h_weight_sfmu_triggerUp = ctx.get_handle<float>("weight_sfmu_trigger_up");
  h_weight_sfmu_triggerDown = ctx.get_handle<float>("weight_sfmu_trigger_down");

  h_murmuf_upup        = ctx.get_handle< float >("weight_murmuf_upup");
  h_murmuf_upnone      = ctx.get_handle< float >("weight_murmuf_upnone");
  h_murmuf_noneup      = ctx.get_handle< float >("weight_murmuf_noneup");
  h_murmuf_nonedown    = ctx.get_handle< float >("weight_murmuf_nonedown");
  h_murmuf_downnone    = ctx.get_handle< float >("weight_murmuf_downnone");
  h_murmuf_downdown    = ctx.get_handle< float >("weight_murmuf_downdown");

  // binning definition
  const int nbins = 34;
  double bins[nbins] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3250, 4000, 6000};
  //const int nbins = 24;
  //double bins[nbins] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 4000, 6000};

  // nominal histogram
  book<TH1F>("pt_ST_nominal", "S_{T} [GeV]", nbins-1, bins);

  // list all variations here
  // for each of them, an up- and down variation hist is created!!!
  variations = {
    "pt_ST_topPt", "pt_ST_pu", "pt_ST_prefiring", "pt_ST_btagging_total", "pt_ST_btagging_hf", "pt_ST_btagging_hfstats1", "pt_ST_btagging_hfstats2", "pt_ST_btagging_lf",
    "pt_ST_btagging_lfstats1", "pt_ST_btagging_lfstats2", "pt_ST_btagging_cferr1", "pt_ST_btagging_cferr2", "pt_ST_sfelec_id", "pt_ST_sfelec_trigger", "pt_ST_sfelec_reco",
    "pt_ST_sfmu_id", "pt_ST_sfmu_iso", "pt_ST_sfmu_trigger", "pt_ST_decorrelation", "pt_ST_btagYield",
  };

  for (const auto variation : variations) {
    book<TH1F>(variation + "Up", "S_{T} [GeV]", nbins-1, bins);
    book<TH1F>(variation + "Down", "S_{T} [GeV]", nbins-1, bins);
  }

  // 100 histograms for the PDF stuff
  for(int i=0; i<100; i++){
    std::stringstream ss_name;
    ss_name << "pt_ST_PDF_" << i+1;

    stringstream ss_title;
    ss_title << "S_T [GeV] for PDF No. "  << i+1 << " out of 100" ;

    std::string s_name = ss_name.str();
    std::string s_title = ss_title.str();
    const char* char_name = s_name.c_str();
    const char* char_title = s_title.c_str();

    hist_names[i] = s_name;

    book<TH1F>(char_name, char_title, nbins-1, bins);

  }

  book<TH1F>("pt_ST_murmuf_upup", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_murmuf_upnone", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_murmuf_noneup", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_murmuf_nonedown", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_murmuf_downnone", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_murmuf_downdown", "S_{T} [GeV]", nbins-1, bins);
}


void TstarTstarSignalRegionHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  if(debug) cout << "Starting Tstar Hists." << endl;

  double st = event.get(h_ST_HOTVR);
  if(st > 6000) st = 5999.9; // handling overflow

  // fill nominal
  hist("pt_ST_nominal")->Fill(st, weight);

  // all others only make sense for MC. Fill them with "nominal" data here
  if(!is_MC) {

    for (const auto variation : variations) {
      hist(variation+"Up")->Fill(st, weight);
      hist(variation+"Down")->Fill(st, weight);
    }
    return;

  }

  // topPt
  hist("pt_ST_topPtDown")->Fill(st, weight);
  try {
    hist("pt_ST_topPtUp")->Fill(st, weight/event.get(h_weight_ttbar));
  } catch (...) { // catching it h_weight_ttbar was not filled...
    hist("pt_ST_topPtUp")->Fill(st, weight);
    if(debug) std::cout << "No top PT weight found, filling with nominal." << std::endl;
  }

  // pu
  hist("pt_ST_puUp")->Fill(st, event.get(h_weight_puUp)*weight/event.get(h_weight_puNominal));
  hist("pt_ST_puDown")->Fill(st, event.get(h_weight_puDown)*weight/event.get(h_weight_puNominal));

  // prefiring
  hist("pt_ST_prefiringUp")->Fill(st, event.get(h_prefiringWeightUp)*weight/event.get(h_prefiringWeightNominal));
  hist("pt_ST_prefiringDown")->Fill(st, event.get(h_prefiringWeightDown)*weight/event.get(h_prefiringWeightNominal));

  if(debug) cout << "Starting btagging..." << endl;
  // b-tagging

  // for the total, we need to quadratically add all >1 and all <1 values, respectively
  // so lets push all of them to a vector and do that afterwards
  std::vector<double> btag_variation_values_Up;
  std::vector<double> btag_variation_values_Down;

  hist("pt_ST_btagging_hfUp")->Fill(st, event.get(h_weight_btagdisc_hfUp)*weight/event.get(h_weight_btagdiscNominal));
  hist("pt_ST_btagging_hfDown")->Fill(st, event.get(h_weight_btagdisc_hfDown)*weight/event.get(h_weight_btagdiscNominal));
  btag_variation_values_Up.push_back(event.get(h_weight_btagdisc_hfUp) / event.get(h_weight_btagdiscNominal) - 1);
  btag_variation_values_Down.push_back(event.get(h_weight_btagdisc_hfDown) / event.get(h_weight_btagdiscNominal) - 1);

  hist("pt_ST_btagging_hfstats1Up")->Fill(st, event.get(h_weight_btagdisc_hfstats1Up)*weight/event.get(h_weight_btagdiscNominal));
  hist("pt_ST_btagging_hfstats1Down")->Fill(st, event.get(h_weight_btagdisc_hfstats1Down)*weight/event.get(h_weight_btagdiscNominal));
  btag_variation_values_Up.push_back(event.get(h_weight_btagdisc_hfstats1Up) / event.get(h_weight_btagdiscNominal) - 1);
  btag_variation_values_Down.push_back(event.get(h_weight_btagdisc_hfstats1Down) / event.get(h_weight_btagdiscNominal) - 1);

  hist("pt_ST_btagging_hfstats2Up")->Fill(st, event.get(h_weight_btagdisc_hfstats2Up)*weight/event.get(h_weight_btagdiscNominal));
  hist("pt_ST_btagging_hfstats2Down")->Fill(st, event.get(h_weight_btagdisc_hfstats2Down)*weight/event.get(h_weight_btagdiscNominal));
  btag_variation_values_Up.push_back(event.get(h_weight_btagdisc_hfstats2Up) / event.get(h_weight_btagdiscNominal) - 1);
  btag_variation_values_Down.push_back(event.get(h_weight_btagdisc_hfstats2Down) / event.get(h_weight_btagdiscNominal) - 1);

  hist("pt_ST_btagging_lfUp")->Fill(st, event.get(h_weight_btagdisc_lfUp)*weight/event.get(h_weight_btagdiscNominal));
  hist("pt_ST_btagging_lfDown")->Fill(st, event.get(h_weight_btagdisc_lfDown)*weight/event.get(h_weight_btagdiscNominal));
  btag_variation_values_Up.push_back(event.get(h_weight_btagdisc_lfUp) / event.get(h_weight_btagdiscNominal) - 1);
  btag_variation_values_Down.push_back(event.get(h_weight_btagdisc_lfDown) / event.get(h_weight_btagdiscNominal) - 1);

  hist("pt_ST_btagging_lfstats1Up")->Fill(st, event.get(h_weight_btagdisc_lfstats1Up)*weight/event.get(h_weight_btagdiscNominal));
  hist("pt_ST_btagging_lfstats1Down")->Fill(st, event.get(h_weight_btagdisc_lfstats1Down)*weight/event.get(h_weight_btagdiscNominal));
  btag_variation_values_Up.push_back(event.get(h_weight_btagdisc_lfstats1Up) / event.get(h_weight_btagdiscNominal) - 1);
  btag_variation_values_Down.push_back(event.get(h_weight_btagdisc_lfstats1Down) / event.get(h_weight_btagdiscNominal) - 1);

  hist("pt_ST_btagging_lfstats2Up")->Fill(st, event.get(h_weight_btagdisc_lfstats2Up)*weight/event.get(h_weight_btagdiscNominal));
  hist("pt_ST_btagging_lfstats2Down")->Fill(st, event.get(h_weight_btagdisc_lfstats2Down)*weight/event.get(h_weight_btagdiscNominal));
  btag_variation_values_Up.push_back(event.get(h_weight_btagdisc_lfstats2Up) / event.get(h_weight_btagdiscNominal) - 1);
  btag_variation_values_Down.push_back(event.get(h_weight_btagdisc_lfstats2Down) / event.get(h_weight_btagdiscNominal) - 1);

  hist("pt_ST_btagging_cferr1Up")->Fill(st, event.get(h_weight_btagdisc_cferr1Up)*weight/event.get(h_weight_btagdiscNominal));
  hist("pt_ST_btagging_cferr1Down")->Fill(st, event.get(h_weight_btagdisc_cferr1Down)*weight/event.get(h_weight_btagdiscNominal));
  btag_variation_values_Up.push_back(event.get(h_weight_btagdisc_cferr1Up) / event.get(h_weight_btagdiscNominal) - 1);
  btag_variation_values_Down.push_back(event.get(h_weight_btagdisc_cferr1Down) / event.get(h_weight_btagdiscNominal) - 1);

  hist("pt_ST_btagging_cferr2Up")->Fill(st, event.get(h_weight_btagdisc_cferr2Up)*weight/event.get(h_weight_btagdiscNominal));
  hist("pt_ST_btagging_cferr2Down")->Fill(st, event.get(h_weight_btagdisc_cferr2Down)*weight/event.get(h_weight_btagdiscNominal));
  btag_variation_values_Up.push_back(event.get(h_weight_btagdisc_cferr2Up) / event.get(h_weight_btagdiscNominal) - 1);
  btag_variation_values_Down.push_back(event.get(h_weight_btagdisc_cferr2Down) / event.get(h_weight_btagdiscNominal) - 1);

  // now lets do the total
  double total_greater = 0;
  double total_smaller = 0;
  for (uint i = 0; i < btag_variation_values_Up.size(); i++) {
    double valueUp = btag_variation_values_Up.at(i);
    double valueDown = btag_variation_values_Down.at(i);
    if(valueUp >= valueDown) {
      total_greater += valueUp*valueUp;
      total_smaller += valueDown*valueDown;
    } else {
      total_smaller += valueUp*valueUp;
      total_greater += valueDown*valueDown;
    }
  }
  hist("pt_ST_btagging_totalUp")->Fill(st, weight * (1 + sqrt(total_greater)));
  hist("pt_ST_btagging_totalDown")->Fill(st, weight * (1 - sqrt(total_smaller)));

  // b-tagging yield
  if(event.get(h_flag_muonevent)) {
    // fetching the event weight based on muon pt
    double this_weight = btagyield_mu->GetBinContent( btagyield_mu->GetXaxis()->FindBin(event.muons->at(0).pt()) );
    hist("pt_ST_btagYieldUp")->Fill(st, weight * this_weight);
    hist("pt_ST_btagYieldDown")->Fill(st, weight / this_weight);
  } else {
    double this_weight = btagyield_ele->GetBinContent( btagyield_ele->GetXaxis()->FindBin(event.electrons->at(0).pt()) );
    hist("pt_ST_btagYieldUp")->Fill(st, weight * this_weight);
    hist("pt_ST_btagYieldDown")->Fill(st, weight / this_weight);
  }

  if(debug) cout << "Starting ele..." << endl;
  // sfelec_id
  hist("pt_ST_sfelec_idUp")->Fill(st, event.get(h_weight_sfelec_idUp)*weight/event.get(h_weight_sfelec_idNominal));
  hist("pt_ST_sfelec_idDown")->Fill(st, event.get(h_weight_sfelec_idDown)*weight/event.get(h_weight_sfelec_idNominal));

  // sfelec_trigger
  if( event.get(h_weight_sfelec_triggerNominal) != 0) {
    // some entries in h_weight_sfelec_triggerNominal can be 0 -> we do not want to include those! 
    hist("pt_ST_sfelec_triggerUp")->Fill(st, event.get(h_weight_sfelec_triggerUp)*weight/event.get(h_weight_sfelec_triggerNominal));
    hist("pt_ST_sfelec_triggerDown")->Fill(st, event.get(h_weight_sfelec_triggerDown)*weight/event.get(h_weight_sfelec_triggerNominal));
  }
  
  // sfelec_reco
  hist("pt_ST_sfelec_recoUp")->Fill(st, event.get(h_weight_sfelec_recoUp)*weight/event.get(h_weight_sfelec_recoNominal));
  hist("pt_ST_sfelec_recoDown")->Fill(st, event.get(h_weight_sfelec_recoDown)*weight/event.get(h_weight_sfelec_recoNominal));

  if(debug) cout << "Starting mu..." << endl;
  // sfmu_id
  hist("pt_ST_sfmu_idUp")->Fill(st, event.get(h_weight_sfmu_idUp)*weight/event.get(h_weight_sfmu_idNominal));
  hist("pt_ST_sfmu_idDown")->Fill(st, event.get(h_weight_sfmu_idDown)*weight/event.get(h_weight_sfmu_idNominal));

  // sfmu_iso
  hist("pt_ST_sfmu_isoUp")->Fill(st, event.get(h_weight_sfmu_isoUp)*weight/event.get(h_weight_sfmu_isoNominal));
  hist("pt_ST_sfmu_isoDown")->Fill(st, event.get(h_weight_sfmu_isoDown)*weight/event.get(h_weight_sfmu_isoNominal));

  // sfmu_iso
  hist("pt_ST_sfmu_triggerUp")->Fill(st, event.get(h_weight_sfmu_triggerUp)*weight/event.get(h_weight_sfmu_triggerNominal));
  hist("pt_ST_sfmu_triggerDown")->Fill(st, event.get(h_weight_sfmu_triggerDown)*weight/event.get(h_weight_sfmu_triggerNominal));

  // decorrelation uncertainty
  if(debug) cout << "Starting decorrelation..." << endl;
  double decorr_uncertainty = decorrelationUncertainty->GetBinContent( decorrelationUncertainty->FindBin(st) );
  hist("pt_ST_decorrelationUp")->Fill(st, weight + (decorr_uncertainty * weight));
  hist("pt_ST_decorrelationDown")->Fill(st, weight - (decorr_uncertainty * weight));

  if(debug) cout << "Starting systweights..." << endl;
  float orig_weight = event.genInfo->originalXWGTUP();
  int MY_FIRST_INDEX = 9;
  if(needsOtherMCweightHandling) MY_FIRST_INDEX = 47;
  for(int i=0; i<100; i++){
    double pdf_weight = 0;
    if(event.genInfo->systweights().size() > 0) pdf_weight =event.genInfo->systweights().at(i+MY_FIRST_INDEX);
    const char* name = hist_names[i].c_str();
    hist(name)->Fill(st, weight * pdf_weight / orig_weight);
  }

  try {
    float murmuf_upup        = event.get(h_murmuf_upup);
    float murmuf_upnone      = event.get(h_murmuf_upnone);
    float murmuf_noneup      = event.get(h_murmuf_noneup);
    float murmuf_nonedown    = event.get(h_murmuf_nonedown);
    float murmuf_downnone    = event.get(h_murmuf_downnone);
    float murmuf_downdown    = event.get(h_murmuf_downdown);

    hist("pt_ST_murmuf_upup")->Fill(st, weight * murmuf_upup);
    hist("pt_ST_murmuf_upnone")->Fill(st, weight * murmuf_upnone);
    hist("pt_ST_murmuf_noneup")->Fill(st, weight * murmuf_noneup);
    hist("pt_ST_murmuf_nonedown")->Fill(st, weight * murmuf_nonedown);
    hist("pt_ST_murmuf_downnone")->Fill(st, weight * murmuf_downnone);
    hist("pt_ST_murmuf_downdown")->Fill(st, weight * murmuf_downdown);
  } catch (...) { }

  if(debug) cout << "Done with plotting" << endl;

}

TstarTstarSignalRegionHists::~TstarTstarSignalRegionHists(){}
