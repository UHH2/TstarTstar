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
  h_ST = ctx.get_handle<double>("STHOTVR");
  is_MC = ctx.get("dataset_type") == "MC";

  if(ctx.get("debug", "<not set>") == "true") debug = true;

  needsOtherMCweightHandling = ctx.get("dataset_version").find("TstarTstar") != std::string::npos;
  if(needsOtherMCweightHandling) std::cout << "We are signal so we need other idices" << std::endl;

  // Histogram containing decorrelation uncertainty
  TFile *decorrelationUncertaintyFile = new TFile("/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/decorrelationComparison_total_smooth.root");
  decorrelationUncertainty = (TH1*)decorrelationUncertaintyFile->Get("decorrelation_uncertainty");

  // weight handles

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
  double bins[nbins] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500,
    2600, 2700, 2800, 2900, 3000, 3250, 4000, 6000};

  // nominal histogram
  book<TH1F>("pt_ST_nominal", "S_{T} [GeV]", nbins-1, bins);

  // variation histograms
  book<TH1F>("pt_ST_puUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_puDown", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_prefiringUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_prefiringDown", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_btagging_totalUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagging_totalDown", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_btagging_hfUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagging_hfDown", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_btagging_hfstats1Up", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagging_hfstats1Down", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_btagging_hfstats2Up", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagging_hfstats2Down", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_btagging_lfUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagging_lfDown", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_btagging_lfstats1Up", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagging_lfstats1Down", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_btagging_lfstats2Up", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagging_lfstats2Down", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_btagging_cferr1Up", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagging_cferr1Down", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_btagging_cferr2Up", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagging_cferr2Down", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_sfelec_idUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_sfelec_idDown", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_sfelec_triggerUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_sfelec_triggerDown", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_sfelec_recoUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_sfelec_recoDown", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_sfmu_idUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_sfmu_idDown", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_sfmu_isoUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_sfmu_isoDown", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_sfmu_triggerUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_sfmu_triggerDown", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_ST_decorrelationUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_decorrelationDown", "S_{T} [GeV]", nbins-1, bins);

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

  double st = event.get(h_ST);
  if(st > 6000) st = 5999.9; // handling overflow

  // fill nominal
  hist("pt_ST_nominal")->Fill(st, weight);

  if(!is_MC) return;

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
  for (int i = 0; i < btag_variation_values_Up.size(); i++) {
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

  if(debug) cout << "Starting btagging..." << endl;
  float orig_weight = event.genInfo->originalXWGTUP();
  int MY_FIRST_INDEX = 9;
  if(needsOtherMCweightHandling) MY_FIRST_INDEX = 41;
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



}

TstarTstarSignalRegionHists::~TstarTstarSignalRegionHists(){}
