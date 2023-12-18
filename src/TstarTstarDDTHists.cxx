#include "UHH2/TstarTstar/include/TstarTstarDDTHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

TstarTstarDDTHists::TstarTstarDDTHists(Context & ctx, const string & dirname, const std::vector<TString> points): Hists(ctx, dirname){

  // ST handle
  h_ST_HOTVR = ctx.get_handle<double>("ST_HOTVR");

  points_ = points;

  // binning definition
  const int nbins = 34;
  double bins[nbins] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500,
    2600, 2700, 2800, 2900, 3000, 3250, 4000, 6000};

  // nominal histogram
  for (auto point : points_) {
    book<TH1F>("pt_ST_"+point, "S_{T} [GeV]", nbins-1, bins);
    book<TH1F>("not_pt_ST_"+point, "S_{T} [GeV]", nbins-1, bins);
  }

}


void TstarTstarDDTHists::fill(const Event & event, std::vector<double> taggerScores){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  bool debug = false;

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  if(debug) cout << "Starting Tstar Hists." << endl;

  double st = event.get(h_ST_HOTVR);
  if(st > 6000) st = 5999.9; // handling overflow

  int i = 0;
  for (auto point : points_) {
      if( taggerScores.at(i) > 0 ) hist("pt_ST_"+point)->Fill(st, weight);
      else hist("not_pt_ST_"+point)->Fill(st, weight);
      i++;
  }

}

TstarTstarDDTHists::~TstarTstarDDTHists(){}
