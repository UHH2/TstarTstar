#include "UHH2/TstarTstar/include/TstarTstarPDFNormHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

TstarTstarPDFNormHists::TstarTstarPDFNormHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  needsOtherMCweightHandling = ctx.get("dataset_version").find("TstarTstar") != std::string::npos;
  //if(needsOtherMCweightHandling) std::cout << "We are signal so we need other idices" << std::endl;
  
  // weight handles
  book<TH1F>("nominal", "nominal", 1, 0, 1);

  // 100 histograms for the PDF stuff
  for(int i=0; i<100; i++){
    std::stringstream ss_name;
    ss_name << "PDF_" << i+1;

    stringstream ss_title;
    ss_title << "S_T [GeV] for PDF No. "  << i+1 << " out of 100" ;

    std::string s_name = ss_name.str();
    std::string s_title = ss_title.str();
    const char* char_name = s_name.c_str();
    const char* char_title = s_title.c_str();

    hist_names[i] = s_name;

    book<TH1F>(char_name, char_title, 1, 0, 1);

  }

  h_murmuf_upup        = ctx.get_handle< float >("weight_murmuf_upup");
  h_murmuf_upnone      = ctx.get_handle< float >("weight_murmuf_upnone");
  h_murmuf_noneup      = ctx.get_handle< float >("weight_murmuf_noneup");
  h_murmuf_nonedown    = ctx.get_handle< float >("weight_murmuf_nonedown");
  h_murmuf_downnone    = ctx.get_handle< float >("weight_murmuf_downnone");
  h_murmuf_downdown    = ctx.get_handle< float >("weight_murmuf_downdown");

  book<TH1F>("murmuf_upup", "murmuf_upup", 1, 0, 1);
  book<TH1F>("murmuf_upnone", "murmuf_upnone", 1, 0, 1);
  book<TH1F>("murmuf_noneup", "murmuf_noneup", 1, 0, 1);
  book<TH1F>("murmuf_nonedown", "murmuf_nonedown", 1, 0, 1);
  book<TH1F>("murmuf_downnone", "murmuf_downnone", 1, 0, 1);
  book<TH1F>("murmuf_downdown", "murmuf_downdown", 1, 0, 1);


}


void TstarTstarPDFNormHists::fill(const Event & event){

  bool debug = false;

  double weight = event.weight;

  hist("nominal")->Fill(0.5,weight);

  float orig_weight = event.genInfo->originalXWGTUP();
  int MY_FIRST_INDEX = 9;
  if(needsOtherMCweightHandling) MY_FIRST_INDEX = 47;
  for(int i=0; i<100; i++){
    double pdf_weight = 0;
    if(event.genInfo->systweights().size() > 0) pdf_weight =event.genInfo->systweights().at(i+MY_FIRST_INDEX);
    const char* name = hist_names[i].c_str();
    hist(name)->Fill(0.5, weight * pdf_weight / orig_weight);
  }

  // fill for MC scale :)
  float murmuf_upup        = event.get(h_murmuf_upup);
  float murmuf_upnone      = event.get(h_murmuf_upnone);
  float murmuf_noneup      = event.get(h_murmuf_noneup);
  float murmuf_nonedown    = event.get(h_murmuf_nonedown);
  float murmuf_downnone    = event.get(h_murmuf_downnone);
  float murmuf_downdown    = event.get(h_murmuf_downdown);
  hist("murmuf_upup")->Fill(0.5, weight * murmuf_upup);
  hist("murmuf_upnone")->Fill(0.5, weight * murmuf_upnone);
  hist("murmuf_noneup")->Fill(0.5, weight * murmuf_noneup);
  hist("murmuf_nonedown")->Fill(0.5, weight * murmuf_nonedown);
  hist("murmuf_downnone")->Fill(0.5, weight * murmuf_downnone);
  hist("murmuf_downdown")->Fill(0.5, weight * murmuf_downdown);

}

TstarTstarPDFNormHists::~TstarTstarPDFNormHists(){}
