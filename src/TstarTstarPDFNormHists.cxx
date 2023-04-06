#include "UHH2/TstarTstar/include/TstarTstarPDFNormHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

TstarTstarPDFNormHists::TstarTstarPDFNormHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  
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


}


void TstarTstarPDFNormHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  bool debug = false;

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  hist("nominal")->Fill(0.5,weight);

  float orig_weight = event.genInfo->originalXWGTUP();
  int MY_FIRST_INDEX = 9;
  for(int i=0; i<100; i++){
    double pdf_weight = 0;
    if(event.genInfo->systweights().size() > 0) pdf_weight =event.genInfo->systweights().at(i+MY_FIRST_INDEX);
    const char* name = hist_names[i].c_str();
    hist(name)->Fill(0.5, weight * pdf_weight / orig_weight);
  }

}

TstarTstarPDFNormHists::~TstarTstarPDFNormHists(){}
