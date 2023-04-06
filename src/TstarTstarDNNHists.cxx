#include "UHH2/TstarTstar/include/TstarTstarDNNHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <string>

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

TstarTstarDNNHists::TstarTstarDNNHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  h_DNN_output = ctx.get_handle<double>("DNN_output");
  h_newTagger = ctx.get_handle<double>("newTagger");
  h_DNN_Inputs = ctx.get_handle<std::vector<double>>("DNN_Inputs");
  h_ST = ctx.get_handle<double>("STHOTVR");

  // book all histograms here
  book<TH1F>("DNN_output", "DNN output", 20, 0, 1);
  book<TH1F>("newTagger", "DDT score", 20, -1, 1);
  book<TH1F>("DNN_output_noWeights", "DNN output NO WEIGHTS", 20, 0, 1);
  DNN_2D_ST = book<TH2D>("2D_DNN_ST", "DNN output against ST", 40, 0, 4000, 50, 0, 1);
  newTagger_2D_ST = book<TH2D>("2D_newTagger_ST", "DNN_{DDT} against ST", 40, 0, 4000, 50, -1, 1);
  book<TH1F>("DNN_output_nolowST", "DNN output nolowST", 20, 0, 1);

  DNN_2D_1 = book<TH2D>("2D_DNN_1", "DNN output against lepton pt", 50, -1, 1, 50, 0, 1);
  DNN_2D_2 = book<TH2D>("2D_DNN_2", "DNN output against lepton eta", 50, -1, 1, 50, 0, 1);
  DNN_2D_3 = book<TH2D>("2D_DNN_3", "DNN output against lepton phi", 50, -1, 1, 50, 0, 1);
  DNN_2D_4 = book<TH2D>("2D_DNN_4", "DNN output against lepton relIso", 50, -1, 1, 50, 0, 1);
  DNN_2D_5 = book<TH2D>("2D_DNN_5", "DNN output against HOTVR-1 pt", 50, -1, 1, 50, 0, 1);
  DNN_2D_6 = book<TH2D>("2D_DNN_6", "DNN output against HOTVR-1 eta", 50, -1, 1, 50, 0, 1);
  DNN_2D_7 = book<TH2D>("2D_DNN_7", "DNN output against HOTVR-1 phi", 50, -1, 1, 50, 0, 1);
  DNN_2D_8 = book<TH2D>("2D_DNN_8", "DNN output against HOTVR-1 tau1", 50, -1, 1, 50, 0, 1);
  DNN_2D_9 = book<TH2D>("2D_DNN_9", "DNN output against HOTVR-1 tau2", 50, -1, 1, 50, 0, 1);
  DNN_2D_10 = book<TH2D>("2D_DNN_10", "DNN output against HOTVR-1 tau3", 50, -1, 1, 50, 0, 1);
  DNN_2D_11 = book<TH2D>("2D_DNN_11", "DNN output against HOTVR-1 subjets", 50, -1, 1, 50, 0, 1);
  DNN_2D_12 = book<TH2D>("2D_DNN_12", "DNN output against HOTVR-2 pt", 50, -1, 1, 50, 0, 1);
  DNN_2D_13 = book<TH2D>("2D_DNN_13", "DNN output against HOTVR-2 eta", 50, -1, 1, 50, 0, 1);
  DNN_2D_14 = book<TH2D>("2D_DNN_14", "DNN output against HOTVR-2 phi", 50, -1, 1, 50, 0, 1);
  DNN_2D_15 = book<TH2D>("2D_DNN_15", "DNN output against HOTVR-2 tau1", 50, -1, 1, 50, 0, 1);
  DNN_2D_16 = book<TH2D>("2D_DNN_16", "DNN output against HOTVR-2 tau2", 50, -1, 1, 50, 0, 1);
  DNN_2D_17 = book<TH2D>("2D_DNN_17", "DNN output against HOTVR-2 tau3", 50, -1, 1, 50, 0, 1);
  DNN_2D_18 = book<TH2D>("2D_DNN_18", "DNN output against HOTVR-2 subjets", 50, -1, 1, 50, 0, 1);
  DNN_2D_19 = book<TH2D>("2D_DNN_19", "DNN output against HOTVR-3 pt", 50, -1, 1, 50, 0, 1);
  DNN_2D_20 = book<TH2D>("2D_DNN_20", "DNN output against HOTVR-3 eta", 50, -1, 1, 50, 0, 1);
  DNN_2D_21 = book<TH2D>("2D_DNN_21", "DNN output against HOTVR-3 phi", 50, -1, 1, 50, 0, 1);
  DNN_2D_22 = book<TH2D>("2D_DNN_22", "DNN output against HOTVR-3 tau1", 50, -1, 1, 50, 0, 1);
  DNN_2D_23 = book<TH2D>("2D_DNN_23", "DNN output against HOTVR-3 tau2", 50, -1, 1, 50, 0, 1);
  DNN_2D_24 = book<TH2D>("2D_DNN_24", "DNN output against HOTVR-3 tau3", 50, -1, 1, 50, 0, 1);
  DNN_2D_25 = book<TH2D>("2D_DNN_25", "DNN output against HOTVR-3 subjets", 50, -1, 1, 50, 0, 1);
  DNN_2D_26 = book<TH2D>("2D_DNN_26", "DNN output against b-jet pt", 50, -1, 1, 50, 0, 1);
  DNN_2D_27 = book<TH2D>("2D_DNN_27", "DNN output against b-jet eta", 50, -1, 1, 50, 0, 1);
  DNN_2D_28 = book<TH2D>("2D_DNN_28", "DNN output against b-jet phi", 50, -1, 1, 50, 0, 1);
  DNN_2D_29 = book<TH2D>("2D_DNN_29", "DNN output against b-jet b-tag", 50, -1, 1, 50, 0, 1);
  DNN_2D_30 = book<TH2D>("2D_DNN_30", "DNN output against MET pt", 50, -1, 1, 50, 0, 1);
  DNN_2D_31 = book<TH2D>("2D_DNN_31", "DNN output against MET phi", 50, -1, 1, 50, 0, 1);
  DNN_2D_32 = book<TH2D>("2D_DNN_32", "DNN output against N AK4", 50, -1, 1, 50, 0, 1);
  DNN_2D_33 = book<TH2D>("2D_DNN_33", "DNN output against N HOTVR", 50, -1, 1, 50, 0, 1);

}


void TstarTstarDNNHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  double DNNoutput = event.get(h_DNN_output);
  std::vector<double> inputs = event.get(h_DNN_Inputs);

  hist("DNN_output")->Fill(DNNoutput, weight);
  hist("DNN_output_noWeights")->Fill(DNNoutput, 1);
  double st = event.get(h_ST);
  DNN_2D_ST->Fill(st, DNNoutput, weight);
  if(st > 500) hist("DNN_output_nolowST")->Fill(DNNoutput, weight);

  DNN_2D_1->Fill(inputs.at(0), DNNoutput, weight);
  DNN_2D_2->Fill(inputs.at(1), DNNoutput, weight);
  DNN_2D_3->Fill(inputs.at(2), DNNoutput, weight);
  DNN_2D_4->Fill(inputs.at(3), DNNoutput, weight);
  DNN_2D_5->Fill(inputs.at(4), DNNoutput, weight);
  DNN_2D_6->Fill(inputs.at(5), DNNoutput, weight);
  DNN_2D_7->Fill(inputs.at(6), DNNoutput, weight);
  DNN_2D_8->Fill(inputs.at(7), DNNoutput, weight);
  DNN_2D_9->Fill(inputs.at(8), DNNoutput, weight);
  DNN_2D_10->Fill(inputs.at(9), DNNoutput, weight);
  DNN_2D_11->Fill(inputs.at(10), DNNoutput, weight);
  DNN_2D_12->Fill(inputs.at(11), DNNoutput, weight);
  DNN_2D_13->Fill(inputs.at(12), DNNoutput, weight);
  DNN_2D_14->Fill(inputs.at(13), DNNoutput, weight);
  DNN_2D_15->Fill(inputs.at(14), DNNoutput, weight);
  DNN_2D_16->Fill(inputs.at(15), DNNoutput, weight);
  DNN_2D_17->Fill(inputs.at(16), DNNoutput, weight);
  DNN_2D_18->Fill(inputs.at(17), DNNoutput, weight);
  DNN_2D_19->Fill(inputs.at(18), DNNoutput, weight);
  DNN_2D_20->Fill(inputs.at(19), DNNoutput, weight);
  DNN_2D_21->Fill(inputs.at(20), DNNoutput, weight);
  DNN_2D_22->Fill(inputs.at(21), DNNoutput, weight);
  DNN_2D_23->Fill(inputs.at(22), DNNoutput, weight);
  DNN_2D_24->Fill(inputs.at(23), DNNoutput, weight);
  DNN_2D_25->Fill(inputs.at(24), DNNoutput, weight);
  DNN_2D_26->Fill(inputs.at(25), DNNoutput, weight);
  DNN_2D_27->Fill(inputs.at(26), DNNoutput, weight);
  DNN_2D_28->Fill(inputs.at(27), DNNoutput, weight);
  DNN_2D_29->Fill(inputs.at(28), DNNoutput, weight);
  DNN_2D_30->Fill(inputs.at(29), DNNoutput, weight);
  DNN_2D_31->Fill(inputs.at(30), DNNoutput, weight);
  DNN_2D_32->Fill(inputs.at(31), DNNoutput, weight);
  DNN_2D_33->Fill(inputs.at(32), DNNoutput, weight);

  try {
    double newTagger = event.get(h_newTagger);
    hist("newTagger")->Fill(newTagger, weight);
    newTagger_2D_ST->Fill(st, newTagger, weight);
  } catch(...) {}

}

TstarTstarDNNHists::~TstarTstarDNNHists(){}
