#include "UHH2/TstarTstar/include/TstarTstarDNNInputHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <string>
#include <math.h>

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

TstarTstarDNNInputHists::TstarTstarDNNInputHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  h_DNN_Inputs = ctx.get_handle<std::vector<double>>("DNN_Inputs");
  h_DNN_AddInputs = ctx.get_handle<std::vector<double>>("DNN_AddInputs");


  // book all histograms here
  DNN_2D_1 = book<TH1D>("DNN_Input_1", "lepton p_{T}",  20, 0, 2000);
  DNN_2D_2 = book<TH1D>("DNN_Input_2", "lepton #eta",  20, -4, 4);
  DNN_2D_3 = book<TH1D>("DNN_Input_3", "lepton #phi",  20, -4, 4);
  DNN_2D_4 = book<TH1D>("DNN_Input_4", "lepton Iso_{rel}",  20, 0, 2);
  DNN_2D_5 = book<TH1D>("DNN_Input_5", "HOTVR-1 p_{T}",  20, 0, 3000);
  DNN_2D_6 = book<TH1D>("DNN_Input_6", "HOTVR-1 #eta",  20, -4, 4);
  DNN_2D_7 = book<TH1D>("DNN_Input_7", "HOTVR-1 #phi",  20, -4, 4);
  DNN_2D_8 = book<TH1D>("DNN_Input_8", "HOTVR-1 #tau_{1}",  20, 0, 1);
  DNN_2D_9 = book<TH1D>("DNN_Input_9", "HOTVR-1 #tau_{2}",  20, 0, 1);
  DNN_2D_10 = book<TH1D>("DNN_Input_10", "HOTVR-1 #tau_{3}",  20, 0, 1);
  DNN_2D_11 = book<TH1D>("DNN_Input_11", "HOTVR-1 N_{subjets}",  20, 0, 5);
  DNN_2D_12 = book<TH1D>("DNN_Input_12", "HOTVR-2 p_{T}",  20, 0, 3000);
  DNN_2D_13 = book<TH1D>("DNN_Input_13", "HOTVR-2 #eta",  20, -4, 4);
  DNN_2D_14 = book<TH1D>("DNN_Input_14", "HOTVR-2 #phi",  20, -4, 4);
  DNN_2D_15 = book<TH1D>("DNN_Input_15", "HOTVR-2 #tau_{1}",  20, 0, 1);
  DNN_2D_16 = book<TH1D>("DNN_Input_16", "HOTVR-2 #tau_{2}",  20, 0, 1);
  DNN_2D_17 = book<TH1D>("DNN_Input_17", "HOTVR-2 #tau_{3}",  20, 0, 1);
  DNN_2D_18 = book<TH1D>("DNN_Input_18", "HOTVR-2 N_{subjets}",  20, 0, 5);
  DNN_2D_19 = book<TH1D>("DNN_Input_19", "HOTVR-3 p_{T}",  20, 0, 3000);
  DNN_2D_20 = book<TH1D>("DNN_Input_20", "HOTVR-3 #eta",  20, -4, 4);
  DNN_2D_21 = book<TH1D>("DNN_Input_21", "HOTVR-3 #phi",  20, -4, 4);
  DNN_2D_22 = book<TH1D>("DNN_Input_22", "HOTVR-3 #tau_{1}",  20, 0, 1);
  DNN_2D_23 = book<TH1D>("DNN_Input_23", "HOTVR-3 #tau_{2}",  20, 0, 1);
  DNN_2D_24 = book<TH1D>("DNN_Input_24", "HOTVR-3 #tau_{3}",  20, 0, 1);
  DNN_2D_25 = book<TH1D>("DNN_Input_25", "HOTVR-3 N_{subjets}",  20, 0, 5);
  DNN_2D_26 = book<TH1D>("DNN_Input_26", "b-jet p_{T}",  20, 0, 1000);
  DNN_2D_27 = book<TH1D>("DNN_Input_27", "b-jet #eta",  20, -4, 4);
  DNN_2D_28 = book<TH1D>("DNN_Input_28", "b-jet #phi",  20, -4, 4);
  DNN_2D_29 = book<TH1D>("DNN_Input_29", "b-jet DeepJet score",  20, 0, 1);
  DNN_2D_30 = book<TH1D>("DNN_Input_30", "MET p_{T}",  20, 0, 2000);
  DNN_2D_31 = book<TH1D>("DNN_Input_31", "MET #phi",  20, -4, 4);
  DNN_2D_32 = book<TH1D>("DNN_Input_32", "N_{AK4}",  20, 0, 20);
  DNN_2D_33 = book<TH1D>("DNN_Input_33", "N_{HOTVR}",  10, 0, 10);

  AddDNN_2D_1 = book<TH1D>("DNN_AddInput_1", "p_{T} asymmetry",  25, 0, 1);
  AddDNN_2D_2 = book<TH1D>("DNN_AddInput_2", "#Delta R (HOTVR 1, HOTVR 2)",  30, 0, 6);
  AddDNN_2D_3 = book<TH1D>("DNN_AddInput_3", "#Delta R (lepton, closest HOTVR)",  30, 0, 6);

}


void TstarTstarDNNInputHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  std::vector<double> inputs = event.get(h_DNN_Inputs);
  std::vector<double> addInputs = event.get(h_DNN_AddInputs);

  DNN_2D_1->Fill(inputs.at(0), weight);
  DNN_2D_2->Fill(inputs.at(1), weight);
  DNN_2D_3->Fill(inputs.at(2), weight);
  DNN_2D_4->Fill(inputs.at(3), weight);
  DNN_2D_5->Fill(inputs.at(4), weight);
  DNN_2D_6->Fill(inputs.at(5), weight);
  DNN_2D_7->Fill(inputs.at(6), weight);
  DNN_2D_8->Fill(inputs.at(7), weight);
  DNN_2D_9->Fill(inputs.at(8), weight);
  DNN_2D_10->Fill(inputs.at(9), weight);
  DNN_2D_11->Fill(inputs.at(10), weight);
  DNN_2D_12->Fill(inputs.at(11), weight);
  DNN_2D_13->Fill(inputs.at(12), weight);
  DNN_2D_14->Fill(inputs.at(13), weight);
  DNN_2D_15->Fill(inputs.at(14), weight);
  DNN_2D_16->Fill(inputs.at(15), weight);
  DNN_2D_17->Fill(inputs.at(16), weight);
  DNN_2D_18->Fill(inputs.at(17), weight);
  DNN_2D_19->Fill(inputs.at(18), weight);
  DNN_2D_20->Fill(inputs.at(19), weight);
  DNN_2D_21->Fill(inputs.at(20), weight);
  DNN_2D_22->Fill(inputs.at(21), weight);
  DNN_2D_23->Fill(inputs.at(22), weight);
  DNN_2D_24->Fill(inputs.at(23), weight);
  DNN_2D_25->Fill(inputs.at(24), weight);
  DNN_2D_26->Fill(inputs.at(25), weight);
  DNN_2D_27->Fill(inputs.at(26), weight);
  DNN_2D_28->Fill(inputs.at(27), weight);
  DNN_2D_29->Fill(inputs.at(28), weight);
  DNN_2D_30->Fill(inputs.at(29), weight);
  DNN_2D_31->Fill(inputs.at(30), weight);
  DNN_2D_32->Fill(inputs.at(31), weight);
  DNN_2D_33->Fill(inputs.at(32), weight);

  AddDNN_2D_1->Fill(addInputs.at(0), weight);
  AddDNN_2D_2->Fill(addInputs.at(1), weight);
  AddDNN_2D_3->Fill(addInputs.at(2), weight);

}

TstarTstarDNNInputHists::~TstarTstarDNNInputHists(){}
