#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>
#include <string>

// header
#include "UHH2/TstarTstar/include/NeuralNetworkModules.h"

// UHH2 stuff
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/AdditionalSelections.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TTbarReconstruction.h"


using namespace std;

// quick inv mass method
float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

// A quick functions to get the DNN mean and std for normalisation
// a vector containing the values is returned
std::vector<double> get_DNN_vals(const std::string& file) {
  std::vector<double> vals;
  std::ifstream infile(file);
  std::string line;
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    vals.push_back(std::stod(line));
  }
  return vals;
}

NeuralNetworkModule::NeuralNetworkModule(Context& ctx, const std::string & ModelName, const std::string& ConfigName): NeuralNetworkBase(ctx, ModelName, ConfigName) {}

bool NeuralNetworkModule::process(Event & event, std::vector<double> values){
  CreateInputs(event, values);
  std::vector<std::string> inputs;
  inputs.reserve(LayerInputs.size());
  for (auto const& imap: LayerInputs) inputs.push_back(imap.first);

  tensorflow::run(NNsession, inputs, NNInputs, LayerOutputs, &NNoutputs);
  return true;
}

void NeuralNetworkModule::CreateInputs(Event & event, std::vector<double> values){
  // clear possible previous values
  NNInputs.clear();
  NNoutputs.clear();
  // add tensors to vector
  NNInputs.push_back( tensorflow::Tensor(tensorflow::DT_FLOAT, {1, values.size()}));
  for(uint i = 0; i < values.size(); i++){
    NNInputs.at(0).tensor<float, 2>()(0,i) = values.at(i);
  }

  // sanity check when including
  if (NNInputs.size()!=LayerInputs.size()) throw logic_error("NeuralNetworkIncluder.cxx: Create a number of inputs diffetent wrt. LayerInputs.size()="+to_string(LayerInputs.size())+" NNInputs.size()="+to_string(NNInputs.size()));
}

// ######################################
// ######################################
// ######################################

NeuralNetworkInputCreator::NeuralNetworkInputCreator(Context& ctx) {
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_do_masspoint = ctx.get_handle<bool>("do_masspoint");
  h_is_muevt = ctx.get_handle<bool>("is_muevt");

  h_CHS_matched = ctx.get_handle<vector<Jet>>("CHS_matched");

  // Check if MC or Data
  is_MC = ctx.get("dataset_type") == "MC";
  // Init random for masspoint
  std::srand(std::time(nullptr)); // Initialize random
}

std::vector<double> NeuralNetworkInputCreator::getInputs() {
  return DNNInputs;
}

std::vector<double> NeuralNetworkInputCreator::getAddInputs() {
  return AddDNNInputs;
}

bool NeuralNetworkInputCreator::setInputs(std::vector<double> vec) {
  DNNInputs = vec;
  return true;
}

bool NeuralNetworkInputCreator::setAddInputs(std::vector<double> vec) {
  AddDNNInputs = vec;
  return true;
}

bool NeuralNetworkInputCreator::createInputs(Event& event) {

  // push values to vector
  const FlavorParticle& lepton = event.get(h_primlep);
  std::vector<double> values;

  // Lepton
  values.push_back(lepton.pt());
  values.push_back(lepton.eta());
  values.push_back(lepton.phi());
  if(event.get(h_is_muevt)) values.push_back(event.muons->at(0).relIso());
  else values.push_back(event.electrons->at(0).relIso());

  // HOTVR jets
  uint i = 0;
  for (const auto & topjet : *event.topjets){
    if(i == 3) break;
    values.push_back(topjet.pt());
    values.push_back(topjet.eta());
    values.push_back(topjet.phi());
    values.push_back(topjet.tau1_groomed());
    values.push_back(topjet.tau2_groomed());
    values.push_back(topjet.tau3_groomed());
    values.push_back(topjet.subjets().size());
    i++;
  }
  for(;i < 3;i++){
    values.push_back(0); // pt
    values.push_back(3.5); // eta
    values.push_back(4); // phi
    values.push_back(0); // tau1
    values.push_back(0); // tau2
    values.push_back(0); // tau3
    values.push_back(0); // subjets
  }


  // we'll need to redo the matching here as indices have changed in the selection
  // this should be fixed there in the next iteration
  // but it is not fixed right now, so this fix is needed
  double maxbtag = 0.;
  Jet btaggedjet_PUPPI;
  Jet btaggedjet_CHS;
  for (int i = 0; i < event.jets->size(); i++) { // loop over PUPPI jets (as we are using these)

    // find matching CHS jet
    Jet this_jet = event.jets->at(i);
    Jet matched_jet;

    double deltaR_min = 99;
    for(const Jet & CHSjet : event.get(h_CHS_matched)){ // CHS jets
      double deltaR_CHS = deltaR(this_jet, CHSjet);
      if(deltaR_CHS < deltaR_min) {
        deltaR_min = deltaR_CHS;
        matched_jet = CHSjet;
      }
    } // end CHS loop

    double this_btag = matched_jet.btag_DeepJet();
    if(this_btag > maxbtag) {
      btaggedjet_PUPPI = this_jet;
      btaggedjet_CHS = matched_jet;
    }
    
  }
  values.push_back(btaggedjet_PUPPI.pt());
  values.push_back(btaggedjet_PUPPI.eta());
  values.push_back(btaggedjet_PUPPI.phi());
  values.push_back(btaggedjet_CHS.btag_DeepJet());

  // MET
  values.push_back(event.met->pt());
  values.push_back(event.met->phi());

  // Event variables
  values.push_back(event.jets->size());
  values.push_back(event.topjets->size());

  if(values.size() != 33) {
    std::cout << "You are an idiot that can not count! " << values.size() << endl;
    throw "You are an idiot that can not count!";
  }

  DNNInputs = values;

  return true;
}

bool NeuralNetworkInputCreator::createAddInputs(Event& event) {

  std::vector<double> values;

  // first add input: p_T asymmetry jet 1 jet 4
  double pt_asym = event.jets->at(0).pt()-event.jets->at(3).pt();
  double pt_sum = 0.;
  for (const auto & jet : *event.jets) pt_sum += jet.pt();
  pt_asym /= pt_sum;
  values.push_back(pt_asym);

  // second add input: deltaR leading 2 HOTVR jets
  if(event.topjets->size() > 1) values.push_back(deltaR(event.topjets->at(0), event.topjets->at(1)));
  else values.push_back(-1);

  // third add input: deltaR lepton closest HOTVR jet
  const FlavorParticle& lepton = event.get(h_primlep);
  double min_deltaR = 9999;
  for (const auto & topjet : *event.topjets) {
    double cur_deltaR = deltaR(lepton, topjet);
    if(cur_deltaR < min_deltaR) min_deltaR = cur_deltaR;
  }
  values.push_back(min_deltaR);

  AddDNNInputs = values;

  return true;
}

// ######################################
// ######################################
// ######################################

NeuralNetworkInputNormalizer::NeuralNetworkInputNormalizer(Context& ctx, const std::string& path) {
  // Define the path for the normalisations
  string DNN_path = path;
  string MeansName = DNN_path+"/means.txt";
  string StdsName = DNN_path+"/stds.txt";

  // get vector vector
  DNNInputs_mean = get_DNN_vals(MeansName);
  DNNInputs_std = get_DNN_vals(StdsName);
}

std::vector<double> NeuralNetworkInputNormalizer::normalizeInputs(uhh2::Event& event, std::vector<double> DNNInputs) {

  // the means and stds fpr AddInputs are expected to be contained in the same files as for the "regular" inputs
  // if the DNNInputs-Vector is not correctly matching this, the code will crash here.
  // std::cout << "Inputs: " << DNNInputs.size() << " Means: " << DNNInputs_mean.size() << " Stds: " << DNNInputs_std.size() << std::endl;
  assert(DNNInputs.size() == DNNInputs_mean.size());
  assert(DNNInputs.size() == DNNInputs_std.size());


  for(uint i = 0; i < DNNInputs.size(); i++){
    DNNInputs.at(i) = (DNNInputs.at(i)-DNNInputs_mean.at(i))/DNNInputs_std.at(i);
  }

  return DNNInputs;

}

// ######################################
// ######################################
// ######################################

NeuralNetworkIncluder::NeuralNetworkIncluder(Context& ctx, bool parametrized) {
  is_parametrized = parametrized;
  //path = "/nfs/dust/cms/user/flabe/MLCorner/TstarNN/reweightingApproach/NonParametric/";
  // the old one: path = "/nfs/dust/cms/user/flabe/TstarTstar/Jupyter/TstarTstar_NN/models/chosenModel/main__layers_25_25_25_25_1__activation_tanh__loss_mean_squared_error__lr_1e-05__batchsize_4096__dropout_0__batchnorm_0/frozen/";
  path = "/nfs/dust/cms/user/flabe/TstarTstar/Jupyter/TstarTstar_NN/models_retraining_2024/UL16postVFP__layers_25_25_25_25_1__activation_tanh__loss_mean_squared_error__lr_1e-05__batchsize_4096__dropout_0__batchnorm_0/frozen/";
  NNInputCreator.reset(new NeuralNetworkInputCreator(ctx));
  //NNInputNormalizer.reset(new NeuralNetworkInputNormalizer(ctx, path+"/data/"));
  NNInputNormalizer.reset(new NeuralNetworkInputNormalizer(ctx, path));
  //NNModule.reset(new NeuralNetworkModule(ctx, path+"/model.pb", path+"/model.config.pbtxt"));
  NNModule.reset(new NeuralNetworkModule(ctx, path+"/frozen_graph.pb", path+"/frozen_graph.config.pbtxt"));
  h_masspoint = ctx.get_handle<double>("masspoint");
  h_DNN_output = ctx.declare_event_output<double>("DNN_output");
  //h_DoAddInputs = ctx.get_handle<bool>("doAddInputs");
}

bool NeuralNetworkIncluder::process(Event& event) {
  NNInputCreator->createInputs(event);
  std::vector<double> inputs = NNInputCreator->getInputs();
  /**
  if(event.get(h_DoAddInputs)) {
    NNInputCreator->createAddInputs(event);
    std::vector<double> Addinputs = NNInputCreator->getAddInputs();
    inputs.insert(inputs.end(), Addinputs.begin(), Addinputs.end());
  }
  **/
  if(is_parametrized) inputs.push_back(event.get(h_masspoint));

  // Normalizing
  inputs = NNInputNormalizer->normalizeInputs(event, inputs);

  NNInputCreator->setInputs(inputs);
  NNModule->process(event, NNInputCreator->getInputs());
  std::vector<tensorflow::Tensor> NNOutputs = NNModule->GetOutputs();
  event.set(h_DNN_output, (double)(NNOutputs.at(0).tensor<float, 2>()(0,0)));
  return true;
}

// ######################################
// ######################################
// ######################################

NeuralNetworkInputWriter::NeuralNetworkInputWriter(Context& ctx) {
  NNInputCreator.reset(new NeuralNetworkInputCreator(ctx));

  for(uint i = 0; i < 33; i++){
    handles.push_back(ctx.declare_event_output<double>("DNN_Input_"+std::to_string(i)));
  }

  for(uint i = 0; i < 3; i++){
    Addhandles.push_back(ctx.declare_event_output<double>("DNN_AddInput_"+std::to_string(i)));
  }

  // Define Handle for DNN
  h_DNN_Inputs = ctx.declare_event_output<std::vector<double>>("DNN_Inputs");
  h_DNN_AddInputs = ctx.declare_event_output<std::vector<double>>("DNN_AddInputs");
  h_masspoint = ctx.declare_event_output<double>("masspoint");

  is_MC = ctx.get("dataset_type") == "MC";
}

bool NeuralNetworkInputWriter::process(uhh2::Event& event) {

  NNInputCreator->createInputs(event);
  NNInputCreator->createAddInputs(event);
  std::vector<double> Inputs = NNInputCreator->getInputs();
  std::vector<double> AddInputs = NNInputCreator->getAddInputs();
  event.set(h_DNN_Inputs, Inputs);
  event.set(h_DNN_AddInputs, AddInputs);

  int i = 0;
  for(auto & handle : handles){
    event.set(handle, Inputs.at(i));
    i++;
  }

  i = 0;
  for(auto & handle : Addhandles){
    event.set(handle, AddInputs.at(i));
    i++;
  }

  double masspoint = -1;
  if(is_MC){
    assert(event.genparticles);
    for(const GenParticle & gp : *event.genparticles){
      if((gp.pdgId() == 9000005 || gp.pdgId() == -9000005) && (gp.status()==23 || gp.status()==22)){
        masspoint = inv_mass(gp.v4());
      }
    }
  }
  if(masspoint == -1){ // if no TstarFound (eg Data or Background) set random between 200 and 2000, inclusively!
    masspoint = std::rand() % 1801 + 200;
  }
  event.set(h_masspoint, masspoint);


  return true;
}
