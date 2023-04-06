#pragma once

#include "UHH2/common/include/NeuralNetworkBase.hpp"
#include <string>

class NeuralNetworkModule : public NeuralNetworkBase {
public:
  explicit NeuralNetworkModule(uhh2::Context&, const std::string & ModelName, const std::string& ConfigName);
  virtual bool process(uhh2::Event&, std::vector<double> values);
  virtual void CreateInputs(uhh2::Event & event, std::vector<double> values);
};

class NeuralNetworkInputCreator: uhh2::AnalysisModule {
public:
  explicit NeuralNetworkInputCreator(uhh2::Context&);
  virtual bool process(uhh2::Event&) override {return false;}
  virtual bool createInputs(uhh2::Event&);
  virtual bool createAddInputs(uhh2::Event&);
  virtual std::vector<double> getInputs();
  virtual std::vector<double> getAddInputs();
  virtual bool setInputs(std::vector<double>);
  virtual bool setAddInputs(std::vector<double>);
private:
  uhh2::Event::Handle<bool> h_do_masspoint;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
  uhh2::Event::Handle<bool> h_is_muevt;
  std::string DNN_path;
  bool is_MC;
  std::vector<double> DNNInputs;
  std::vector<double> AddDNNInputs;
};

class NeuralNetworkInputNormalizer: uhh2::AnalysisModule {
public:
  explicit NeuralNetworkInputNormalizer(uhh2::Context&, const std::string&);
  virtual bool process(uhh2::Event&) override {return false;}
  virtual std::vector<double> normalizeInputs(uhh2::Event&, std::vector<double>);
private:
  std::vector<double> DNNInputs_mean;
  std::vector<double> DNNInputs_std;
};

class NeuralNetworkIncluder: uhh2::AnalysisModule {
public:
  explicit NeuralNetworkIncluder(uhh2::Context&, bool parametrized);
  virtual bool process(uhh2::Event&) override;
private:
  bool is_parametrized;
  std::string DNN_path;
  std::unique_ptr<NeuralNetworkModule> NNModule;
  std::unique_ptr<NeuralNetworkInputCreator> NNInputCreator;
  std::unique_ptr<NeuralNetworkInputNormalizer> NNInputNormalizer;
  std::string path;
  uhh2::Event::Handle<double> h_DNN_output;
  uhh2::Event::Handle<double> h_masspoint;
  uhh2::Event::Handle<bool> h_DoAddInputs;
};


class NeuralNetworkInputWriter: uhh2::AnalysisModule {
public:
  explicit NeuralNetworkInputWriter(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
private:
  bool is_MC;
  std::unique_ptr<NeuralNetworkInputCreator> NNInputCreator;

  uhh2::Event::Handle<double> h_masspoint;
  uhh2::Event::Handle<std::vector<double>> h_DNN_Inputs;
  uhh2::Event::Handle<std::vector<double>> h_DNN_AddInputs;

  std::vector<uhh2::Event::Handle<double>> handles;
  std::vector<uhh2::Event::Handle<double>> Addhandles;


};
