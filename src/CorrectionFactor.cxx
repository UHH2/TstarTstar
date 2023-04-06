#include <UHH2/TstarTstar/include/CorrectionFactor.h>

double CorrectionFactor::get_factor(double pt, double eta){

  // first transform eta to absolute value
  if(eta < 0) eta *= -1;

  // make sure that there is not jet with eta > 2.5
  // (this should not happen because in CombineXCone this should be sorted out)
  if(eta > 2.5) eta = 2.5;

  // get eta bin
  int etabin = 0;
  int N_eta_bins = sizeof(eta_bins)/sizeof(eta_bins[0]) - 1;
  for(int i = 0; i < N_eta_bins; i++){
    if(eta > eta_bins[i]) etabin = i;
  }

  // get factor from function
  if(pt > 430) pt = 430;
  if(pt < 30)  pt = 30;


  double factor = 1;
  if(CorUp || CorDown){
    // first get central factor
    double f_c = CentralCorrectionFunctions[etabin]->Eval(pt);
    // get up/down factor
    double f_ud = UpDownCorrectionGraphs[etabin]->Eval(pt);
    // get difference
    double df = f_ud - f_c;
    // get additional sys (this already is a difference)
    double dg = AdditionalSys->Eval(pt);
    // to get the total factor one has to:
    // 1. add df and dg in quadrature (this gives the total diff to central factor)
    // 2. add this to central factor
    if(CorUp)        factor = f_c + sqrt(df*df + dg*dg);
    else if(CorDown) factor = f_c - sqrt(df*df + dg*dg);

  }
  else factor = CentralCorrectionFunctions[etabin]->Eval(pt);

  return factor;
}


void CorrectionFactor::get_function(const TString & year){

  TString dir = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/CorrectionFile/";
  TString filename;
  filename = dir + "Correction_allHad_"+year+".root";
  TFile *file = new TFile(filename);

  TString histname_updown;
  TString histname_central;
  int N_eta_bins = sizeof(eta_bins)/sizeof(eta_bins[0]) - 1;
  for(int i = 0; i < N_eta_bins; i++){
    stringstream BinString;
    BinString << i;
    histname_updown = BinString.str();
    histname_central = BinString.str();
    if(CorUp){
      histname_updown += "/Up";
      UpDownCorrectionGraphs.push_back((TGraph*)file->Get(histname_updown));
    }
    else if(CorDown){
      histname_updown += "/Down";
      UpDownCorrectionGraphs.push_back((TGraph*)file->Get(histname_updown));
    }
    histname_central += "/Central";
    CentralCorrectionFunctions.push_back((TF1*)file->Get(histname_central));
  }
}

void CorrectionFactor::get_additionalSYS(){

  TString dir = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/CorrectionFile/";
  TString filename;
  filename = dir + "Correction_SysFromResolution_"+year+".root";
  TFile *file = new TFile(filename);

  if(CorUp){
    TString histname = "Up";
    AdditionalSys = (TGraph*)file->Get(histname);
  }
  else if(CorDown){
    TString histname = "Down";
    AdditionalSys = (TGraph*)file->Get(histname);
  }
  return;
}

Particle GetLepton(uhh2::Event & event){

  Particle lepton;
  bool muonislep = false;
  bool elecislep = false;

  if(event.muons->size() > 0 && event.electrons->size() > 0){
    if(event.muons->at(0).pt() > event.electrons->at(0).pt()) muonislep = true;
    else elecislep = true;
  }
  else if(event.muons->size() > 0) muonislep = true;
  else if(event.electrons->size() > 0) elecislep = true;

  if(muonislep) lepton = event.muons->at(0);
  if(elecislep) lepton = event.electrons->at(0);

  return lepton;
}


CorrectionFactor::CorrectionFactor(uhh2::Context & ctx, const std::string & name, const std::string & corvar, bool allHad_, const TString & year_):
h_oldjets(ctx.get_handle<std::vector<TopJet>>("xconeCHS")),
h_newjets(ctx.declare_event_output<std::vector<TopJet>>(name)),
h_cor_factor_had(ctx.declare_event_output<std::vector<float>>("cor_factor_had")),
h_cor_factor_lep(ctx.declare_event_output<std::vector<float>>("cor_factor_lep")),
allHad(allHad_),
year(year_)
{
  CorUp = false;
  CorDown=false;

  if(corvar == "up")   CorUp = true;
  else if(corvar == "down") CorDown = true;

  get_function(year); // read file here, not in every event
  get_additionalSYS();
}

bool CorrectionFactor::process(uhh2::Event & event){
  std::vector<TopJet> oldjets = event.get(h_oldjets);
  std::vector<TopJet> newjets;
  for(unsigned int i=0; i < oldjets.size(); i++) newjets.push_back(oldjets.at(i));

  vector<float> factors_had;
  vector<float> factors_lep;

  // now correct had subjets
  std::vector<Jet> oldsubjets = oldjets.at(0).subjets();
  std::vector<Jet> newsubjets;
  LorentzVector newjet_v4;
  Jet newjet;
  float factor;
  for(unsigned int i=0; i < oldsubjets.size(); i++){
    factor = get_factor(oldsubjets.at(i).pt(), oldsubjets.at(i).eta());
    newjet_v4 = oldsubjets.at(i).v4() * factor;
    newjet = oldsubjets.at(i); // this acutally is not required
    newjet.set_v4(newjet_v4);  // Because of this, all JEC factors are set to the defaul value of 1
    newsubjets.push_back(newjet);
    factors_had.push_back(factor);
  }
  newjets.at(0).set_subjets(newsubjets);

  // also correct second jet
  std::vector<Jet> oldsubjets2 = oldjets.at(1).subjets();
  std::vector<Jet> newsubjets2;
  LorentzVector newjet2_v4;
  Jet newjet2;
  float factor2;
  for(unsigned int i=0; i < oldsubjets2.size(); i++){
    factor2 = get_factor(oldsubjets2.at(i).pt(), oldsubjets2.at(i).eta());
    newjet2_v4 = oldsubjets2.at(i).v4() * factor2;
    newjet2 = oldsubjets2.at(i); // this acutally is not required
    newjet2.set_v4(newjet2_v4);  // Because of this, all JEC factors are set to the defaul value of 1
    newsubjets2.push_back(newjet2);
    factors_lep.push_back(factor2);
  }
  newjets.at(1).set_subjets(newsubjets2);

  event.set(h_cor_factor_had, factors_had);
  event.set(h_cor_factor_lep, factors_lep);
  event.set(h_newjets, newjets);

  return true;
}
