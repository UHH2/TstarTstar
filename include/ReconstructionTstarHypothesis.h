#pragma once

#include "UHH2/core/include/Particle.h"
#include "UHH2/core/include/Photon.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/common/include/ReconstructionHypothesis.h"
#include <map>

/**
 *  @short container class to store the results of the Tstar-Tstar (=top quark + gluon/gamma + top quark + gluon) reconstruction
 * 
 * The ttbar (semileptonic) reconstruction stored as ReconstructionHypothesis 
 * ReconstructionTstarHypothesis then consists of  assigning the jets (and/or photon) of the event 
 * to either Tstar with the the leptonically decaying top or Tstar with the hadronically decaying top. In addition
 * to accessing these information (i.e., which jets are assigned to which side, etc.), each
 * hypothesis can have a number of associated *discriminators*. A discriminator is identified
 * by name and is a floating point value meant to measure how 'good' the hypothesis is according to some criterion;
 * see ReconstructionTstarHypothesisDiscriminators.h for different criteria and to fill the discriminators.
 */
class ReconstructionTstarHypothesis {
public:
  explicit ReconstructionTstarHypothesis(){}

  ReconstructionHypothesis ttbar_hyp() const{ return m_ttbar_hyp;}
  LorentzVector tstarlep_v4() const{return m_tstarlep_v4;} // Tstar with top decay leptonically
  LorentzVector tstarhad_v4() const{return m_tstarhad_v4;}  // Tstar with top decay hadronically

  LorentzVector tstar1gamma_v4() const{return m_tstar1gamma_v4;} //1st Tstar -> t + gamma 
  LorentzVector tstar1gluon_v4() const{return m_tstar1gluon_v4;} //1st Tstar -> t + gluon

  LorentzVector tstar2gamma_v4() const{return m_tstar2gamma_v4;} //2nd Tstar -> t + gamma 
  LorentzVector tstar2gluon_v4() const{return m_tstar2gluon_v4;} //2nd Tstar -> t + gluon

  LorentzVector gluon1_v4() const{return m_gluon1_v4;} //1st gluon from Tstar -> t + gluon
  LorentzVector gluon2_v4() const{return m_gluon2_v4;} //2nd gluon from Tstar -> t + gluon

  double chi2() const{return m_chi2;} // return designated chi2 value TODO maybe replace by structure as in ReconstructionHypothesis

  const std::vector<Jet>& tstarlep_jets() const{return m_tstarlep_jets;}
  const std::vector<Jet>& tstarhad_jets() const{return m_tstarhad_jets;}
  const std::vector<Photon>& tstarlep_photons() const{return m_tstarlep_photons;}
  const std::vector<Photon>& tstarhad_photons() const{return m_tstarhad_photons;}


  /* /// get the discriminator value for this hypothesis; thows a runtime_error if it does not exist. */
  /* float discriminator(const std::string & l) const { */
  /*     auto it = m_discriminators.find(l); */
  /*     if(it == m_discriminators.end()){ */
  /*         throw std::runtime_error("ReconstructionTstarHypothesis::discriminator: discriminator with label '" + l + "' not set"); */
  /*     } */
  /*     return it->second; */
  /* } */
  
  /* /// test if a discriminator value with a certian label has already been added */
  /* bool has_discriminator(const std::string & label) const { */
  /*     return m_discriminators.find(label) != m_discriminators.end(); */
  /* } */
  
  void set_tstarlep_v4(LorentzVector v4){m_tstarlep_v4=v4;}
  void set_tstarhad_v4(LorentzVector v4){m_tstarhad_v4=v4;} 
  void set_ttbar(ReconstructionHypothesis ttbar_hyp){m_ttbar_hyp = ttbar_hyp;}

  void set_tstar1gamma_v4(LorentzVector v4){m_tstar1gamma_v4 = v4;}
  void set_tstar2gamma_v4(LorentzVector v4){m_tstar2gamma_v4 = v4;}

  void set_tstar1gluon_v4(LorentzVector v4){m_tstar1gluon_v4 = v4;}
  void set_tstar2gluon_v4(LorentzVector v4){m_tstar2gluon_v4 = v4;}

  void set_gluon1_v4(LorentzVector v4){m_gluon1_v4 = v4;}//FixME: assumes gluon = 1 jet
  void set_gluon2_v4(LorentzVector v4){m_gluon2_v4 = v4;}//FixME: assumes gluon = 1 jet

  void set_used_ttbarjet_flags(const std::vector<bool> flags){m_used_ttbarjet_flags = flags;}
  std::vector<bool> used_ttbarjet_flags(){return m_used_ttbarjet_flags;}
  void add_tstarlep_jet(const Jet& j){m_tstarlep_jets.push_back(j);}
  void add_tstarhad_jet(const Jet& j){m_tstarhad_jets.push_back(j);}

  void add_tstarlep_photon(const Photon& j){m_tstarlep_photons.push_back(j);}
  void add_tstarhad_photon(const Photon& j){m_tstarhad_photons.push_back(j);}

  void set_chi2(double chi2){m_chi2 = chi2;}

  /* void set_discriminator(const std::string & label, float discr){ */
  /*     m_discriminators[label] = discr; */
  /* } */
  
private:

  ReconstructionHypothesis m_ttbar_hyp;

  LorentzVector m_tstarlep_v4;
  LorentzVector m_tstarhad_v4;
  LorentzVector m_tstargluon_v4;
  LorentzVector m_tstargamma_v4;
  LorentzVector m_tstar1gamma_v4, m_tstar2gamma_v4;
  LorentzVector m_tstar1gluon_v4, m_tstar2gluon_v4;
  LorentzVector m_gluon1_v4, m_gluon2_v4;
  double m_chi2;

  std::vector<Jet> m_tstarlep_jets;
  std::vector<Jet> m_tstarhad_jets;

  std::vector<Photon> m_tstarlep_photons;
  std::vector<Photon> m_tstarhad_photons;

  std::vector<bool> m_used_ttbarjet_flags;

  //  std::map<std::string, float> m_discriminators;
};

