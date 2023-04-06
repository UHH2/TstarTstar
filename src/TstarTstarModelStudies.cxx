#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/PhotonIds.h"
#include <UHH2/common/include/MuonIds.h>
#include "UHH2/common/include/MCWeight.h"

#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarModelHists.h"


using namespace std;
using namespace uhh2;

namespace uhh2 {

/** bla
 *
 * blub
 *
 */
class TstarTstarModelStudies: public AnalysisModule {
public:

  explicit TstarTstarModelStudies(Context & ctx);
  virtual bool process(Event & event) override;

private:

  // ##### Histograms #####
  // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_ModelHists, h_ModelHists_ptReweighted;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;

  TH1D* hist_Tstar_pt_12;
  TH1D* hist_Tstar_pt_32;

};


TstarTstarModelStudies::TstarTstarModelStudies(Context & ctx){

  for(auto & kv : ctx.get_all()){
    cout << " " << kv.first << " = " << kv.second << endl;
  }

  ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
  h_ModelHists.reset(new TstarTstarModelHists(ctx, "ModelHists"));
  h_ModelHists_ptReweighted.reset(new TstarTstarModelHists(ctx, "ModelHists_ptreweighted"));

  TFile *f = new TFile("/nfs/dust/cms/user/flabe/TstarTstar/CMSSW_10_2_17/src/UHH2/TstarTstar/macros/ModelStudies_pt_weights.root");
  hist_Tstar_pt_12 = (TH1D*)f->Get("spin12_Tstar_pt");
  hist_Tstar_pt_32 = (TH1D*)f->Get("spin32_Tstar_pt");

}


bool TstarTstarModelStudies::process(Event & event) {

  ttgenprod->process(event);
  h_ModelHists->fill(event);

  // trying pt reweighted scheme
  // finding Tstars and determining which model is used
  std::vector<GenParticle> Tstars;
  bool isSpin12 = false;
  for(const GenParticle & gp : *event.genparticles){
    if((gp.pdgId() == 600) || (gp.pdgId() == -600)) { // 600 or 25001
      isSpin12 = true;
      Tstars.push_back(gp);
    }
    else if((gp.pdgId() == 9000005 && (gp.status()==23 || gp.status()==22)) || (gp.pdgId() == -9000005 && (gp.status()==23 || gp.status()==22))) {
      if(isSpin12) std::cout << "Error, models mixed?" << std::endl;
      Tstars.push_back(gp);
    }
  }
  double Tstarpt = 0;
  for (const auto & Tstar : Tstars) {
    Tstarpt+=Tstar.pt();
  }
  Tstarpt/=Tstars.size();

  double reweight = 1;
  if(isSpin12) reweight = hist_Tstar_pt_12->GetBinContent(hist_Tstar_pt_12->GetXaxis()->FindBin(Tstarpt));
  else reweight = hist_Tstar_pt_32->GetBinContent(hist_Tstar_pt_32->GetXaxis()->FindBin(Tstarpt));
  if(reweight == 0) {
    std::cout << "Weight was 0 for pt = " << Tstarpt <<", set to 1!" << std::endl;
    reweight = 1;
  }

  double weight = event.weight;
  event.weight = event.weight/reweight;
  h_ModelHists_ptReweighted->fill(event);
  event.weight = weight;

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarPreselectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarModelStudies)

}
