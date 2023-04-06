// Efficiency histograms for trigger related studies
// author: A.Karavdina
// date: 24.09.2019
// Run it with following command:
// root -l -b -q TriggerEffPlots.C

void weightPlots(){
  
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.06,"x");  
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.05,"x");  
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleYOffset(1.20);
  gStyle->SetTitleXOffset(1.0);
  gStyle->SetPalette(55);

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.17);

  gROOT->ForceStyle();
  Double_t w = 600;
  Double_t h = 600;

  //Files after selection
  //We expect histograms filled with and without trigger selection stored in the same file
  TString path = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Analysis/2016/hadded/";
  std::vector<TString> filenames_signal = {"uhh2.AnalysisModuleRunner.MC.TstarTstar_M-700.root", "uhh2.AnalysisModuleRunner.MC.TstarTstar_M-1000.root", "uhh2.AnalysisModuleRunner.MC.TstarTstar_M-1600.root", "uhh2.AnalysisModuleRunner.MC.TstarTstar_M-2000.root"};
  TString filename_background = "uhh2.AnalysisModuleRunner.MC.TTbar.root";
  
  TFile *input_background = TFile::Open(path+filename_background);
  TTree *tree_background = (TTree*)input_background->Get("AnalysisTree");

  TH1D *hist = new TH1D("evt weight", "evt weight", 100, 0, 1000);

  double weight;
  tree_background->SetBranchAddress("evt_weight", &weight);
  uint entries = tree_background->GetEntries();
  for (uint i = 0; i < entries; i++){
    std::cout << "Filling entries... " << i << "/" << entries << " (" << i*100/entries << "%)" << "\r";
    tree_background->GetEntry(i);
    hist->Fill(weight);
  }

  TCanvas *c = new TCanvas("c", "c", w, h);
  hist->Draw();
  c->SaveAs("test.pdf");
  

  // TH1D *hist_denom = (TH1D*)input_denom->Get(subpath+"/"+histname);//histogram before photon selection
  // if(!hist_eff || !hist_denom) cout<<"Hists are empty"<<endl;;
  // if(!hist_eff || !hist_denom) return;
  // if(hist_denom->GetEntries()>0 && hist_eff->GetEntries()>0) hist_eff->Divide(hist_denom);
  // TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  // hist_eff->GetXaxis()->SetTitle(hist_eff->GetTitle());
  // hist_eff->GetYaxis()->SetTitle("Efficiency");
  // hist_eff->GetYaxis()->SetRangeUser(0,1.5);
  // hist_eff->SetTitle("");
  // hist_eff->SetMarkerStyle(20);
  // hist_eff->SetMarkerColor(1);
  // hist_eff->Draw();


  // c1_hist->SaveAs("PhoIDEff_"+label+"_"+photonID[i]+"_"+subpath+"_"+histname+".pdf");

}
