// Some hist macro
// author: F.Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q costhetaPlots.C

void costhetaPlots(bool reweighted = false){
  
  Double_t w = 800;
  Double_t h = 600;

  std::vector<TString> histnames={"rest_costheta_mother_gluon","rest_costheta_lowpt_mother_gluon", "rest_costheta_medpt_mother_gluon", "rest_costheta_highpt_mother_gluon", "CS_costheta_mother_gluon","CS_costheta_lowpt_mother_gluon", "CS_costheta_medpt_mother_gluon", "CS_costheta_highpt_mother_gluon", "rest_costheta_top_gluon","rest_costheta_lowpt_top_gluon", "rest_costheta_medpt_top_gluon", "rest_costheta_highpt_top_gluon", "CS_costheta_top_gluon","CS_costheta_lowpt_top_gluon", "CS_costheta_medpt_top_gluon", "CS_costheta_highpt_top_gluon"};

  TString path = "/nfs/dust/cms/user/flabe/CMSSW/TstarTstar/102X_v1/Analysis/2016/hadded/";
  TString filename_base = "uhh2.AnalysisModuleRunner.MC.";
  std::vector<TString> filenames_signal = {"TstarTstar_M-700", "TstarTstar_M-1600"};
  std::vector<TString> filenames_bkg = {"TTbar"};
  TString subpath="Top_check";
  if(reweighted) subpath+="_reweighted";

  for(const auto & histname : histnames){
    TCanvas *canvas = new TCanvas("chist", "c", w, h);
    TLegend *leg = new TLegend(0.2,0.2);

    canvas->SetLogy();

    // draw background stack
    THStack *BG_stack = new THStack("BG","");
    for(const auto filename : filenames_bkg){
      TFile *input = TFile::Open(path+filename_base+filename+".root");
      TH1D *hist = (TH1D*)input->Get(subpath+"/"+histname); 
      if(!hist) cout<<"Hist is empty"<<endl;;
      hist->SetFillColor(810);
      hist->Scale(1/(filenames_bkg.size()*hist->Integral()));
      leg->AddEntry(hist, filename, "f");
      BG_stack->Add(hist);
    }
    BG_stack->Draw("hist");
    BG_stack->SetTitle("");
    BG_stack->SetMaximum(0.2);
    BG_stack->SetMinimum(7e-3);
    BG_stack->GetXaxis()->SetTitle("cos(\theta)");
    BG_stack->GetYaxis()->SetTitle("events / total events");
    canvas->Update();

    // draw signals
    int i = 2;
    for(const auto & filename : filenames_signal){
      TFile *input = TFile::Open(path+filename_base+filename+".root");
      TH1D *hist = (TH1D*)input->Get(subpath+"/"+histname); 
      if(!hist) cout<<"Hist is empty"<<endl;;
      hist->Scale(1/hist->Integral());
      hist->GetXaxis()->SetTitle("cos(\theta)");
      hist->GetYaxis()->SetTitle("events / total events");
      hist->SetTitle(histname);
      hist->SetLineColor(1);
      hist->SetLineStyle(i);
      leg->AddEntry(hist, filename, "l");
      hist->Draw("hist same");
      i++;
    }


    leg->Draw("same");
    canvas->SaveAs("plot_"+subpath+"_"+histname+".pdf");
  }

}

