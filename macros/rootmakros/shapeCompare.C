// Just a basic macro to create simple plots from given hists.
// author: F. Labe
// Run it with following command:
// root -l -b -q shapeCompare.C

void shapeCompare(){

  gStyle->SetOptStat(0);

  Double_t w = 800;
  Double_t h = 600;

  TString histname = "pt_ST_rebinned";
  TString sampleType = "MC";
  TString sample = "TTbar";

  std::vector<TString> files = {
    "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/UL18/hadded/uhh2.AnalysisModuleRunner."+sampleType+"."+sample+".root",
    "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/UL18/hadded/uhh2.AnalysisModuleRunner."+sampleType+"."+sample+".root",
  };
  std::vector<TString> folders = {"AfterNNLO", "Aftertopptreweighting"};
  assert(file.size() == folder.size());

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->SetBorderSize(0);

  TPad *pad1 = new TPad("pad1", "The pad 80% of the height", 0.0, 0.35, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "The pad 20% of the height", 0.0, 0.0,  1.0, 0.35);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  pad1->SetLogy();

  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);

  leg->SetTextFont(42);
  leg->SetTextSize(0.035);

  int i = 0;
  for (auto file : files) {
    cout << "Processing " << file << "." << endl;
    TString folder = folders.at(i);
    TFile *input = TFile::Open(file);
    if(!input) cout<<"File is empty"<<endl;
    TH1D *hist = (TH1D*)input->Get(folder+"/"+histname); //histogram
    if(!hist) cout<<"Hist is empty"<<endl;
    //hist->Scale(1/hist->Integral());
    hist->GetYaxis()->SetTitle("events");
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(i+1);
    hist->SetLineColor(i+1);
    if(i == 0) hist->Draw("hist");
    else hist->Draw("hist same");
    leg->AddEntry(hist, folder, "l");
    i++;
  }

  leg->Draw("same");

  pad2->cd();

  TString base_filename = files.at(0);
  TString base_folder = folders.at(0);
  TFile *base_file = TFile::Open(base_filename);
  TH1D *base_hist = (TH1D*)base_file->Get(base_folder+"/"+histname);
  //base_hist->Scale(1/base_hist->Integral());

  i = 0;
  for (auto file : files) {
    if (i > 0) {
      cout << "Processing " << file << "." << endl;
      TString folder = folders.at(i);
      TFile *input = TFile::Open(file);
      if(!input) cout<<"File is empty"<<endl;
      TH1D *hist = (TH1D*)input->Get(folder+"/"+histname); //histogram
      if(!hist) cout<<"Hist is empty"<<endl;
      //hist->Scale(1/hist->Integral());
      hist->Divide(base_hist);
      hist->GetYaxis()->SetTitle("events");
      hist->GetXaxis()->SetTitle(hist->GetTitle());
      hist->SetTitle("");
      hist->SetMarkerStyle(20);
      hist->SetMarkerColor(i+1);
      hist->SetLineColor(i+1);
      hist->GetYaxis()->SetRangeUser(0.5, 1.5);
      if(i == 1) hist->Draw("hist");
      else hist->Draw("hist same");
    }

    i++;
  }

  c1_hist->SaveAs("plots/shape_"+histname+"_"+sample+".pdf");
}
