

void JECJERcheck(){

  // some style options
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
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.05);

  gROOT->ForceStyle();
  Double_t w = 800;
  Double_t h = 600;

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  c1_hist->SetLogz();

  // get histograms
  TString filename = "TTbarToSemiLeptonic_UL18";
  TString subpath = "Aftertopptreweighting";
  TString histname = "pt_ST_rebinned";
  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/UL18/";

  std::vector<TString> folders = {"nominal", "JEC_up", "JEC_down", "JER_up", "JER_down"};
  std::vector<int> colors = {1, 2, 2, 4, 4};
  std::vector<int> styles = {1, 2, 3, 2, 3};

  TString fileprefix = "uhh2.AnalysisModuleRunner.MC.";


  TCanvas *canvas = new TCanvas("can_crosscheck", "c", 500, 500);
  TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.35,1.0,1.0);
  TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.35);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  pad1->SetLogy();

  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);

  auto legend = new TLegend(0.55,0.7,0.78,0.88);
  gStyle->SetLegendTextSize(0.05);

  int i = 1;
  for (auto folder : folders) {
    TFile *input = TFile::Open(path+"/"+folder+"/"+fileprefix+filename+".root");
    if(!input) cout << "Empty file" << endl;
    TH1D *hist = (TH1D*)input->Get(subpath+"/"+histname);
    if(!hist) cout << "Empty hist" << endl;

    hist->SetLineColor(colors.at(i-1));
    hist->SetLineStyle(styles.at(i-1));
    if(i == 1) {
      hist->SetTitle("");
      hist->GetYaxis()->SetTitle("events");
      hist->Draw("hist");
    }
    else hist->Draw("hist same");

    legend->AddEntry(hist, folder, "l");

    i++;
  }

  legend->SetBorderSize(0);
  legend->Draw();

  // #########
  pad2->cd();
  // #########

  TFile *input_ref = TFile::Open(path+"/"+folders.at(0)+"/"+fileprefix+filename+".root");
  if(!input_ref) cout << "Empty file" << endl;
  TH1D *hist_ref = (TH1D*)input_ref->Get(subpath+"/"+histname);
  if(!hist_ref) cout << "Empty hist" << endl;

  i = 1;
  for (auto folder : folders) {

    TFile *input = TFile::Open(path+"/"+folder+"/"+fileprefix+filename+".root");
    if(!input) cout << "Empty file" << endl;
    TH1D *hist = (TH1D*)input->Get(subpath+"/"+histname);
    if(!hist) cout << "Empty hist" << endl;

    hist->SetLineColor(colors.at(i-1));
    hist->SetLineStyle(styles.at(i-1));
    hist->Divide(hist_ref);

    if(i == 2) {
      hist->GetXaxis()->SetTitle( hist->GetTitle() );
      hist->GetYaxis()->SetRangeUser(0.8, 1.2);
      hist->SetTitle("");
      hist->Draw("hist");
    }
    else hist->Draw("hist same");


    i++;
  }

  canvas->SaveAs("plots/JECJERcheck.pdf");


}
