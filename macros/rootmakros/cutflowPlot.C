// Efficiency histograms for trigger related studies
// author: F.Labe
// date: 24.09.2019
// Run it with following command:
// root -l -b -q cutflowPlots.C

void cutflowPlot(TString suffix = ""){

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

  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.17);

  gROOT->ForceStyle();
  Double_t w = 800;
  Double_t h = 600;

  // Reuseable buffers
  TFile *input;
  TH1D *hist;

  // Drawing Definitions
  TCanvas *canvas = new TCanvas("chist", "c", w, h);
  //TLegend *leg = new TLegend(0.22,0.175,0.5,0.33);
  TLegend *leg = new TLegend(0.53,0.755,0.81,0.88);
  canvas->SetLogy();
  leg->SetNColumns(2);

  // Defining paths
  TString pathPresel = "/nfs/dust/cms/user/flabe/TstarTstar/data/Preselection/UL16postVFP/hadded/";
  TString pathSel = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/UL16postVFP/hadded/";
  TString pathReco = "/nfs/dust/cms/user/flabe/TstarTstar/data/Analysis/UL16postVFP/hadded/";
  TString fileprefix = "uhh2.AnalysisModuleRunner.MC.";
  TString histname = "N_jets";

  // Defining Steps
  std::vector<TString> preselSteps = {"AfterCommon", "AfterLep", "AfterJets", "AfterST"};
  std::vector<TString> selSteps = {"AfterBtag", "After2D", "AfterdR"};
  std::vector<TString> recoSteps = {};
  int stepcount = preselSteps.size() + selSteps.size() + recoSteps.size();

  // Adding suffix
  int is_first = true;
  for (auto & step : preselSteps){
    if(is_first) is_first = false;
    else step.Append(suffix);
  }
  for (auto & step : selSteps){
    step.Append(suffix);
  }
  for (auto & step : recoSteps){
    step.Append(suffix);
  }


  // Defining Samples
  std::vector<TString> signalSamples = {"TstarTstar_M-700.root", "TstarTstar_M-1000.root", "TstarTstar_M-1600.root"};
  std::vector<TString> BGSamples = {"TTbar.root", "WJets.root", "ST.root", "VV.root"};
  //std::vector<TString> BGSamples = {};

  // Defining Drawing options
  std::vector<int> colors_Signal = {632, 820, 432};
  std::vector<int> colors_BG = {810, 800, 600, 416};
  //std::vector<TString> labels = {"No cuts", "N lepton","N HOTVR jet = 1", "MET", "Trigger", "2D", "N HOTVR jet = 3", "non-overlap AK4"};
  std::vector<TString> labels = {"Trigger", "N_{lep}","N_{jet}", "MET","b-tag", "2D","dR", "should not be visible"};

  // ########################
  // ## Finish Definitions ##
  // ########################

  // Saving initial values
  std::vector<double> initial_signal;
  std::vector<double> initial_BG;
  double initial_BG_sum = 0;
  for(const auto & sample : signalSamples){
    input = TFile::Open(pathPresel+fileprefix+sample);
    hist = (TH1D*)input->Get(preselSteps.at(0)+"/"+histname);
    initial_signal.push_back(hist->Integral());
  }
  for(const auto & sample : BGSamples){
    input = TFile::Open(pathPresel+fileprefix+sample);
    hist = (TH1D*)input->Get(preselSteps.at(0)+"/"+histname);
    initial_BG.push_back(hist->Integral());
    initial_BG_sum += hist->Integral();
  }
  // if no presel is present:
  /**
  for(const auto & sample : signalSamples){
    input = TFile::Open(pathSel+fileprefix+sample);
    hist = (TH1D*)input->Get(selSteps.at(0)+"/"+histname);
    initial_signal.push_back(hist->Integral());
  }
  for(const auto & sample : BGSamples){
    input = TFile::Open(pathSel+fileprefix+sample);
    hist = (TH1D*)input->Get(selSteps.at(0)+"/"+histname);
    initial_BG.push_back(hist->Integral());
    initial_BG_sum += hist->Integral();
  }
  **/

  // Defining Hists
  std::vector<TH1*> cutflow_Signal;
  std::vector<TH1*> cutflow_BG;
  for(const auto & sample : signalSamples){cutflow_Signal.push_back(new TH1D(sample, sample, stepcount, 0, stepcount));}    // TODO remove ".root" in the end
  for(const auto & sample : BGSamples){cutflow_BG.push_back(new TH1D(sample, sample, stepcount, 0, stepcount));}

  // Filling hists
  int index_step = 1;

  // Preselection
  std::cout << "Start preselSteps." << endl;
  for(const auto & step : preselSteps){
    std::cout << "Start preselStep: " << step << endl;
    int index_sample = 0;
    for(const auto & sample : signalSamples){
      input = TFile::Open(pathPresel+fileprefix+sample);
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_Signal.at(index_sample)->SetBinContent(index_step, val/initial_signal.at(index_sample));
      index_sample++;
    }
    index_sample = 0;
    for(const auto & sample : BGSamples){
      input = TFile::Open(pathPresel+fileprefix+sample);
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_BG.at(index_sample)->SetBinContent(index_step, val/initial_BG_sum);
      index_sample++;
    }
    index_step++;
  }

  // Selection
  std::cout << "Start selSteps." << endl;
  for(const auto & step : selSteps){
    std::cout << "Start selStep: " << step << endl;
    int index_sample = 0;
    for(const auto & sample : signalSamples){
      input = TFile::Open(pathSel+fileprefix+sample);
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_Signal.at(index_sample)->SetBinContent(index_step, val/initial_signal.at(index_sample));
      index_sample++;
    }
    index_sample = 0;
    for(const auto & sample : BGSamples){
      input = TFile::Open(pathSel+fileprefix+sample);
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_BG.at(index_sample)->SetBinContent(index_step, val/initial_BG_sum);
      index_sample++;
    }
    index_step++;
  }

  // Reconstruction
  std::cout << "Start recoSteps." << endl;
  for(const auto & step : recoSteps){
    std::cout << "Start recoStep: " << step << endl;
    int index_sample = 0;
    for(const auto & sample : signalSamples){
      input = TFile::Open(pathReco+fileprefix+sample);
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_Signal.at(index_sample)->SetBinContent(index_step, val/initial_signal.at(index_sample));
      index_sample++;
    }
    index_sample = 0;
    for(const auto & sample : BGSamples){
      input = TFile::Open(pathReco+fileprefix+sample);
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_BG.at(index_sample)->SetBinContent(index_step, val/initial_BG_sum);
      index_sample++;
    }
    index_step++;
  }

  // ##########################
  // ## Finish Filling Hists ##
  // ##########################

  std::cout << "Start drawing." << endl;
  // drawing background
  if(cutflow_BG.size() == 4){
    THStack *BG_stack = new THStack("BG","");
    std::vector<TString> BG_labels = {"TTbar", "ST", "WJets", "VV"};
    for(int index_draw_BG = cutflow_BG.size()-1; index_draw_BG >= 0; index_draw_BG--){
      int index_label = 1;
      for(const auto & label : labels) {cutflow_BG.at(index_draw_BG)->GetXaxis()->SetBinLabel(index_label, label); index_label++;}
      cutflow_BG.at(index_draw_BG)->SetFillColor(colors_BG.at(index_draw_BG));
      cutflow_BG.at(index_draw_BG)->SetLineColor(colors_BG.at(index_draw_BG));
      leg->AddEntry(cutflow_BG.at(index_draw_BG), BG_labels.at(index_draw_BG), "f");
      BG_stack->Add(cutflow_BG.at(index_draw_BG));
    }
    BG_stack->SetMinimum(5e-4);
    BG_stack->Draw("");
    BG_stack->SetTitle("");
    BG_stack->GetYaxis()->SetTitleSize(0.05);
    BG_stack->GetYaxis()->SetTitle("Efficiency");
    canvas->Update();
  }

  // drawing signal
  int index_draw_signal = 0;
  std::vector<TString> signal_labels = {"Tstar M-700", "Tstar M-1000", "Tstar M-1600"};
  for(const auto & hist : cutflow_Signal){
    int index_label = 1;
    for(const auto & label : labels) {hist->GetXaxis()->SetBinLabel(index_label, label); index_label++;}
    hist->SetLineWidth(3);
    hist->SetTitle("");
    hist->SetLineColor(colors_Signal.at(index_draw_signal)); // This can be improved by use of iterator
    hist->GetYaxis()->SetTitle("Efficiency");
    hist->Draw("same");
    leg->AddEntry(hist, signal_labels.at(index_draw_signal), "l");
    index_draw_signal++;
  }

  // Draw legend
  leg->SetBorderSize(0);
  leg->Draw("same");

  // draw Lumi text
  TString infotext = TString::Format("35.9 fb^{-1} (13 TeV)");
  TLatex *text = new TLatex(3.5, 24, infotext);
  text->SetNDC();
  text->SetTextAlign(33);
  text->SetX(0.83);
  text->SetTextFont(42);
  text->SetY(0.93);
  text->SetTextSize(0.025);
  text->Draw();

  // Draw line after certain bin to split presel, sel,
  bool doLine = false;
  if(doLine){
    double pos1 = cutflow_BG.at(0)->GetBinLowEdge(1+preselSteps.size());
    TLine *line1 = new TLine(pos1,1.25,pos1,1.9);
    //line1->Draw("same");
    double pos2 = cutflow_BG.at(0)->GetBinLowEdge(1+preselSteps.size()+selSteps.size());
    TLine *line2 = new TLine(pos2,1.25,pos2,1.8);
    line2->Draw("same");
  }
  canvas->SaveAs("Cutflow"+suffix+".pdf");

  int index_sample = 0;
  for(const auto & sample : signalSamples){
    double pre_ttag = cutflow_Signal.at(index_sample)->GetBinContent(6);
    double post_ttag = cutflow_Signal.at(index_sample)->GetBinContent(8);
    cout << "For " << sample << " the efficiency is " << post_ttag/pre_ttag << endl;
    index_sample++;
  }

}
