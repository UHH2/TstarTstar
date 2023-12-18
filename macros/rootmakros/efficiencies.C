// Efficiency histograms for trigger related studies
// author: F.Labe
// date: 24.09.2019
// Run it with following command:
// root -l -q efficiencies.C

void efficiencies(TString suffix = ""){

  gROOT->SetBatch(kTRUE); // to not open canvas and get XQuartz in the way

  Double_t w = 450;
  Double_t h = 400;

  // Reuseable buffers
  TFile *input;
  TH1D *hist;

  // Drawing Definitions
  TCanvas *canvas = new TCanvas("chist", "c", w, h);
  TLegend *leg = new TLegend(0.22,0.2,0.4,0.60);
  //TLegend *leg = new TLegend(0.45,0.755,0.825,0.88);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);

  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.16);  gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.09);
  gPad->SetTicky();

  canvas->SetLogy();
  //leg->SetNColumns(4);

  // Defining paths
  TString pathPresel = "/nfs/dust/cms/user/flabe/TstarTstar/data/Preselection/UL18/hadded/";
  TString pathSel = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/UL18/hadded/";
  TString pathReco = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/UL18/hadded/";
  TString fileprefix = "uhh2.AnalysisModuleRunner.MC.";
  TString histname = "N_jets";

  // Defining Steps
  std::vector<TString> preselSteps = {"AfterLep",  "AfterJets", "AfterST"};
  std::vector<TString> selSteps = {"METcut", "AK4cut", "HOTVRcut", "btagcut", "2Dcut", "STcut"};
  std::vector<TString> recoSteps = {};
  int stepcount = preselSteps.size() + selSteps.size() + recoSteps.size();
  std::vector<TString> labels = {"Trigger and N_{lep} = 1", "jet pre-cut", "ST pre-cut", "MET cut", "AK4 cut", "HOTVR cut", "b-tag cut", "2D cut", "ST cut", "not visible"};

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
  std::vector<TString> signalSamples = {"TstarTstar_M-800", "TstarTstar_M-1800"};
  std::vector<TString> BGSamples = {"TTbar", "WJets", "ST", "QCD"};
  //std::vector<TString> BGSamples = {};

  std::vector<TString> signal_labels = {"t* 800 GeV",  "t* 1800 GeV", };
  std::vector<TString> BG_labels = {"TTbar", "WJets", "ST", "QCD"};

  // Defining Drawing options
  std::vector<int> colors_Signal = {1,  1,};
  std::vector<int> line_Signal = {1, 2};
  std::vector<int> colors_BG = {810, 600,800,  867};

  // ########################
  // ## Finish Definitions ##
  // ########################

  std::cout << "initials" << std::endl;

  // Saving initial values
  std::vector<double> initial_signal;
  std::vector<double> initial_BG;
  double initial_BG_sum = 0;
  for(const auto & sample : signalSamples){
    input = TFile::Open(pathPresel+fileprefix+sample+".root");
    hist = (TH1D*)input->Get(preselSteps.at(0)+"/"+histname);
    initial_signal.push_back(hist->Integral());
  }
  for(const auto & sample : BGSamples){
    input = TFile::Open(pathPresel+fileprefix+sample+".root");
    hist = (TH1D*)input->Get(preselSteps.at(0)+"/"+histname);
    double val = hist->Integral();
    initial_BG.push_back(val);
    initial_BG_sum += val;
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
      input = TFile::Open(pathPresel+fileprefix+sample+".root");
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_Signal.at(index_sample)->SetBinContent(index_step, val/initial_signal.at(index_sample));
      index_sample++;
    }
    index_sample = 0;
    for(const auto & sample : BGSamples){
      input = TFile::Open(pathPresel+fileprefix+sample+".root");
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_BG.at(index_sample)->SetBinContent(index_step, val/initial_BG.at(index_sample));
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
      input = TFile::Open(pathSel+fileprefix+sample+".root");
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_Signal.at(index_sample)->SetBinContent(index_step, val/initial_signal.at(index_sample));
      index_sample++;
    }
    index_sample = 0;
    for(const auto & sample : BGSamples){
      input = TFile::Open(pathSel+fileprefix+sample+".root");
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_BG.at(index_sample)->SetBinContent(index_step, val/initial_BG.at(index_sample));
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
      input = TFile::Open(pathReco+fileprefix+sample+".root");
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_Signal.at(index_sample)->SetBinContent(index_step, val/initial_signal.at(index_sample));
      index_sample++;
    }
    index_sample = 0;
    for(const auto & sample : BGSamples){
      input = TFile::Open(pathReco+fileprefix+sample+".root");
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      cutflow_BG.at(index_sample)->SetBinContent(index_step, val/initial_BG.at(index_sample));
      index_sample++;
    }
    index_step++;
  }

  std::cout << "Done steps" << std::endl;

  // Output QCD efficiency
  /**
  int indexQCD = 1;
  double val_pre = cutflow_Signal.at(indexQCD)->GetBinContent(cutflow_Signal.at(indexQCD)->FindLastBinAbove()-1);
  double val_post = cutflow_Signal.at(indexQCD)->GetBinContent(cutflow_Signal.at(indexQCD)->FindLastBinAbove());
  std::cout << "Reconstruction efficiency for " << BG_labels.at(indexQCD) << ": " << val_post/val_pre << std::endl;
  **/

  // ##########################
  // ## Finish Filling Hists ##
  // ##########################

  std::cout << "drawing BG" << endl;

  // drawing BG
  int index_draw_BG = 0;
  for(const auto & hist : cutflow_BG){
    int index_label = 1;
    for(const auto & label : labels) {hist->GetXaxis()->SetBinLabel(index_label, label); index_label++;}
    hist->SetLineWidth(4);
    hist->SetTitle("");
    hist->SetLineColor(colors_BG.at(index_draw_BG));
    hist->GetYaxis()->SetTitle("Efficiency");
    if(index_draw_BG == 0){
      hist->SetStats(0);
      hist->GetXaxis()->SetTitleFont(42);
      hist->GetXaxis()->SetLabelFont(42);
      hist->GetYaxis()->SetTitleFont(42);
      hist->GetYaxis()->SetLabelFont(42);

      hist->GetXaxis()->SetLabelSize(0.045); // 0.055
      hist->GetXaxis()->SetLabelOffset(0.01);
      hist->GetXaxis()->SetTickLength(0.03);
      //hist->GetXaxis()->SetTitleSize(0.05);
      //hist->GetXaxis()->SetTitleOffset(1.2);

      hist->GetYaxis()->SetTitleOffset(1.5); // 1.8
      hist->GetYaxis()->SetTitleSize(0.06); // 0.05
      hist->GetYaxis()->SetLabelSize(0.05); // 0.045
      hist->GetYaxis()->SetTickLength(0.02);
      hist->GetYaxis()->SetLabelOffset(0.011);

      hist->SetMinimum(1e-4);
      hist->Draw();
    }
    else hist->Draw("same");
    leg->AddEntry(hist, BG_labels.at(index_draw_BG), "l");
    index_draw_BG++;
  }

  std::cout << "Drawing signal" << endl;

  // drawing signal
  int index_draw_signal = 0;

  for(const auto & hist : cutflow_Signal){
    int index_label = 1;
    for(const auto & label : labels) {hist->GetXaxis()->SetBinLabel(index_label, label); index_label++;}
    hist->SetLineWidth(4);
    hist->SetTitle("");
    hist->SetLineColor(colors_Signal.at(index_draw_signal)); // This can be improved by use of iterator
    hist->SetLineStyle(line_Signal.at(index_draw_signal)); // This can be improved by use of iterator
    hist->GetYaxis()->SetTitle("Efficiency");
    hist->Draw("same");
    leg->AddEntry(hist, signal_labels.at(index_draw_signal), "l");
    index_draw_signal++;
  }

  // Draw legend
  leg->SetBorderSize(0);
  leg->Draw("same");

  // draw Lumi text
  TString infotext = TString::Format("%3.1f fb^{-1} (%d TeV)",  59.83, 13);
  TLatex *text = new TLatex(3.5, 24, infotext);
  text->SetNDC();
  text->SetTextAlign(33);
  text->SetX(0.93);
  text->SetTextFont(42);
  text->SetY(1);
  text->SetTextSize(0.045);
  text->Draw();

  // draw CMS Work in Progress text
  TString cmstext = "CMS";
  TLatex *text2 = new TLatex(3.5, 24, cmstext);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.5);
  text2->SetTextFont(62);
  text2->SetTextSize(0.05);
  text2->SetY(0.3);
  text2->Draw();
  TString preltext = "Work in Progress";
  TLatex *text3 = new TLatex(3.5, 24, preltext);
  text3->SetNDC();
  text3->SetTextAlign(13);
  text3->SetX(0.5);
  text3->SetTextFont(52);
  text3->SetTextSize(0.035);
  text3->SetY(0.25);
  text3->Draw();

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
  canvas->SaveAs("plots/efficiency"+suffix+".pdf");

  // signal effs
  int index_sample = 0;
  for(const auto & sample : signalSamples){
    double pre_ttag = cutflow_Signal.at(index_sample)->GetBinContent(8);
    double post_ttag = cutflow_Signal.at(index_sample)->GetBinContent(9);
    cout << "For " << sample << " the efficiency is " << post_ttag/pre_ttag << endl;
    index_sample++;
  }

  // background effs
  index_sample = 0;
  for(const auto & sample : BGSamples){
    double pre_ttag = cutflow_BG.at(index_sample)->GetBinContent(8);
    double post_ttag = cutflow_BG.at(index_sample)->GetBinContent(9);
    cout << "For " << sample << " the efficiency is " << post_ttag/pre_ttag << endl;
    index_sample++;
  }

}
