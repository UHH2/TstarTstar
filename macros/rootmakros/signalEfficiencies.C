// Efficiency histograms for trigger related studies
// author: F.Labe
// date: 24.09.2019
// Run it with following command:
// root -l -q efficiencies.C

void signalEfficiencies(TString suffix = ""){

  gROOT->SetBatch(kTRUE); // to not open canvas and get XQuartz in the way

  Double_t w = 420;
  Double_t h = 400;

  // Reuseable buffers
  TFile *input;
  TH1D *hist;

  // Drawing Definitions
  TCanvas *canvas = new TCanvas("chist", "c", w, h);
  TLegend *leg = new TLegend(0.22,0.175,0.4,0.60);
  //TLegend *leg = new TLegend(0.45,0.755,0.825,0.88);

  leg->SetTextFont(42);
  leg->SetTextSize(0.035);

  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.16);  gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.065);

  //canvas->SetLogy();
  //leg->SetNColumns(4);

  // Defining paths
  TString pathPresel = "/nfs/dust/cms/user/flabe/TstarTstar/data/Preselection/hadded/";
  TString pathSel = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/hadded/";
  TString pathReco = "/nfs/dust/cms/user/flabe/TstarTstar/data/Analysis/hadded/";
  TString fileprefix = "uhh2.AnalysisModuleRunner.MC.";
  TString histname = "N_jets";

  // Defining Steps
  std::vector<TString> preselSteps = {"AfterTrigger", "AfterLep", "AfterJets", "AfterFatJets", "AfterMET"};
  std::vector<TString> selSteps = {"AfterBtag", "AfterdR"};
  int stepcount = preselSteps.size() + selSteps.size();

  // Adding suffix
  int is_first = true;
  for (auto & step : preselSteps){
    if(is_first) is_first = false;
    else step.Append(suffix);
  }
  for (auto & step : selSteps){
    step.Append(suffix);
  }


  // Defining Samples
  std::vector<double> signalMasses = {700, 800, 900, 1000, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
  // attention 1100 is missing!
  TString signalNameBase = "TstarTstar_M-";

  // Defining Drawing options
  std::vector<int> colors = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<TString> labels = {"Trigger", "N_{lep}", "N_{AK4}", "N_{HOTVR}", "MET", "b-tag", "kin.", "reco", "should not be visible"};

  // ########################
  // ## Finish Definitions ##
  // ########################

  // Saving initial values
  std::vector<double> initial;
  for(const auto & mass : signalMasses){
    TString sample = signalNameBase+std::to_string((int)mass);
    input = TFile::Open(pathPresel+fileprefix+sample+".root");
    hist = (TH1D*)input->Get(preselSteps.at(0)+"/"+histname);
    initial.push_back(hist->Integral());
  }

  // vector to store graphs
  std::vector<TGraph*> graphs;

  // Filling hists
  int index_step = 1;

  // Preselection
  std::cout << "Start preselSteps." << endl;
  for(const auto & step : preselSteps){
    std::cout << "Start preselStep: " << step << endl;
    std::vector<double> vals;
    int index_sample = 0;
    for(const auto & mass : signalMasses){
      TString sample = signalNameBase+std::to_string((int)mass);
      input = TFile::Open(pathPresel+fileprefix+sample+".root");
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      vals.push_back(val/initial.at(index_sample));
      index_sample++;
    }
    graphs.push_back(new TGraph(signalMasses.size(), &signalMasses[0], &vals[0]));
    index_step++;
  }

  // Selection
  std::cout << "Start selSteps." << endl;
  for(const auto & step : selSteps){
    std::cout << "Start selStep: " << step << endl;
    std::vector<double> vals;
    int index_sample = 0;
    for(const auto & mass : signalMasses){
      auto sample = signalNameBase+std::to_string((int)mass);
      input = TFile::Open(pathSel+fileprefix+sample+".root");
      hist = (TH1D*)input->Get(step+"/"+histname);
      double val = hist->Integral();
      vals.push_back(val/initial.at(index_sample));
      index_sample++;
    }
    graphs.push_back(new TGraph(signalMasses.size(), &signalMasses[0], &vals[0]));
    index_step++;
  }

  // ###########################
  // ## Finish Filling Graphs ##
  // ###########################

  std::cout << "drawing..." << endl;

  // drawing BG
  int index_draw = 0;
  for(const auto & graph : graphs){
    graph->SetLineWidth(4);
    graph->SetTitle("");
    graph->SetLineColor(colors.at(index_draw));
    graph->GetYaxis()->SetTitle("Efficiency");
    graph->GetXaxis()->SetTitle("S_{T} [GeV]");
    graph->SetMarkerStyle(kPlus);
    if(index_draw == 0){
      //graph->SetStats(0);
      graph->GetXaxis()->SetTitleFont(42);
      graph->GetXaxis()->SetLabelFont(42);
      graph->GetYaxis()->SetTitleFont(42);
      graph->GetYaxis()->SetLabelFont(42);

      graph->GetXaxis()->SetLabelSize(0.055); // 0.045
      graph->GetXaxis()->SetLabelOffset(0.01);
      graph->GetXaxis()->SetTickLength(0.03);
      graph->GetXaxis()->SetTitleSize(0.05);
      graph->GetXaxis()->SetTitleOffset(1.2);

      graph->GetYaxis()->SetTitleOffset(1.5); // 1.8
      graph->GetYaxis()->SetTitleSize(0.06); // 0.05
      graph->GetYaxis()->SetLabelSize(0.05); // 0.045
      graph->GetYaxis()->SetTickLength(0.02);
      graph->GetYaxis()->SetLabelOffset(0.011);

      graph->SetMinimum(0);
      graph->SetMaximum(1.01);
      graph->Draw();
      graph->GetXaxis()->SetLimits(700,2000);
      graph->Draw("PC");
      canvas->Update();
    }
    else graph->Draw("same PC");
    leg->AddEntry(graph, labels.at(index_draw), "l");
    index_draw++;
  }

  // Draw legend
  leg->SetBorderSize(0);
  leg->Draw("same");

  // draw Lumi text
  TString infotext = TString::Format("%3.1f fb^{-1} (%d TeV)", 137., 13);
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
  text2->SetY(0.27);
  text2->Draw();
  TString preltext = "Work in Progress";
  TLatex *text3 = new TLatex(3.5, 24, preltext);
  text3->SetNDC();
  text3->SetTextAlign(13);
  text3->SetX(0.5);
  text3->SetTextFont(52);
  text3->SetTextSize(0.035);
  text3->SetY(0.22);
  text3->Draw();

  // saving output
  canvas->SaveAs("plots/signalEfficiency"+suffix+".pdf");

}
