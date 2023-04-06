#include <TString.h>
#include <iostream>
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TText.h>
#include <TPaveText.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TROOT.h>
#include <TKey.h>
#include <TLatex.h>
#include <TClass.h>
#include <fstream>
#include <string>

using namespace std;

void getScaleEnvelope(){

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

  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.07);

  TString filename_base = "";
  TString year = "UL18";
  filename_base += "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/"+year+"/hadded/uhh2.AnalysisModuleRunner.MC.";

  TString region = "SignalRegion";

  vector<TString> samples = {"TTbar", "ST"};
  vector<bool> isSignal (samples.size(), false);
  vector<TString> masspoints = {"800", "1600", "1800", "2500"};
  for (auto mass : masspoints) {
    samples.push_back("TstarTstar_M-" + mass);
    isSignal.push_back(true);
    samples.push_back("TstarTstar_M-" + mass + "_32");
    isSignal.push_back(true);
  }

  TString channel = "mu";

  for(unsigned int i=0; i<samples.size(); i++){

    cout << "sample " << samples.at(i) << endl;

    // check if the sample is a signal or a background
    isSignal.at(i) = (samples.at(i).Index("TstarTstar") == 0 );

    TString filename = filename_base + samples.at(i) + ".root";
    TFile* f_in = new TFile(filename, "READ");

    TH1F *h_nominal = (TH1F*)f_in->Get(region + "_" + channel + "/pt_ST_nominal");
    TH1F *h_scale_up = (TH1F*)h_nominal->Clone();
    TH1F *h_scale_down = (TH1F*)h_nominal->Clone();

    std::vector<TH1F*> variations;
    std::vector<TString> variation_names = {
      "murmuf_upup",
      "murmuf_upnone",
      "murmuf_noneup",
      "murmuf_nonedown",
      "murmuf_downnone",
      "murmuf_downdown"
    };
    for (auto name : variation_names) {
      TString basename = region + "_" + channel + "/pt_ST_";
      variations.push_back((TH1F*)f_in->Get(basename + name));
    }

    // Loop over each bin of the ST histogram and take envelope
    for(int j=1; j < h_nominal->GetXaxis()->GetNbins()+1; j++){

      float highest = 0;
      float lowest = 9999999;

      for (auto variation : variations){
        float value = variation->GetBinContent(j);
        if(value > highest) highest = value;
        if(value < lowest) lowest = value;
      }

      h_scale_up->SetBinContent(j, highest);
      h_scale_down->SetBinContent(j, lowest);

    }

    // lets plot how we make the envelope
    std::vector<int> colors = {2,3,4,6,7,8,9, 11, 12};
    TCanvas *can = new TCanvas("chist", "c", 600, 600);
    TLegend *leg = new TLegend(0.3,0.8,0.9,0.9);
    leg->SetNColumns(2);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    auto line = TLine(500.1,1,6000,1);

    int j = 0;
    for (auto variation : variations){
      auto variation_ratio = (TH1F*)variation->Clone();
      variation_ratio->Divide(h_nominal);

      variation_ratio->SetLineColor(colors.at(j));
      variation_ratio->SetLineWidth(3);
      variation_ratio->SetLineStyle(2);

      variation_ratio->GetYaxis()->SetRangeUser(0.4, 1.8);
      variation_ratio->GetXaxis()->SetRangeUser(500, 6000);
      variation_ratio->GetXaxis()->SetTitle( variation_ratio->GetTitle() );
      variation_ratio->GetYaxis()->SetTitle( "variation / nominal" );
      variation_ratio->SetTitle("");

      if(j == 0) variation_ratio->Draw("hist");
      else variation_ratio->Draw("hist same");

      leg->AddEntry(variation_ratio, variation_names.at(j), "l");

      j++;
    }

    line.Draw("same");
    leg->Draw();

    // Save the histo with the up/down variations in root file
    TFile* f_out = new TFile("/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/" + region + "_scale_" + year + "_" + channel + "_" + samples.at(i) + ".root", "RECREATE");
    h_scale_up->SetName(samples.at(i)+"_scale_up");
    h_scale_down->SetName(samples.at(i)+"_scale_down");
    h_scale_up->Write();
    h_scale_down->Write();
    delete f_out;

    // lets also draw the envelope
    auto h_scale_up_ratio = (TH1F*)h_scale_up->Clone();
    h_scale_up_ratio->Divide(h_nominal);
    auto h_scale_down_ratio = (TH1F*)h_scale_down->Clone();
    h_scale_down_ratio->Divide(h_nominal);

    h_scale_up_ratio->SetLineColor(1);
    h_scale_up_ratio->SetLineWidth(4);
    h_scale_up_ratio->SetLineStyle(1);
    h_scale_down_ratio->SetLineColor(1);
    h_scale_down_ratio->SetLineWidth(4);
    h_scale_down_ratio->SetLineStyle(1);

    can->SaveAs("plots/envelope_MCscale_noenvelope_" + region + "_" + year + "_" + channel + "_" + samples.at(i) +".pdf");

    h_scale_up_ratio->Draw("hist same");
    h_scale_down_ratio->Draw("hist same");

    can->SaveAs("plots/envelope_MCscale_" + region + "_" + year + "_" + channel + "_" + samples.at(i) +".pdf");

    delete f_in;
  }
}
