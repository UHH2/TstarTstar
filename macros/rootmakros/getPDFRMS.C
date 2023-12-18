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

void getPDFRMS( TString year = "UL18", TString channel = "mu", TString region = "SignalRegion"){

  // this assumes that all signals are replicas, and all others are hessians (only relevant for ttbar anyway atm.)
  // ST will be fixed later, need to check that sample then!
  TString filename_base = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/" + year + "/hadded/uhh2.AnalysisModuleRunner.MC.";

  std::cout << "Working on " << year << " for the " << channel << " channel in the " << region << "." << std::endl;

  vector<TString> samples = {"TTbar", "ST"};
  vector<bool> isSignal (samples.size(), false);
  vector<TString> masspoints = {"700", "800", "900", "1000", "1100", "1200", "1300", "1400", "1500", "1600", "1700", "1800", "1900", "2000", "2250", "2500", "2750"};
  for (auto mass : masspoints) {
    samples.push_back("TstarTstar_M-" + mass);
    isSignal.push_back(true);
    samples.push_back("TstarTstar_Spin32_M-" + mass);
    isSignal.push_back(true);
  }


  for(unsigned int i=0; i<samples.size(); i++){

    cout << "sample " << samples.at(i) << endl;

    TString filename = filename_base + samples.at(i) + ".root";
    TFile* f_in = new TFile(filename, "READ");

    // For each Signal sample read the normalization values - one for each of the 100 PDF replicas
    vector<double> pdf_norm (100, 1.); //For bkgs set normalization value to 1
    if( isSignal.at(i) ){
      string pdf_numb[100];
      ifstream normfile("/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/signalnorm/SignalNorm_" + year + "_" + samples.at(i) + ".txt", ios::in);
      if (normfile.is_open()){
        for(int j = 0; j < 100; ++j){
          normfile >> pdf_numb[j] >> pdf_norm[j];
        }
        normfile.close();
      } else {
        throw std::runtime_error("Signal norm file could not be opened.");
      }
    }

    TH1F *h_nominal = (TH1F*)f_in->Get(region + "_" + channel + "/pt_ST_nominal");
    TH1F *h_PDF_up = (TH1F*)h_nominal->Clone();
    TH1F *h_PDF_down = (TH1F*)h_nominal->Clone();

    float sum_bins = 0;
      // Loop over each bin of the Mtt histograms
     for(int j=1; j < h_nominal->GetXaxis()->GetNbins()+1; j++){

      float nominal = h_nominal->GetBinContent(j);
      float sum_bins = 0;

       // Loop over each of the 100 Histogrmas reweighted with the PDF replicas
       for(int i=0; i<100; i++){

           stringstream ss_name;
           ss_name << region + "_" + channel + "/pt_ST_PDF_" << i + 1;
           string s_name = ss_name.str();
           const char* char_name = s_name.c_str();

           float bin =  ((TH1F*)(f_in->Get(char_name)))->GetBinContent(j);
           float norm_bin = bin * pdf_norm[i];

           sum_bins += pow(norm_bin - nominal, 2);
      }

      float rms;
      if( isSignal.at(i) ){
        rms = sqrt( sum_bins / 100  );
      } else {
        rms = sqrt( sum_bins );
      }

      h_PDF_up->SetBinContent(j, nominal + rms);
      h_PDF_down->SetBinContent(j, nominal - rms);

    }

    // Save the histo with the up/down variations in root file
    TFile* f_out = new TFile("/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/PDF/" + region + "_PDF_" + year + "_" + channel + "_" + samples.at(i) + ".root", "RECREATE");
    h_PDF_up->SetName(samples.at(i)+"_PDF_up");
    h_PDF_down->SetName(samples.at(i)+"_PDF_down");
    h_PDF_up->Write();
    h_PDF_down->Write();
    delete f_out;

    // Plots nominal hist + up/down variations
    TCanvas* c = new TCanvas("c", "c", 1200, 800);
    c->Divide(1,1);
    c->cd(1);
    gPad->SetTopMargin(0.07);
    gPad->SetBottomMargin(0.17);
    gPad->SetLeftMargin(0.2);
    gPad->SetRightMargin(0.1);
    //gPad->SetLogy();
    gStyle->SetOptStat(0);

    h_PDF_up->SetLineWidth(1);
    h_PDF_up->SetLineColor(kRed);
    h_PDF_down->SetLineWidth(1);
    h_PDF_down->SetLineColor(kBlue);

    h_PDF_up->Divide(h_nominal);
    h_PDF_down->Divide(h_nominal);

    h_PDF_up->Draw("HIST");
    h_PDF_up->SetTitle("");
    h_PDF_up->GetXaxis()->SetTitle("S_{T} [GeV]");
    h_PDF_up->GetXaxis()->SetRangeUser(600.,6000.);
    h_PDF_up->GetXaxis()->SetTitleSize(0.055);
    h_PDF_up->GetXaxis()->SetLabelSize(0.05);
    h_PDF_up->GetYaxis()->SetTitle("Events");
    h_PDF_up->GetYaxis()->SetRangeUser(0.5,1.5);
    h_PDF_up->GetYaxis()->SetTitleSize(0.055);
    h_PDF_up->GetYaxis()->SetLabelSize(0.05);
    h_PDF_down->Draw("HIST same");

    // draw horizontal line
    auto line = TLine(600.1,1,6000,1);
    line.Draw("same");

    auto legend = new TLegend(0.7,0.67,0.85,0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.035);
    legend->SetTextFont(42);
    legend->AddEntry(h_PDF_up,"PDF up","le");
    legend->AddEntry(h_PDF_down,"PDF down","le");
    legend->Draw();

    TString cmstext = "CMS";
    TLatex *text2 = new TLatex(3.5, 24, cmstext);
    text2->SetNDC();
    text2->SetTextAlign(13);
    text2->SetX(0.24);
    text2->SetTextFont(62);
    text2->SetTextSize(0.06825);
    text2->SetY(0.895);
    text2->Draw();


    TString supptext = "Simulation";
    TLatex *text4 = new TLatex(3.5, 24, supptext);
    text4->SetNDC();
    text4->SetTextAlign(13);
    text4->SetX(0.24);
    text4->SetTextFont(52);
    text4->SetTextSize(0.55*0.06825);
    text4->SetY(0.8312);
    text4->Draw();

    TString supptext2 = "Work in progress";
    TLatex *text5 = new TLatex(3.5, 24, supptext2);
    text5->SetNDC();
    text5->SetTextAlign(13);
    text5->SetX(0.24);
    text5->SetTextFont(52);
    text5->SetTextSize(0.55*0.06825);
    text5->SetY(0.79);
    text5->Draw();

    c->SaveAs("/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/plots/" + region + "_PDF_" + year + "_" + channel + "_" + samples.at(i) + ".pdf");
    c->Close();

    delete f_in;
  }
}
