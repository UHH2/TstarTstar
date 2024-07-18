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

void combinePDFRMSforST( TString year = "UL18", TString channel = "mu", TString region = "SignalRegion"){

    // this will combine the result histograms for the PDF variations in ST, assuming they can just be "hadded"

    // TFile* f_out = new TFile("/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/PDF/" + region + "_PDF_" + year + "_" + channel + "_" + samples.at(i) + ".root", "RECREATE");

    TString base_path = "/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/PDF/" + region + "_PDF_" + year + "_" + channel + "_"; 

    // getting the first histograms
    std::cout << "Getting first..." << std::endl;
    TFile* f_in_A = new TFile(base_path + "ST_s-channel.root", "READ");
    TH1F *h_PDF_up = (TH1F*)f_in_A->Get("ST_s-channel_PDF_up")->Clone();
    TH1F *h_PDF_down = (TH1F*)f_in_A->Get("ST_s-channel_PDF_down")->Clone();

    // adding the second histograms
    std::cout << "Adding second..." << std::endl;
    TFile* f_in_B = new TFile(base_path + "ST_others.root", "READ");
    h_PDF_up->Add( (TH1F*)f_in_B->Get("ST_others_PDF_up") );
    h_PDF_down->Add( (TH1F*)f_in_B->Get("ST_others_PDF_down") );

    // outputting
    std::cout << "Outputting..." << std::endl;
    TFile* f_out = new TFile(base_path + "ST" + ".root", "RECREATE");
    h_PDF_up->SetName("ST_PDF_up");
    h_PDF_down->SetName("ST_PDF_down");
    h_PDF_up->Write();
    h_PDF_down->Write();
    delete f_out;

}