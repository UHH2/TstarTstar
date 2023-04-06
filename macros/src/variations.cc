#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/Utils.h"

#include "TSystem.h"

using namespace std;

void SetupCanvas(bool bPlotRatio=false);
void HistCosmetics(TH1* hist);
void PortraitCosmetics(TH1* hist, bool bPlotRatio);
void YieldCosmetics(TH1* hist);
void RatioCosmetics(TH1* hist,TString xtitle, int color, int style, int marker, int width);
void Cosmetics(TH1* hist, TString xtitle, int color, int style, int marker, int width, bool bPlotRatio);
void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre);
void DrawLegend(vector<TH1F*> hists, TString hname);
void CompareHistStructureLocal(TH1F* h1, TH1F* h2);
void DrawJMS(TF2* chi2, VecDD &points, VecDD &extrema, VecD &minimum, TString hist);
TH1F* GetRatioLocal(TH1F* h1, TH1F* h2, bool isData = false);

TCanvas* m_can;
TPad* m_pad1;
TPad* m_pad2;

TPad* m_rp1_top;
TPad* m_rp1;
TPad* m_rp2_top;
TPad* m_rp2;

bool norm = false;

int main(int argc, char* argv[]){

  norm = stob(argv[1]);

  SetupGlobalStyle();
  gStyle->SetErrorX(0);
  gStyle->SetOptTitle(0);

  vector<TString> years = {"UL18"};
  vector<vector<TString>> all  = {
    { "nominal",       "prefiringUp",    "prefiringDown"},
    { "nominal",       "sfelec_idUp",    "sfelec_idDown"},
    { "nominal",       "sfmu_idUp",      "sfmu_idDown"},
    { "nominal",       "sfmu_isolationUp",      "sfmu_isolationDown"},
  };
  vector<TString> LegendEntries = { "nominal", "up", "down"};
  vector<Color_t> colors =        {    kBlack, kAzure+7,   kAzure+7};
  vector<int>     styles =        {         1,        1,          2};

  int rebin = 5;

  TFile *file = new TFile("/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/UL18/hadded/uhh2.AnalysisModuleRunner.MC.TTbar.root");

  for(auto year: years){
    double lumi = 65.;
    for(unsigned int v=0; v<all.size(); v++){

      vector<TString> variations = all[v];
      cout << __LINE__ << endl;
      vector<TH1F*> hists, ratio_hists;
      for(TString var: variations){
        TH1F* h = (TH1F*) file->Get("SignalRegion_total/pt_ST_"+var);
        double integral = h->Integral();
        if(norm) h->Scale(1/integral);
        hists.push_back(h);
      }
      cout << __LINE__ << endl;

      for(auto h: hists) ratio_hists.push_back(GetRatioLocal(h, hists[0]));

      SetupCanvas(true);
      cout << __LINE__ << endl;

      for(unsigned int i=0; i<hists.size(); i++)            Cosmetics(hists[i],       "#tau_{32}", colors[i],  styles[i],  0, 2, true);
      for(unsigned int i=0; i<ratio_hists.size(); i++) RatioCosmetics(ratio_hists[i], "#tau_{32}", colors[i],  styles[i],  2, 0);

      cout << __LINE__ << endl;
      m_rp1_top->cd();
      gPad->SetLogy();
      TString title = variations[1];
      title.ReplaceAll("Up", "");
      hists[0]->SetTitle(title);
      hists[0]->GetYaxis()->SetRangeUser(0.1, hists[0]->GetMaximum()*1.2);
      hists[0]->Draw("hist ][");
      for(unsigned int i=1; i<hists.size(); i++) hists[i]->Draw("hist same ][");
      gPad->RedrawAxis();
      DrawLegend(hists, "tau32_"+title);

      // void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre);
      DrawLumi(lumi, true, false, false);

      m_rp1->cd();
      gPad->SetLogy();
      ratio_hists[0]->Draw("hist ][");
      for(unsigned int i=1; i<ratio_hists.size(); i++) ratio_hists[i]->Draw("hist same ][");
      gPad->RedrawAxis();

      if(norm) m_can->SaveAs("plots/"+title+"_"+year+"_norm.pdf");
      else     m_can->SaveAs("plots/"+title+"_"+year+".pdf");

    }
  }


  return 0;

}

// ===========================================================================
// ===========================================================================
// ===========================================================================

void SetupCanvas(bool bPlotRatio) // Copied from SPlotter.cxx
{
  // set up a canvas for single EPS files
  // optimised plots for including in theses or publications and documents
  // different possibilities
  // ratio/no ratio plots

  Int_t CanWidth;
  Int_t CanHeight;
  CanWidth = 400;
  CanHeight = 400;

  // set up the canvas
  m_can = new TCanvas("canvas","Control Plots", CanWidth, CanHeight);

  Float_t yplot = 0.65;
  Float_t yratio = 0.34;

  //  coordinates:
  // set up the coordinates of the two pads:    //
  Float_t y1, y2, y3;                           //  y3 +-------------+
  y3 = 0.99;                                    //     |             |
  y2 = y3-yplot;                                //     |     pad1    |
  y1 = y2-yratio;                               //  y2 |-------------|
  Float_t x1, x2;                               //     |     rp1     |
  x1 = 0.01;                                    //  y1 +-------------+
  x2 = 0.99;                                    //     x1            x2
  //
  // No Pad 2!


  m_rp1_top = new TPad("pad1", "Control Plots 2", x1, y2, x2, y3);
  m_rp1 = new TPad("rp1", "Ratio2", x1, y1, x2, y2);
  m_pad1 = new TPad("pad1", "Control Plots 2", x1, y1, x2, y3);

  m_rp2_top = new TPad("pad1", "Control Plots 2", x1, y2, x2, y3);
  m_rp2 = new TPad("rp1", "Ratio2", x1, y1, x2, y2);
  m_pad2 = new TPad("pad1", "Control Plots 2", x1, y1, x2, y3);


  m_pad1->SetTopMargin(0.05); m_pad1->SetBottomMargin(0.16);  m_pad1->SetLeftMargin(0.19); m_pad1->SetRightMargin(0.05);
  m_pad2->SetTopMargin(0.05); m_pad2->SetBottomMargin(0.16);  m_pad2->SetLeftMargin(0.19); m_pad2->SetRightMargin(0.05);

  m_rp1_top->SetTopMargin(0.065); m_rp1_top->SetBottomMargin(0.04);  m_rp1_top->SetLeftMargin(0.19); m_rp1_top->SetRightMargin(0.05);
  m_rp2_top->SetTopMargin(0.065); m_rp2_top->SetBottomMargin(0.04);  m_rp2_top->SetLeftMargin(0.19); m_rp2_top->SetRightMargin(0.05);

  m_rp1->SetTopMargin(0.0);    m_rp1->SetBottomMargin(0.35);  m_rp1->SetLeftMargin(0.19);  m_rp1->SetRightMargin(0.05);
  m_rp2->SetTopMargin(0.0);    m_rp2->SetBottomMargin(0.35);  m_rp2->SetLeftMargin(0.19);  m_rp2->SetRightMargin(0.05);

  m_pad1->Draw();
  m_pad2->Draw();

  m_rp1_top->Draw();
  // m_rp2_top->Draw();

  if (bPlotRatio){
    m_rp1->Draw();
    // m_rp2->Draw();
  }

  return;

}

// =============================================================================
// ===                                                                       ===
// =============================================================================

void YieldCosmetics(TH1* hist)
{
  // cosmetics for the lumi yield histogram
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetLabelOffset(0.008);
  hist->GetXaxis()->SetTickLength(0.03);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.2);

  hist->GetYaxis()->SetTitleOffset(1.8);
  hist->GetYaxis()->SetTitleSize(0.055);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTickLength(0.02);
  hist->GetYaxis()->SetLabelOffset(0.011);

  hist->GetXaxis()->SetTitle("integrated luminosity [fb^{-1}]");
  double dlum = hist->GetXaxis()->GetBinWidth(1);
  TString xtit = TString::Format("events per %3.1f fb^{-1}", dlum);
  hist->GetYaxis()->SetTitle(xtit);
}

// =============================================================================
// ===                                                                       ===
// =============================================================================

void PortraitCosmetics(TH1* hist, bool bPlotRatio)
{

  // top histogram of the ratio plot
  if (bPlotRatio){

    // x-axis
    hist->GetXaxis()->SetTickLength(0.05);
    hist->GetXaxis()->SetLabelSize(0.00); //remove labels from Xaxis of upper plot
    hist->GetXaxis()->SetTitleSize(0.00); // remove title from Xaxis

    // y-axis
    hist->GetYaxis()->SetTitleSize(0.07);
    hist->GetYaxis()->SetLabelSize(0.062);
    hist->GetYaxis()->SetLabelOffset(0.01);
    hist->GetYaxis()->SetTitleOffset(0.8);
    hist->GetYaxis()->SetTickLength(0.02);

    // only this histogram
  } else {

    hist->GetXaxis()->SetLabelSize(0.045);
    hist->GetXaxis()->SetLabelOffset(0.008);
    hist->GetXaxis()->SetTickLength(0.03);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.2);

    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.045);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetLabelOffset(0.011);

  }

}

// =============================================================================
// ===                                                                       ===
// =============================================================================

void Cosmetics(TH1* hist,TString xtitle, int color, int style, int marker, int width, bool bPlotRatio)
{
  PortraitCosmetics(hist, bPlotRatio);

  hist->SetLineColor(color);
  hist->SetLineStyle(style);
  hist->SetLineWidth(width);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(0.7);

  // set Y-axis title
  hist->GetYaxis()->SetTitle(norm?"a.u.":"Events");
  if(xtitle.Contains("JEC")) hist->GetYaxis()->SetTitle("f^{XCone}");

  // set X-axis title
  hist->GetXaxis()->SetTitle(xtitle);
  hist->SetTitle("");

  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetLabelFont(42);

  // top histogram of the ratio plot
  if (bPlotRatio){

    // x-axis
    hist->GetXaxis()->SetTickLength(0.05);

    // y-axis
    hist->GetYaxis()->SetTitleSize(0.07);
    hist->GetYaxis()->SetLabelSize(0.062);
    hist->GetYaxis()->SetLabelOffset(0.01);
    hist->GetYaxis()->SetTitleOffset(1.35);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetNdivisions(505);

    TString hname = hist->GetName();
    hist->GetXaxis()->SetNdivisions(505);

    // only this histogram
  } else {

    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetXaxis()->SetLabelOffset(0.008);
    hist->GetXaxis()->SetTickLength(0.03);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.2);

    hist->GetYaxis()->SetTitleOffset(1.8);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetLabelOffset(0.011);

    hist->GetYaxis()->SetTitleOffset(1.8);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.045);
    hist->GetYaxis()->SetTickLength(0.02);
    hist->GetYaxis()->SetLabelOffset(0.011);
  }
}

// =============================================================================
// ===                                                                       ===
// =============================================================================

void RatioCosmetics(TH1* hist,TString xtitle, int color, int style, int width, int marker)
{
  hist->SetLineColor(color);
  hist->SetLineStyle(style);
  hist->SetLineWidth(width);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(0.7);

  hist->GetYaxis()->SetRangeUser(0.7, 1.5);
  if(xtitle.Contains("#it{m}_{W}")) hist->GetYaxis()->SetRangeUser(0.7, 1.3);
  hist->GetXaxis()->SetNdivisions(505);
  hist->SetMarkerSize(0.7);

  hist->SetTitle("");
  // x-axis
  hist->GetXaxis()->SetTitle(xtitle);
  hist->GetXaxis()->SetLabelSize(0.12);
  hist->GetXaxis()->SetTickLength(0.08);
  hist->GetXaxis()->SetTitleSize(0.12);
  hist->GetXaxis()->SetTitleOffset(1.25);
  hist->GetXaxis()->SetLabelOffset(0.02);
  TString hname = hist->GetName();
  hist->GetXaxis()->SetNdivisions(505);

  // y-axis
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitle("Ratio");
  hist->GetYaxis()->SetTitleSize(0.12);
  hist->GetYaxis()->SetTitleOffset(0.78); //0.66
  hist->GetYaxis()->SetLabelSize(0.11);
  //hist->GetYaxis()->SetNdivisions(210);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetTickLength(0.02);
  hist->GetYaxis()->SetLabelOffset(0.011);

}

// =============================================================================
// ===                                                                       ===
// =============================================================================

void DrawLumi(double m_lumi, bool bPlotRatio, bool DrawCMS, bool forPre)
{
  TString infotext = TString::Format("%3.0f fb^{-1} (13 TeV)", m_lumi);
  TLatex *text1 = new TLatex(3.5, 24, infotext);
  text1->SetNDC();
  text1->SetTextAlign(33);
  text1->SetX(0.95);
  text1->SetTextFont(42);
  if (bPlotRatio){
    text1->SetTextSize(0.06);
    text1->SetY(1.);
  } else {
    text1->SetTextSize(0.045);
    text1->SetX(0.90);
    text1->SetY(0.95);
  }
  text1->Draw();

  if (DrawCMS || forPre){
    TString cmstext = "CMS";
    TLatex *text2 = new TLatex(3.5, 24, cmstext);
    text2->SetNDC();
    text2->SetTextAlign(13);
    text2->SetX(0.22); // standard was 0.24
    text2->SetTextFont(62);
    if (bPlotRatio){
      text2->SetTextSize(0.08);
      text2->SetY(0.87);
    } else {
      text2->SetTextSize(0.05);
      text2->SetY(0.945);
    }
    text2->Draw();
  }

  if (forPre){
    TString preltext = "Work in Progress"; // Preliminary
    TLatex *text3 = new TLatex(3.5, 24, preltext);
    text3->SetNDC();
    text3->SetTextAlign(13);
    text3->SetX(0.22); // standard was 0.24
    text3->SetTextFont(52);
    if (bPlotRatio){
      //  text3->SetTextSize(0.06);
      text3->SetTextSize(0.065);
      text3->SetY(0.78);
    } else {
      // text3->SetTextSize(0.035);
      text3->SetTextSize(0.04);
      text3->SetX(0.33); // standard was 0.24
      text3->SetY(0.94);
    }
    text3->Draw();
  }

}

// =============================================================================
// ===                                                                       ===
// =============================================================================

void DrawLegend(vector<TH1F*> hists, TString hname)
{
  // draw a legend

  int narr = hists.size();
  float yfrac = 0.07;

  float top = 0.85;

  float ysize = yfrac*narr;
  float xleft = 0.64;
  float xright = 0.92;

  xright = xleft + 0.22;
  top += 0.02;

  TLegend *leg = new TLegend(xleft,top-ysize,xright,top, NULL,"brNDC");
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);

  for (Int_t i=0; i<narr; ++i){
    TH1F* sh = hists[i];

    TString legname = TString::Format("leg_entry_%i",i);
    TString legtitle;
    TString var = hname; var.ReplaceAll("tau32_", "");
    if(i==0) legtitle = "t#bar{t}";

    if(var.Contains("qcdbased") || var.Contains("gluonmove")){
      if(i==1) legtitle = var;
    }
    else{
      if(i==1) legtitle = var+" up";
      if(i==2) legtitle = var+" down";
    }


    TLegendEntry* entry = NULL;
    int marker = sh->GetMarkerStyle();
    int lstyle = sh->GetLineStyle();

    if (marker>0){
      entry = leg->AddEntry(sh, legtitle, "pe");
      entry->SetLineWidth(1);
      entry->SetLineColor(sh->GetLineColor());
      entry->SetMarkerColor(sh->GetLineColor());
      entry->SetMarkerStyle(marker);
      entry->SetMarkerSize(1.0);
      entry->SetMarkerSize(0.8);
    } else {

      entry = leg->AddEntry(sh, legtitle, "l");
      entry->SetLineColor(sh->GetLineColor());
      entry->SetMarkerStyle(0);
      entry->SetMarkerSize(0);
      entry->SetMarkerColor(sh->GetLineColor());
      entry->SetLineWidth(2);
      entry->SetLineStyle(lstyle);

      entry->SetTextAlign(12);
    }
  }


  leg->Draw();

}

// =============================================================================
// ===                                                                       ===
// =============================================================================

void CompareHistStructureLocal(TH1F* h1, TH1F* h2){
  // Check if histograms have the same structure
  bool sameBins = h1->GetNbinsX() == h2->GetNbinsX();
  bool sameCenter = true;
  for(unsigned int i=1; i<h1->GetNbinsX(); i++){
    sameCenter = h1->GetBinCenter(i) == h2->GetBinCenter(i);
    if(!sameCenter) break;
  }
  if(!sameBins||!sameCenter) throw runtime_error("GetRatioLocal: Histograms do not have the same structure!");
}

TH1F* GetRatioLocal(TH1F* h1, TH1F* h2, bool isData){
  CompareHistStructureLocal(h1, h2);
  TH1F* ratio = (TH1F*) h1->Clone();
  int Nbins = h1->GetNbinsX();
  for(int i=1; i<=Nbins;i++){
    double N1 = h1->GetBinContent(i);
    double N2 = h2->GetBinContent(i);
    double E1 = h1->GetBinError(i);
    double E2 = h2->GetBinError(i);
    if(N1==0 || N2==0){
      ratio->SetBinContent(i, isData?10:1);
      ratio->SetBinError(i, 0);
    }
    else{
      double r = N1/N2;
      double error = sqrt(E1/N2 * E1/N2 + N1*E2/(N2*N2) * N1*E2/(N2*N2));
      ratio->SetBinContent(i, r);
      ratio->SetBinError(i, error);
    }
  }
  return ratio;
}
