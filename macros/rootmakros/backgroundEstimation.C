// author: F.Labe
// date: 24.09.2021
// Run it with following command:
// root -l -b -q backgroundEstimation.C

template<class T, size_t N>
constexpr size_t size(T (&)[N]) { return N; }

TF1 *fit1, *fit2, *fitMean, *btagUp, *btagDown; // yes these need to be defined out here. Why? because root is absolutely stupid
Double_t meanFunc(Double_t *x, Double_t *par) {
  return ( fit1->EvalPar(x,par) + fit2->EvalPar(x,par) ) / 2;
}

Double_t fitRatio1(Double_t *x, Double_t *par) {
  return ( fit1->EvalPar(x,par) / fitMean->EvalPar(x,par) );
}

Double_t fitRatio2(Double_t *x, Double_t *par) {
  return ( fit2->EvalPar(x,par) / fitMean->EvalPar(x,par) );
}

Double_t btagRatioUp(Double_t *x, Double_t *par) {
  return ( btagUp->EvalPar(x,par) / fitMean->EvalPar(x,par) );
}

Double_t btagRatioDown(Double_t *x, Double_t *par) {
  return ( btagDown->EvalPar(x,par) / fitMean->EvalPar(x,par) );
}

double GetXForHighestY(TGraphAsymmErrors graph) {
    if (graph.GetN() == 0) {
        // Handle invalid input or empty graph
        return 0.0; // Or any appropriate default value
    }

    double maxY = -1e10; // Initialize to a very small value
    double maxX = -1e10; // Initialize to a very small value

    int nPoints = graph.GetN();
    double* xValues = graph.GetX();
    double* yValues = graph.GetY();

    for (int i = 0; i < nPoints; ++i) {
        if (yValues[i] > maxY) {
            maxY = yValues[i];
            maxX = xValues[i];
        }
    }

    return maxX;
}

// no year means full run 2
// total channel means combination of both ele and mu
// JE_string for example "_JECUp" ATTENTION MUST BE HADDED MANUALLY IN DATA FOLDER
// systematic for example btagging_totalUp. empty string means do not do any systematic 
void backgroundEstimation(TString channel = "mu", TString region = "SR", TString systematic = "", bool storeOutputToFile = false, bool plot_other_ratios = false){ // keep the last two to (true, false) when runnning in batch mode!

  TString year = "";
  TString JE_string = "";
  bool plot_stat_unc = !plot_other_ratios;

  // definitions
  std::vector<TString> nontop_backgrounds = {"WJets", "QCD", "VV", "DYJets"};
  std::vector<TString> top_backgrounds = {"ST", "TTbar"};

  TString subpath_SR;
  if (region == "VR") subpath_SR = "ValidationRegion_" + channel;
  else if (region == "SR") subpath_SR = "SignalRegion_" + channel;
  TString subpath_CR = "ControlRegion_" + channel;  
  TString histname = "pt_ST_nominal";
  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/";
  TString fileprefix = "uhh2.AnalysisModuleRunner.";

  path = path + "/" + year + "/hadded" + JE_string + "/";

  std::cout << "Using path: " << path << std::endl;

  TH1D *histSR_nontop;
  TH1D *hist_btagCR_nontop;
  TH1D *hist_btagCR_top;
  TH1D *hist_btagCR;

  double sizeTick = 12;

  if(systematic != "") {
    histname = "pt_ST_" + systematic;
    systematic = "_" + systematic;
  }

  // this is only used if we plot the baseline!
  TString other_ratios_base_path = "/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/bgest/";
  std::vector<TString> other_ratios = {"btagging_total"};

  // open full set of non-top backgrounds in signal region
  bool first = true;
  for (auto filename : nontop_backgrounds) {
    TFile *input = TFile::Open(path+fileprefix+"MC."+filename+".root");
    if(!input) cout << "Empty file" << endl;
    TH1D *hist = (TH1D*)input->Get(subpath_SR+"/"+histname);
    if(!hist) cout << "Empty hist" << endl;
    if (first) {
      first = false;
      histSR_nontop = hist;
    } else {
      histSR_nontop->Add(hist);
    }
  }

  // open full set of non-top backgrounds in W-jets control region
  first = true;
  vector<double> nontop_integrals;
  for (auto filename : nontop_backgrounds) {
    TFile *input = TFile::Open(path+fileprefix+"MC."+filename+".root");
    if(!input) cout << "Empty file" << endl;
    TH1D *hist = (TH1D*)input->Get(subpath_CR+"/"+histname);
    if(!hist) cout << "Empty hist" << endl;
    // outputting some information on the b-tagging CR
    nontop_integrals.push_back(hist->Integral());
    if (first) {
      first = false;
      hist_btagCR_nontop = hist;
    } else {
      hist_btagCR_nontop->Add(hist);
    }
  }

  // open the top backgrounds
  first = true;
  vector<double> top_integrals;
  for (auto filename : top_backgrounds) {
    TFile *input = TFile::Open(path+fileprefix+"MC."+filename+".root");
    if(!input) cout << "Empty file" << endl;
    TH1D *hist = (TH1D*)input->Get(subpath_CR+"/"+histname);
    if(!hist) cout << "Empty hist" << endl;
    // outputting some information on the b-tagging CR
    top_integrals.push_back(hist->Integral());
    if (first) {
      first = false;
      hist_btagCR_top = hist;
    } else {
      hist_btagCR_top->Add(hist);
    }
  }

  // combining for full picture
  hist_btagCR = (TH1D*) hist_btagCR_nontop->Clone();
  hist_btagCR->Add(hist_btagCR_top);

  for(uint i = 0; i < nontop_backgrounds.size(); i++) std::cout << "Contribution of " << nontop_backgrounds.at(i) << ": " << nontop_integrals.at(i)/hist_btagCR->Integral() << std::endl;
  for(uint i = 0; i < top_backgrounds.size(); i++) std::cout << "Contribution of " << top_backgrounds.at(i) << ": " << top_integrals.at(i)/hist_btagCR->Integral() << std::endl;

  TGraphAsymmErrors purity = TGraphAsymmErrors();
  purity.Divide(hist_btagCR_nontop, hist_btagCR, "cl=0.68 b(1,1) mode");

  if(storeOutputToFile) {
    TString filename = "files/bgest/purity_HOTVR_" + year + "_" + channel + systematic + JE_string + ".root";
    TFile *file = new TFile(filename, "RECREATE");
    purity.SetName("purity");
    purity.Write();
    file->Close();
    delete file;
  }
  
  // calculate ratio histogram
  TGraphAsymmErrors ratio = TGraphAsymmErrors();
  ratio.Divide(histSR_nontop, hist_btagCR_nontop, "pois");

  // fitting landau
  fit1 = new TF1("fit1", "landau", 0, 6000);
  ratio.Fit("fit1", "N", "", 500, 6000);

  TH1D *fit1unc = new TH1D("fit1unc", "Fit 1 with conf.band", 500, 0, 6000);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(fit1unc, 0.68);

  TH1D *fit2unc = new TH1D("fit2unc", "Fit 2 with conf.band", 500, 0, 6000);
  fit2 = new TF1("fit2", "[3] + [0] * exp( -0.5 * (( x - [1])/[2])^2) ", 0, 6000);

  // adapting starting parameters
  if (channel == "mu" && region == "VR") {
    fit2->SetParameters(-0.5, 500, 2000, 0.5);
  } else if (channel == "mu" && region == "SR") {
    fit2->SetParameters(-1, 500, 2000, 1);
  } else if (channel == "ele" && region == "VR") {
    fit2->SetParameters(1, 1000, 1000, 0);
  } else {
    fit2->SetParameters(0.1, 3000, 1000, 1);
  }
  ratio.Fit("fit2", "N", "", 500, 6000);

  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(fit2unc, 0.68);

  // defining average
  fitMean = new TF1("mean", meanFunc, 0, 6000);

  // plot everything
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
  gStyle->SetPadBottomMargin(0.3);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.1);

  gROOT->ForceStyle();
  Double_t w = 600;
  Double_t h = 600;

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);

  TPad *pad1 = new TPad("pad1", "The pad 80% of the height", 0.0, 0.345, 1.0, 1.0);
  double pad1W = pad1->GetWw()*pad1->GetAbsWNDC();
	double pad1H = pad1->GetWh()*pad1->GetAbsHNDC();
  TPad *pad2 = new TPad("pad2", "The pad 20% of the height", 0.0, 0.0, 1.0, 0.34);
  double pad2W = pad2->GetWw()*pad2->GetAbsWNDC();
	double pad2H = pad2->GetWh()*pad2->GetAbsHNDC();

  pad2->Draw();
  pad1->Draw();

  pad1->cd();

  pad1->SetBottomMargin(0.02);
  pad2->SetTopMargin(0);

  pad1->SetTickx();
  pad1->SetTicky();

  double ystart = 0.5;
  if (plot_other_ratios) ystart = 0.45;
  auto legend = new TLegend(0.2,ystart ,0.5,0.8);
  legend->SetBorderSize(0);
  gStyle->SetLegendTextSize(0.04);

  const int nbins = 34;
  double bins[nbins] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500,
    2600, 2700, 2800, 2900, 3000, 3250, 4000, 6000};
  TH1F *dummy = new TH1F("dummy", "", nbins-1, bins);

  // draw plot
  // using a histogram to define canvas, as graphs are annoying
  dummy->GetXaxis()->SetLabelSize(0);
  dummy->GetXaxis()->SetTitleSize(0);
  dummy->GetXaxis()->SetNdivisions(505, kTRUE);
  if (region == "SR") dummy->GetYaxis()->SetTitle("SR / CR");
  else dummy->GetYaxis()->SetTitle("VR / SR");
  dummy->SetTitle("");
  dummy->GetYaxis()->SetRangeUser(0, 2);
  dummy->GetXaxis()->SetRangeUser(0, 6000);
  if(region == "SR") dummy->GetYaxis()->SetRangeUser(0, 0.5);
  dummy->Draw("hist");

  pad1->Update();
  double tickScaleX = (pad1->GetUxmax() - pad1->GetUxmin())/(pad1->GetX2()-pad1->GetX1())*pad1H;
	double tickScaleY = (pad1->GetUymax() - pad1->GetUymin())/(pad1->GetY2()-pad1->GetY1())*pad1W;
	dummy->GetXaxis()->SetTickLength(sizeTick/tickScaleX);
	dummy->GetYaxis()->SetTickLength(sizeTick/tickScaleY);
  dummy->Draw("hist");

  ratio.SetTitle("");
  ratio.SetLineColor(1);
  ratio.SetMarkerStyle(20);

  ratio.Draw("PZ same");

  fit1->SetLineColor(1);
  fit1->SetLineStyle(2);
  fit1->Draw("same");
  fit2->SetLineColor(1);
  fit2->SetLineStyle(3);
  fit2->Draw("same");
  fitMean->SetLineColor(1);
  fitMean->Draw("same");

  TH1D* uncertainty_fits_up = (TH1D*) fit1unc->Clone(); // just to get the structure
  TH1D* uncertainty_fits_down = (TH1D*) fit1unc->Clone(); // just to get the structure

  TH1D* uncertainty_fits = (TH1D*) fit1unc->Clone(); // just to get the structure

  if (plot_stat_unc) {

    fit1unc->SetFillColorAlpha(1, 0.3);
    //fit1unc->Draw("e3 same");
    fit2unc->SetFillColorAlpha(1, 0.3);
    //fit2unc->Draw("e3 same");    

    for (Int_t bin = 1; bin <= fit1unc->GetNbinsX(); ++bin) {
        double bin1 = fit1unc->GetBinContent(bin);
        double bin2 = fit2unc->GetBinContent(bin);

        double binError1 = fit1unc->GetBinError(bin);
        double binError2 = fit2unc->GetBinError(bin);

        double mean = (bin1 + bin2) / 2;
        double uncertainty = TMath::Sqrt((binError1 * binError1 + binError2 * binError2) / 4);

        uncertainty_fits->SetBinContent(bin, mean);
        uncertainty_fits->SetBinError(bin, uncertainty);

        uncertainty_fits_up->SetBinContent(bin, mean + uncertainty);
        uncertainty_fits_down->SetBinContent(bin, mean - uncertainty);
    }

    uncertainty_fits->SetFillColorAlpha(1, 0.3);
    uncertainty_fits->Draw("e3 same");

    /**
    uncertainty_fits_up->SetLineColor(4);
    uncertainty_fits_down->SetLineColor(4);
    uncertainty_fits_up->SetFillColorAlpha(0,0);
    uncertainty_fits_down->SetFillColorAlpha(0,0);
    uncertainty_fits_up->Draw("hist same");
    uncertainty_fits_down->Draw("hist same");
    **/

  }

  if(storeOutputToFile) {
    TFile *output;
    TString filename = "files/bgest/alphaFunction_HOTVR_" + year + "_" + region + "_" + channel + systematic + JE_string +"_fitstat.root";
    output = TFile::Open(filename, "RECREATE");
    uncertainty_fits_up->SetName("fitstat_up");
    uncertainty_fits_up->Write();
    uncertainty_fits_down->SetName("fitstat_down");
    uncertainty_fits_down->Write();
  }

  legend->AddEntry(&ratio, "Nontop simulation", "lep");

  if (plot_other_ratios) {

    int color_int = 2;
    for (auto name : other_ratios) {
      TFile *file_up = new TFile(other_ratios_base_path + "alphaFunction_HOTVR_" + year + "_" + region + "_" + channel + "_" + name + "Up.root");
      TGraphAsymmErrors *graph_up = (TGraphAsymmErrors*) file_up->Get("alpha_ratio")->Clone();
      TF1 *func_up = (TF1*) file_up->Get("fit_mean")->Clone();
      TFile *file_down = new TFile(other_ratios_base_path + "alphaFunction_HOTVR_" + year + "_" + region + "_" + channel + "_" + name + "Down.root");
      TGraphAsymmErrors *graph_down = (TGraphAsymmErrors*) file_down->Get("alpha_ratio")->Clone();
      TF1 *func_down = (TF1*) file_down->Get("fit_mean")->Clone();

      graph_up->SetLineColorAlpha(color_int, 0.5);
      graph_down->SetLineColorAlpha(color_int, 0.5);
      graph_up->SetMarkerColorAlpha(color_int, 0.5);
      graph_down->SetMarkerColorAlpha(color_int, 0.5);
      graph_up->SetMarkerStyle(20);
      graph_down->SetMarkerStyle(20);

      func_up->SetLineColorAlpha(color_int, 0.5);
      func_down->SetLineColorAlpha(color_int, 0.5);
      func_up->SetLineStyle(1);
      func_down->SetLineStyle(1);
      
      graph_up->Draw("P same");
      graph_down->Draw("P same");
      if(name == "btagging_total") {
        func_up->Draw("same");
        func_down->Draw("same");
      }

      TString legendname = name;
      if (legendname == "btagging_total") legendname = "b-tagging total";

      legend->AddEntry(graph_up, "#alpha (" + legendname + ")" , "ep");
      color_int++;
    }

  }

  // fit results
  TString fittxt = TString::Format("#chi^{2}/ndf: %3.2f / %3.0d", fit1->GetChisquare(), fit1->GetNDF());

  // fit results
  TString fit2txt = TString::Format("#chi^{2}/ndf: %3.2f / %3.0d", fit2->GetChisquare(), fit2->GetNDF());

  // legend
  
  legend->AddEntry(fit1,"Landau fit (" + fittxt + ")","l");
  legend->AddEntry(fit2,"Gauss + c fit (" + fit2txt + ")","l");
  
  
  //legend->AddEntry(fitMean,"mean fit","l");
  TGraph *dummyGraph = new TGraph();
  dummyGraph->SetLineColor(1);
  dummyGraph->SetLineWidth(2);
  dummyGraph->SetFillColorAlpha(1, 0.3);
  dummyGraph->SetFillStyle(1001);

  // Add the dummy TGraph to the legend with the fill representation
  legend->AddEntry(dummyGraph, "g_{TF}(S_{T})", "lf");
  
  legend->Draw();


  // draw Lumi text
  /**
  TString infotext = TString::Format("%3.0f fb^{-1} (%d TeV)", 138., 13);
  TLatex *text = new TLatex(3.5, 24, infotext);
  text->SetNDC();
  text->SetTextAlign(33);
  text->SetX(0.88);
  text->SetTextFont(42);
  text->SetY(1);
  text->SetTextSize(0.045);
  text->Draw();
  **/

  double CMS_size = 0.08;
  double CMS_ratio = 1.3;

  TString cmstext = "CMS";
  TLatex *text2 = new TLatex(3.5, 24, cmstext);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.2);
  text2->SetTextFont(62);
  text2->SetTextSize(CMS_size);
  text2->SetY(0.9);
  text2->Draw();

  TString preltext = "Simulation";
  TLatex *text3 = new TLatex(3.5, 24, preltext);
  text3->SetNDC();
  text3->SetTextAlign(13);
  text3->SetX(0.315);
  text3->SetTextFont(52);
  text3->SetTextSize(CMS_size / CMS_ratio);
  text3->SetY(0.887);
  text3->Draw();

  pad2->cd();
  pad2->SetTickx();
  pad2->SetTicky();

  // we'll need to loop over each entry of the ratio histogram
  // for each entry, we'll get the two values for each of the functions
  // we'll define the average as 0, and determine an up- and down variation (relative) from the two values
  // in quadrature, I'll add the decorrelation uncertainty to this
  // as well as the statistical uncertainty of the data points

  // code-wise, we'll fill two histograms here:
  // one containing 1 in each bin, but has the appropriate uncertainties
  // a second one containing the "data points", being the ratio of the ratio and the function average

  TH1F *deviation = new TH1F("deviation", "", nbins-1, bins);

  // the main loop
  int bins_before_600 = 0; // need to start at 0 because of overflow bin
  for (int bin = 0; bin <= nbins; bin++) {
    double st = deviation->GetBinCenter(bin);
    if( st < 500 ) {
      bins_before_600++;
      deviation->SetBinContent(bin, 0);
    }
    else {
      double func_val_1 = fit1->Eval(st);
      double func_val_2 = fit2->Eval(st);
      double func_val_avg = fitMean->Eval(st);

      double val_ratio = ratio.GetY()[bin-bins_before_600];
      double val_error = ratio.GetErrorY(bin-bins_before_600);

      // calculating how much the ratio deviates from the central fit
      deviation->SetBinContent(bin, val_ratio / func_val_avg);
      deviation->SetBinError(bin, val_error / func_val_avg);
    }
    
  }

  deviation->GetXaxis()->SetTitleSize(0.12);
  deviation->GetXaxis()->SetLabelSize(0.095);
  deviation->GetYaxis()->SetTitleSize(0.11);
  deviation->GetYaxis()->SetTitleOffset(0.65);
  deviation->GetYaxis()->SetLabelSize(0.095);

  deviation->GetXaxis()->SetNdivisions(505, kTRUE);
  deviation->GetYaxis()->SetNdivisions(505, kTRUE);
  deviation->GetYaxis()->SetRangeUser(0.5, 1.5);
  deviation->GetXaxis()->SetRangeUser(0, 6000);
  deviation->SetMarkerStyle(8);
  deviation->SetLineColor(1);
  deviation->GetXaxis()->SetTitle("S_{T} [GeV]");
  deviation->GetYaxis()->SetTitle("Nontop / g_{TF}(S_{T})");
  deviation->Draw("P");

  pad2->Update();
	tickScaleX = (pad2->GetUxmax() - pad2->GetUxmin())/(pad2->GetX2()-pad2->GetX1())*pad2H;
	tickScaleY = (pad2->GetUymax() - pad2->GetUymin())/(pad2->GetY2()-pad2->GetY1())*pad2W;
	deviation->GetXaxis()->SetTickLength(sizeTick/tickScaleX);
	deviation->GetYaxis()->SetTickLength(sizeTick/tickScaleY);
	deviation->Draw("P");

  ///// large code block only for variation plots!!!

  if (plot_other_ratios) {

    int color_int = 2;
    for (auto name : other_ratios) {
      TH1F *deviationUp = new TH1F("deviation up " + name , "", nbins-1, bins);
      TH1F *deviationDown = new TH1F("deviation down " + name, "", nbins-1, bins);

      TFile *file_up = new TFile(other_ratios_base_path + "alphaFunction_HOTVR_" + year + "_" + region + "_" + channel + "_" + name + "Up.root");
      TGraphAsymmErrors *graph_up = (TGraphAsymmErrors*) file_up->Get("alpha_ratio")->Clone();
      if(name == "btagging_total") btagUp = (TF1*) file_up->Get("fit_mean")->Clone();
      TFile *file_down = new TFile(other_ratios_base_path + "alphaFunction_HOTVR_" + year + "_" + region + "_" + channel + "_" + name + "Down.root");
      TGraphAsymmErrors *graph_down = (TGraphAsymmErrors*) file_down->Get("alpha_ratio")->Clone();
      if(name == "btagging_total") btagDown = (TF1*) file_down->Get("fit_mean")->Clone();

      // plotting the fit functions if we are btagging_total
      deviationUp->SetLineColorAlpha(color_int, 0.5);
      deviationDown->SetLineColorAlpha(color_int, 0.5);
      deviationUp->SetMarkerColorAlpha(color_int, 0.5);
      deviationDown->SetMarkerColorAlpha(color_int, 0.5);
      deviationUp->SetMarkerStyle(20);
      deviationDown->SetMarkerStyle(20);

      TF1 *btagRatioUp_func, *btagRatioDown_func;

      if(name == "btagging_total") {
        btagRatioUp_func = new TF1("btagRatioUp", btagRatioUp, 0, 6000);
        btagRatioDown_func = new TF1("btagRatioDown", btagRatioDown, 0, 6000);

        btagRatioUp_func->SetLineColorAlpha(color_int, 0.5);
        btagRatioDown_func->SetLineColorAlpha(color_int, 0.5);
        btagRatioUp_func->SetLineStyle(1);
        btagRatioDown_func->SetLineStyle(1);
      }

      for (int bin = 0; bin <= nbins; bin++) {
        double st = deviation->GetBinCenter(bin);
        if( st < 500 ) {
          deviationUp->SetBinContent(bin, 0);
          deviationDown->SetBinContent(bin, 0);
        }
        else {
          double func_val_1 = fit1->Eval(st);
          double func_val_2 = fit2->Eval(st);
          double func_val_avg = fitMean->Eval(st);

          double val_ratio_up = graph_up->GetY()[bin-bins_before_600];
          double val_error_up = graph_up->GetErrorY(bin-bins_before_600);

          double val_ratio_down = graph_down->GetY()[bin-bins_before_600];
          double val_error_down = graph_down->GetErrorY(bin-bins_before_600);

          // calculating how much the ratio deviates from the central fit
          deviationUp->SetBinContent(bin, val_ratio_up / func_val_avg);
          deviationUp->SetBinError(bin, val_error_up / func_val_avg);
          deviationDown->SetBinContent(bin, val_ratio_down / func_val_avg);
          deviationDown->SetBinError(bin, val_error_down / func_val_avg);
        }
        
      }
 
      if(name == "btagging_total") {
        btagRatioUp_func->Draw("same");
        btagRatioDown_func->Draw("same");
      }

      deviationUp->Draw("P same");
      deviationDown->Draw("P same");
      color_int++;

    }

  }
  //end of large code block only for variation plots

  // drawing a line through 0
  TLine line = TLine(0,1,6000,1);
  line.Draw();

  // drawing the lines of the two fits
  TF1 *ratiofit1 = new TF1("fit1_ratio", fitRatio1, 0, 6000);
  TF1 *ratiofit2 = new TF1("fit2_ratio", fitRatio2, 0, 6000);
  ratiofit1->SetLineColor(1);
  ratiofit1->SetLineStyle(2);
  ratiofit1->Draw("same");
  ratiofit2->SetLineColor(1);
  ratiofit2->SetLineStyle(3);
  ratiofit2->Draw("same");

  if (plot_stat_unc) {
    TH1D* uncertainty_fits_ratio = (TH1D*) uncertainty_fits->Clone();
    uncertainty_fits_ratio->Divide(fitMean);
    uncertainty_fits_ratio->Draw("e3 same");
  }
  
  if(plot_other_ratios) c1_hist->SaveAs("plots/backgroundEstimation_HOTVR_" + year + "_" + region + "_" + channel + systematic + JE_string + "_withSysts.pdf");
  else c1_hist->SaveAs("plots/backgroundEstimation_HOTVR_" + year + "_" + region + "_" + channel + systematic + JE_string + ".pdf");

  c1_hist->Clear();

  purity.GetXaxis()->SetTitle("S_{T} [GeV]");
  purity.GetXaxis()->SetNdivisions(505);
  purity.GetYaxis()->SetTitle("purity");
  purity.SetTitle("");

  purity.Draw("AP");

  //text->Draw();
  text2->Draw();
  text3->Draw();

  c1_hist->SaveAs("plots/purity_HOTVR_" + year + "_" + channel + systematic + JE_string + ".pdf");

  // saving fit function to output file
  if(storeOutputToFile) {
    TFile *output;
    TString filename = "files/bgest/alphaFunction_HOTVR_" + year + "_" + region + "_" + channel + systematic + JE_string +".root";
    output = TFile::Open(filename, "RECREATE");
    fit2->Write();
    fit1->Write();
    fitMean->SetName("fit_mean");
    fitMean->Write();
    ratio.SetName("alpha_ratio");
    ratio.Write();
  }

}
