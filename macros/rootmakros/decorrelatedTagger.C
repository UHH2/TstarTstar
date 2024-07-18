// Decorrelate DNN tagger, see https://confluence.desy.de/x/0up4DQ
// author: F.Labe
// date: 24.09.2021
// Run it with following command:
// root -l -b -q decorrelatedTagger.C

TF1 *fit1, *fit2, *fitMean;
Double_t meanFunc(Double_t *x, Double_t *par) {
  return ( fit1->EvalPar(x,par) + fit2->EvalPar(x,par) ) / 2;
}

Double_t CrystalBall(Double_t *x, Double_t *par) {
  return ROOT::Math::crystalball_function(x[0], par[0], par[1], par[2], par[3]);
}

void decorrelatedTagger(){

  Double_t efficiency_ttbar = 0.3;
  TString effic_string = "0p3";

  // some style options
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetTitleSize(0.06,"z");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleYOffset(1.20);
  gStyle->SetTitleXOffset(1.0);
  gStyle->SetPalette(55);

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.2);

  // recommended color style
  gStyle->SetPalette(112); // viridis

  gROOT->ForceStyle();
  Double_t w = 600;
  Double_t h = 600;

  TCanvas *c1_hist = new TCanvas("chist", "c", w, h);
  c1_hist->SetLogz();

  // get histogram
  TString filename="TTbar";
  TString subpath="DNN";
  TString histname="2D_DNN_ST";
  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/hadded/";
  TString fileprefix = "uhh2.AnalysisModuleRunner.MC.";
  TFile *input = TFile::Open(path+fileprefix+filename+".root");
  if(!input) cout << "Empty file" << endl;
  TH2D *hist = (TH2D*)input->Get(subpath+"/"+histname);
  if(!hist) cout << "Empty hist" << endl;

  bool store_functions = false;

  // flip
  TH2D *hist2 = new TH2D("oldrebin",hist->GetTitle(), 40, 0, 4000, 50, 0, 1);
  TAxis *xaxis = hist->GetXaxis(); TAxis *yaxis = hist->GetYaxis();
  for (int j=1; j<=yaxis->GetNbins();j++) {
    for (int i=1; i<=xaxis->GetNbins();i++) {
      hist2->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(yaxis->GetNbins()-j),hist->GetBinContent(i,j));
    }
  }
  hist = hist2;

  // calculate DNN cut value for each bin in x
  int binsX = 99;
  int binsY = 60;
  std::vector<double> cutedgesY;
  std::vector<double> errors;
  std::vector<double> errors_down;
  std::vector<double> errors_up;
  std::vector<double> cutedgesX;
  std::vector<double> binwidthX;
  for(uint x = 2; x <= binsX; x++){
    double total = hist->Integral(x, x, 0, binsY);
    Double_t partial = 0;
    Double_t value;
    // nominal value
    for(uint y = 0; y < binsY; y++){
      partial = hist->Integral(x, x, 0, y);
      if(partial/total > efficiency_ttbar) {
        value = hist->GetYaxis()->GetBinCenter(y);
        cutedgesY.push_back(value);
        errors.push_back((sqrt(total)/total));
        cutedgesX.push_back(hist->GetXaxis()->GetBinCenter(x));
        binwidthX.push_back(hist->GetXaxis()->GetBinWidth(x));
        break;
      }
    }

    // up
    for(uint y = binsY; y > 0; y--){
      Double_t error;
      partial = hist->IntegralAndError(x, x, y, binsY, error);
      error = sqrt(partial)/partial;
      partial+=error;
      if(partial/total > efficiency_ttbar) {
        double result = hist->GetYaxis()->GetBinCenter(y)-value;
        if((hist->GetXaxis()->GetBinCenter(x) > 3000) && result == 0) result = 0.02;
        errors_up.push_back(result);
        break;
      }
    }

    // down
    for(uint y = binsY; y > 0; y--){
      Double_t error;
      partial = hist->IntegralAndError(x, x, y, binsY, error);
      error = sqrt(partial)/partial;
      partial-=error;
      if(partial/total > efficiency_ttbar) {
        double result = value-hist->GetYaxis()->GetBinCenter(y);
        errors_down.push_back(result);
        break;
      }
    }

  }

  // draw plot
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  
  TH1* frame = gPad->DrawFrame(500, 0, 4000, 1.3);
  frame->GetXaxis()->SetTitle("S_{T} [GeV]");
  frame->GetYaxis()->SetTitle("1 #minus s_{DNN}");
  frame->GetXaxis()->SetNdivisions(505);

  hist->SetTitle("");
  hist->GetZaxis()->SetTitle("Events / bin");
  hist->GetZaxis()->SetTitleOffset(1.15);
  hist->GetZaxis()->SetRangeUser(0.1, 1e4);

  hist->Draw("colz same");

  // draw Lumi text
  TString infotext = TString::Format("%3.0f fb^{-1} (%d TeV)", 137.6, 13);
  TLatex *text = new TLatex(3.5, 24, infotext);
  text->SetNDC();
  text->SetTextAlign(33);
  text->SetX(0.85);
  text->SetTextFont(42);
  text->SetY(0.9872);
  text->SetTextSize(0.035);
  text->Draw();

  double CMS_size = 0.06;
  double ratio = 1.3;
  // draw CMS Work in Progress text
  TString cmstext = "CMS";
  TLatex *text2 = new TLatex(3.5, 24, cmstext);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.225);
  text2->SetTextFont(62);
  text2->SetTextSize(CMS_size);
  text2->SetY(0.91);
  text2->Draw();

  TString preltext = "Simulation";
  TLatex *text3 = new TLatex(3.5, 24, preltext);
  text3->SetNDC();
  text3->SetTextAlign(13);
  text3->SetX(0.225);
  text3->SetTextFont(52);
  text3->SetTextSize(CMS_size / ratio);
  text3->SetY(0.85);
  text3->Draw();

  TGraph* edgePoints = new TGraphAsymmErrors(cutedgesX.size(),&cutedgesX[0],&cutedgesY[0], 0, 0, &errors[0], &errors[0]);
  edgePoints->SetMarkerStyle(1);
  //edgePoints->Draw("p same");
  //other option: Gauß + exponential
  //fit1 = new TF1("fit1", "[2] + [0] * exp([1] * x)", 500, 4000);
  fit1 = new TF1("fit1", "crystalball", 500, 4000);
  //fit2 = new TF1("fit2", "[2] + [0] * exp([1] * x)", 500, 4000);
  fit2 = new TF1("fit2", "crystalball", 500, 4000);
  // constant, mean, sigma, alpha, N
  fit1->SetParameters(.8,700,400,-0.2,8e+05);
  fit2->SetParameters(.8,700,400,-0.2,8e+05);
  //fit->SetParameters(2, -3.05300e-03, 4.17680e-01);
  edgePoints->Fit("fit1", "N", "", 500, 4000);
  edgePoints->Fit("fit2", "N", "", 500, 4000);

  fitMean = new TF1("mean", meanFunc, 0, 6000);

  fitMean->SetLineWidth(3);

  fitMean->Draw("same");

  // plot some second line
  /**
  double crystal_constant = 4.65951e-01;
  double crystal_mean = 7.13578e+02;
  double crystal_sigma = 2.40880e+02;
  double crystal_alpha = -1.43648e-01;
  double crystal_n = 5.14847e+05;
  TF1 *otherfunc = new TF1("otherfunc", "crystalball", 500, 4000);
  otherfunc->SetParameters(crystal_constant,crystal_mean,crystal_sigma,crystal_alpha,crystal_n);
  otherfunc->SetLineColor(kBlue);
  otherfunc->Draw("L same");
  **/

  // add legend
  auto legend = new TLegend(0.52, 0.71, 0.77, 0.81);
  gStyle->SetLegendTextSize(0.045);
  TString fitlabelnumber = TString::Format("%.0f", efficiency_ttbar * 100);
  legend->AddEntry(fitMean, "f(S_{T}, "+fitlabelnumber+"%)", "l");

  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->Draw();

  gPad->RedrawAxis();

  c1_hist->SaveAs("plots/variableCuts_"+filename+"_"+subpath+"_"+histname+"_"+effic_string+".pdf");

  // saving fit function to output file
  if(store_functions) {
    TFile *output;
    output = TFile::Open("files/DDTfunc_"+effic_string+".root", "RECREATE");
    fitMean->Write();
  }

}
