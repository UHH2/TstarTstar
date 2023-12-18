

void backgroundEstimation_compareFunctions () {

  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/bgest/";
  TString filename_base = "alphaFunction_HOTVR";

  std::vector<TString> years = {"UL16preVFP", "UL16postVFP", "UL17", "UL18"};
  //std::vector<TString> years = {""};
  
  //std::vector<TString> channels = {"mu", "ele"};
  std::vector<TString> channels = {"ele"};

  TString region = "VR";
  TString fitfunc = "fit_mean";

  TString mode = "allCombs"; // only used for output naming
  if (years.size() == 1) mode = "channels";
  if (channels.size() == 1) mode = "years";

  TCanvas *canvas = new TCanvas("canvas", "c", 400, 300);

  std::vector<int> colors = {810, 600, 416, 800};
  std::vector<int> styles = {2, 3};

  auto legend = new TLegend(0.2,0.7,0.48,0.85);
  legend->SetBorderSize(0);
  gStyle->SetLegendTextSize(0.04);

  bool draw_baseline = true;
  if (draw_baseline) {
    // draw currently used as baseline
    TFile *input = TFile::Open(path + filename_base + "__" + region + "_total.root");
    TF1 *function = (TF1*)((TF1*) input->Get(fitfunc))->Clone("baseline"); // give some unique name

    function->SetLineColor( 1 );
    function->SetLineStyle( 0 );

    function->GetYaxis()->SetRangeUser(0, 1.3);

    function->SetTitle("");
    function->GetXaxis()->SetTitle("S_{T} [GeV]");
    function->GetYaxis()->SetTitle("#alpha-ratio");

    function->Draw("");
    legend->AddEntry(function, "total", "l");
  }


  bool first = true;
  if (draw_baseline) first = false;
  int yeari = 0;
  for (auto year : years) {

    int channeli = 0;
    for (auto channel : channels) {

      TFile *input = TFile::Open(path + filename_base + "_" + year + "_" + region + "_" + channel + ".root");
      TF1 *function = (TF1*)((TF1*) input->Get(fitfunc))->Clone(year + "_" + region + "_" + channel); // give some unique name

      function->SetLineColor( colors.at(yeari) );
      function->SetLineStyle( styles.at(channeli) );

      if (first) function->Draw("");
      else function->Draw("same");

      legend->AddEntry(function, year + " " + channel, "l");

      first = false;

      channeli++;
    }

    yeari++;
  }

  legend->Draw();

  canvas->SaveAs("plots/bgest_" + mode + "_" + region + ".pdf");

}
