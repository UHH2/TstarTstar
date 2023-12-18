

void backgroundEstimation_comparePurities () {

  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/bgest/";
  TString filename_base = "purity_HOTVR";

  TString mode = "years"; // options are allCombs, channels, years and systs

  std::vector<TString> years;
  if ((mode == "allCombs") || (mode == "years")) years = {"UL16preVFP", "UL16postVFP", "UL17", "UL18"};
  else years = {""};
  
  std::vector<TString> channels;
  if ((mode == "allCombs") || (mode == "channels")) channels = {"mu", "ele"};
  else channels = {"ele"};

  std::vector<TString> systs;
  if ((mode == "allCombs") || (mode == "systs")) systs = {"btagging_totalUp", "btagging_totalDown"};
  else systs = {""};

  TString purityname = "purity";

  TCanvas *canvas = new TCanvas("canvas", "c", 600, 500);

  std::vector<int> colors = {810, 600, 416, 800};
  std::vector<int> styles = {2, 3};

  auto legend = new TLegend(0.6,0.2,0.8,0.4);
  legend->SetBorderSize(0);
  gStyle->SetLegendTextSize(0.04);

  bool draw_baseline = true;
  if (draw_baseline) {
    // draw currently used as baseline
    TFile *input = TFile::Open(path + filename_base + "__ele.root");
    TH1D *function = (TH1D*)((TH1D*) input->Get(purityname))->Clone("baseline"); // give some unique name

    function->SetLineColor( 1 );
    function->SetLineStyle( 0 );

    function->Draw("");
    legend->AddEntry(function, "total", "f");
  }

  bool first = true;
  if (draw_baseline) first = false;
  int yeari = 0;
  for (auto year : years) {

    int channeli = 0;
    for (auto channel : channels) {

      int systi = 0;
      for (auto syst : systs) {

        TFile *input;
        TH1D *function;  // give some unique name

        if (syst == "") {
          input = TFile::Open(path + filename_base + "_" + year + "_" + channel + ".root");
          function = (TH1D*)((TH1D*) input->Get(purityname))->Clone(year + "_" + channel);
        } else {
          input = TFile::Open(path + filename_base + "_" + year + "_" + channel + "_" + syst + ".root");
          function = (TH1D*)((TH1D*) input->Get(purityname))->Clone(year + "_" + channel);
        }
        

        function->SetLineColor( colors.at(yeari + systi) );
        function->SetLineStyle( styles.at(channeli) );

        if (first) function->Draw("");
        else function->Draw("same");

        legend->AddEntry(function, year + " " + channel + " " + syst, "f");

        first = false;

        systi++;

      }

      channeli++;
    }

    yeari++;
  }

  legend->Draw();

  canvas->SaveAs("plots/bgest_purity_" + mode + ".pdf");

}
