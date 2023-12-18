

void checkNonTopFracs() {

    // goal: create plot showing relative contributions of non-top backgrounds in the different regions.
    // will contain three histograms per color (three regions), with line styles representing the regions
    
    // loop over samples
    // loop over regions
    // get all histograms of all samples & yields
    
    // load total histogram for all regions from hadded file (could also construct here, but this is easier!)
    // divide all others to be relative

    // do the plotting

    /////

    // settings & definitions
    std::vector<TString>    samples =   {   "WJets",    "VV",       "QCD",      "DYJets"    };
    std::vector<int>        colors  =   {   600,        416,        867,        400         };

    std::vector<TString>    regions =   {   "SignalRegion",         "ValidationRegion",         "ControlRegion"     };
    std::vector<int>        styles  =   {   1,                      2,                          3,                  };

    TString year                    =   "";           // empty means full UL
    TString channel                 =   "total";

    // constructed variables
    TString path                    =   "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/" + year + "/hadded/";
    TString sample_base             =   "uhh2.AnalysisModuleRunner.MC.";

    TString histname               =   "pt_ST_nominal";


    // load baseline histograms
    TFile *file_baseline = TFile::Open(path + sample_base + "nontop.root");
    if(!file_baseline) std::cout << "Baseline does not exist, please hadd all relevant files to create it!" << std::endl;

    std::vector<TH1D*> hists_baseline;
    for (auto region : regions) {
        TH1D *hist_baseline = (TH1D*)file_baseline->Get(region + "_" + channel + "/" + histname);
        if(!hist_baseline) std::cout << "Baseline hist does not exist for " << region << "." << std::endl;
        hists_baseline.push_back(hist_baseline);

    }

    // plot creation and general style
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

    gROOT->ForceStyle();

    TCanvas *canvas = new TCanvas("canvas", "c", 600, 500);
    auto legend = new TLegend(0.75, 0.62, 0.88, 0.88);

    gPad->SetBottomMargin(0.16);
    gPad->SetLeftMargin(0.16);

    bool first = true;

    // main part, do the loading, division and plotting
    int samplei = 0;
    for (auto sample : samples) {

        TFile *file = TFile::Open(path + sample_base + sample + ".root");
        if(!file) std::cout << "File for sample " << sample << " does not exist." << std::endl;

        int regioni = 0;
        for (auto region : regions) {

            TH1D *hist_baseline = hists_baseline.at(regioni);

            TH1D *hist = (TH1D*)file->Get(region + "_" + channel + "/" + histname);
            if(!hist) std::cout << "Hist does not exist for " << region << " in " << sample << "." << std::endl;

            // dividing by baseline
            hist->Divide(hist_baseline);

            // styling
            hist->SetLineColor( colors.at(samplei) );
            hist->SetLineStyle( styles.at(regioni) );
            hist->SetLineWidth( 2 );
            hist->GetYaxis()->SetRangeUser(0, 1);
            hist->GetXaxis()->SetRangeUser(600, 6000);
            hist->GetXaxis()->SetTitle( hist->GetTitle() );
            hist->SetTitle("");
            hist->GetYaxis()->SetTitle("event fraction");

            // plotting
            if (first) hist->Draw("hist");
            else hist->Draw("hist same");
            first = false;

            if(regioni == 0) legend->AddEntry(hist, sample, "l");

            regioni++;
        }

        samplei++;
    }


    // drawing legends and saving plot
    auto legend2 = new TLegend(0.5, 0.69, 0.72, 0.88);
    int regioni = 0;
    for (auto region : regions) {
        TGraph *dummyLine = new TGraph(); // Empty graph for solid line
        dummyLine->SetLineColor(kGray);
        dummyLine->SetLineStyle( styles.at(regioni) );
        dummyLine->SetLineWidth( 2 );
        legend2->AddEntry(dummyLine, region, "l");
        regioni++;
    }

    legend->SetBorderSize(0);
    legend->Draw();

    legend2->SetBorderSize(0);
    legend2->Draw();

    if (year == "") canvas->SaveAs("plots/nontop_fractions.pdf");
    else canvas->SaveAs("plots/nontop_fractions_" + year + ".pdf");


}