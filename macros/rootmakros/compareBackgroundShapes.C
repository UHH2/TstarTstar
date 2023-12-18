

// helper function to compare sample shapes between years

void compareBackgroundShapes() {

    // definitions / settings
    std::vector<TString>    samples =   {   "TTbar",    "WJets",    "QCD",      "DYJets"    }; // VV dropped for the moment, stats too low
    std::vector<int>        colors  =   {   810,        600,         867,        400         };

    std::vector<TString>    years   =   {   "UL16preVFP",   "UL16postVFP",  "UL17",     "UL18"      };
    std::vector<int>        styles  =   {   1,              2,              3,          4           };

    TString channel = "mu";
    TString region = "ControlRegion";

    TString histname = "pt_ST_nominal";
    TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/";
    TString filebase = "uhh2.AnalysisModuleRunner.MC.";

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
    auto legend = new TLegend(0.75, 0.6, 0.88, 0.88);
    bool first = true;

    gPad->SetLogy();
    gPad->SetBottomMargin(0.16);
    gPad->SetLeftMargin(0.16);

    // main lööp
    int samplei = 0;
    for (auto sample : samples) {

        int yeari = 0;
        for (auto year : years) {

            // open file
            TFile *file = TFile::Open(path + "/" + year + "/hadded/" + filebase + sample + ".root");
            if(!file) std::cout << "File for sample " << sample << " in year " + year + " does not exist." << std::endl;

            // get histogram
            TH1D *hist = (TH1D*)file->Get(region + "_" + channel + "/" + histname);
            if(!hist) std::cout << "Hist does not exist for " << sample << " in year " << year << "." << std::endl;

            // normalizing
            hist->Scale( 1/hist->Integral() );

            // styling
            hist->SetLineColor( colors.at(samplei) );
            hist->SetLineStyle( styles.at(yeari) );
            hist->SetLineWidth( 2 );
            hist->GetXaxis()->SetRangeUser(600, 6000);
            hist->GetXaxis()->SetTitle( hist->GetTitle() );
            hist->SetTitle("");
            hist->GetYaxis()->SetTitle("events");

            // plotting
            if (first) hist->Draw("hist");
            else hist->Draw("hist same");
            first = false;

            if(yeari == 0) legend->AddEntry(hist, sample, "l");

            yeari++;
        }

        samplei++;
    }

    // drawing legends and saving plot
    auto legend2 = new TLegend(0.5, 0.67, 0.72, 0.88);
    int yeari = 0;
    for (auto year : years) {
        TGraph *dummyLine = new TGraph(); // Empty graph for solid line
        dummyLine->SetLineColor(kGray);
        dummyLine->SetLineStyle( styles.at(yeari) );
        dummyLine->SetLineWidth( 2 );
        legend2->AddEntry(dummyLine, year, "l");
        yeari++;
    }

    legend->SetBorderSize(0);
    legend->Draw();

    legend2->SetBorderSize(0);
    legend2->Draw();

    canvas->SaveAs("plots/backgroundShapes_" + region + "_" + channel + ".pdf");

}
