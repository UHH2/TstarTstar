
// function to propagate normalization uncertainties to the alpha ratio, to show that this is small :)

// helper functions to square and root a histogram
void squareHistogram(TH1D * hist) {
    for (int bin=0; bin<=hist->GetNcells(); ++bin) {
        double content = hist->GetBinContent(bin);
        hist->SetBinContent(bin, content*content);
    }
}

void rootHistogram(TH1D * hist) {
    for (int bin=0; bin<=hist->GetNcells(); ++bin) {
        double content = hist->GetBinContent(bin);
        hist->SetBinContent(bin, sqrt(content));
    }
}

void backgroundEstimation_normUncPropagation() {

    

    // relative uncertainties on the xsec in all nontop samples
    std::vector<double> nontop_uncertainties = {0.2, 0.2, 0.2, 0.2};

    // #############################################################################
    // first, we'll calculate the alpha ratio in the same way as in the main file. also some plotting.
    // #############################################################################

    ///////// PLOTTING //////////
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
    gStyle->SetPadBottomMargin(0.2);
    gStyle->SetPadLeftMargin(0.18);
    gStyle->SetPadRightMargin(0.1);

    gROOT->ForceStyle();
    Double_t w = 600;
    Double_t h = 600;

    TCanvas *canvas = new TCanvas("chist", "c", w, h);

    auto legend = new TLegend(0.22,0.68,.7,0.94);
    legend->SetBorderSize(0);
    gStyle->SetLegendTextSize(0.04);

    TString region = "VR";
    TString year = "";              // no year means full run 2
    TString channel = "total";      // total channel means combination of both

    // definitions
    std::vector<TString> nontop_backgrounds = {"WJets",  "VV", "QCD", "DYJets"};
    std::vector<TString> top_backgrounds = {"ST", "TTbar"};

    std::vector<int>        colors  =   {   600,        416,        867,        400         };

    TString subpath_SR;
    if (region == "VR") subpath_SR = "ValidationRegion_" + channel;
    else if (region == "SR") subpath_SR = "SignalRegion_" + channel;
    TString subpath_CR = "ControlRegion_" + channel;  
    TString histname = "pt_ST_nominal";
    TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/";
    TString fileprefix = "uhh2.AnalysisModuleRunner.";

    path = path + "/" + year + "/hadded/";

    std::cout << "Using path: " << path << std::endl;

    // first get the baseline
    TH1D *histSR_nontop;
    TH1D *hist_btagCR_nontop;

    // open full set of non-top backgrounds in signal region
    int nontopi = 0;
    for (auto filename : nontop_backgrounds) {

        TFile *input = TFile::Open(path+fileprefix+"MC."+filename+".root");
        if(!input) cout << "Empty file" << endl;
        TH1D *hist = (TH1D*)input->Get(subpath_SR+"/"+histname);
        if(!hist) cout << "Empty hist" << endl;

        if (nontopi == 0) {
            histSR_nontop = hist;
        } else {
            histSR_nontop->Add(hist);
        }

        nontopi++;
    }

    // open full set of non-top backgrounds in W-jets control region
    nontopi = 0;
    for (auto filename : nontop_backgrounds) {
        TFile *input = TFile::Open(path+fileprefix+"MC."+filename+".root");
        if(!input) cout << "Empty file" << endl;
        TH1D *hist = (TH1D*)input->Get(subpath_CR+"/"+histname);
        if(!hist) cout << "Empty hist" << endl;

        // outputting some information on the b-tagging CR
        if (nontopi == 0) {
            hist_btagCR_nontop = hist;
        } else {
            hist_btagCR_nontop->Add(hist);
        }

        nontopi++;
    }

    // calculate ratio histogram
    TH1D* ratio = (TH1D*)histSR_nontop->Clone();
    ratio->Divide(hist_btagCR_nontop);

    // drawing ratio
    ratio->GetXaxis()->SetTitle("S_{T} [GeV]");
    ratio->GetXaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitle("ratio");
    ratio->GetYaxis()->SetRangeUser(0, 2);
    ratio->GetXaxis()->SetRangeUser(0, 6000);
    if(region == "SR") ratio->GetYaxis()->SetRangeUser(0, 0.5);
    ratio->SetTitle("");
    ratio->SetLineColor(1);
    ratio->SetMarkerColor(1);
    ratio->SetMarkerStyle(20);
    ratio->Draw("");

    legend->AddEntry(ratio,"#alpha (nominal)","elp");


    //////////// SOURCE LÖÖÖÖP ////////////

    // then looping over sources, as they are independent!
    for (int uncertainti = 0; uncertainti < nontop_uncertainties.size(); uncertainti++) {
        
        // up variation
        TH1D *histSR_nontop_up;
        TH1D *hist_btagCR_nontop_up;

        // open full set of non-top backgrounds in signal region
        nontopi = 0;
        for (auto filename : nontop_backgrounds) {

            TFile *input = TFile::Open(path+fileprefix+"MC."+filename+".root");
            if(!input) cout << "Empty file" << endl;
            TH1D *hist = (TH1D*)input->Get(subpath_SR+"/"+histname);
            if(!hist) cout << "Empty hist" << endl;

            // scale for up- and down-variation
            if (nontopi == uncertainti) hist->Scale(1 + nontop_uncertainties.at(nontopi));

            if (nontopi == 0) {
                histSR_nontop_up = hist;
            } else {
                histSR_nontop_up->Add(hist);
            }

            nontopi++;
        }

        // open full set of non-top backgrounds in W-jets control region
        nontopi = 0;
        for (auto filename : nontop_backgrounds) {
            TFile *input = TFile::Open(path+fileprefix+"MC."+filename+".root");
            if(!input) cout << "Empty file" << endl;
            TH1D *hist = (TH1D*)input->Get(subpath_CR+"/"+histname);
            if(!hist) cout << "Empty hist" << endl;

            // scale for up- and down-variation
            if (nontopi == uncertainti) hist->Scale(1 + nontop_uncertainties.at(nontopi));

            // outputting some information on the b-tagging CR
            if (nontopi == 0) {
                hist_btagCR_nontop_up = hist;
            } else {
                hist_btagCR_nontop_up->Add(hist);
            }

            nontopi++;
        }

        // calculate ratio histogram
        TH1D* ratio_up = (TH1D*)histSR_nontop_up->Clone();
        ratio_up->Divide(hist_btagCR_nontop_up);

        // down variation
        TH1D *histSR_nontop_down;
        TH1D *hist_btagCR_nontop_down;

        // open full set of non-top backgrounds in signal region
        nontopi = 0;
        for (auto filename : nontop_backgrounds) {

            TFile *input = TFile::Open(path+fileprefix+"MC."+filename+".root");
            if(!input) cout << "Empty file" << endl;
            TH1D *hist = (TH1D*)input->Get(subpath_SR+"/"+histname);
            if(!hist) cout << "Empty hist" << endl;

            // scale for down- and down-variation
            if (nontopi == uncertainti) hist->Scale(1 - nontop_uncertainties.at(nontopi));

            if (nontopi == 0) {
                histSR_nontop_down = hist;
            } else {
                histSR_nontop_down->Add(hist);
            }

            nontopi++;
        }

        // open full set of non-top backgrounds in W-jets control region
        nontopi = 0;
        for (auto filename : nontop_backgrounds) {
            TFile *input = TFile::Open(path+fileprefix+"MC."+filename+".root");
            if(!input) cout << "Empty file" << endl;
            TH1D *hist = (TH1D*)input->Get(subpath_CR+"/"+histname);
            if(!hist) cout << "Empty hist" << endl;

            // scale for down- and down-variation
            if (nontopi == uncertainti) hist->Scale(1 - nontop_uncertainties.at(nontopi));

            // outputting some information on the b-tagging CR
            if (nontopi == 0) {
                hist_btagCR_nontop_down = hist;
            } else {
                hist_btagCR_nontop_down->Add(hist);
            }

            nontopi++;
        }

        // calculate ratio histogram
        TH1D* ratio_down = (TH1D*)histSR_nontop_down->Clone();
        ratio_down->Divide(hist_btagCR_nontop_down);

        legend->AddEntry(ratio_down, nontop_backgrounds.at(uncertainti) + " #sigma #pm 20%", "l");

        // drawing errors
        ratio_up->SetLineColor(colors.at(uncertainti));
        ratio_up->SetMarkerColor(colors.at(uncertainti));
        ratio_up->SetLineWidth(2);
        //ratio_up->SetMarkerStyle(20);
        ratio_up->Draw("same");
        ratio_down->SetLineColor(colors.at(uncertainti));
        ratio_down->SetMarkerColor(colors.at(uncertainti));
        ratio_down->SetLineWidth(2);
        //ratio_down->SetMarkerStyle(20);
        ratio_down->Draw("same");

    }

    legend->Draw();

    canvas->SaveAs("plots/norm_uncertainty_propagation_" + region + "_" + year + "_" + channel + ".pdf");


}
