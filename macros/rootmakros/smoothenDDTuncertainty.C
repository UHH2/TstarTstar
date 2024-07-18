

void smoothenDDTuncertainty() {

    TString sample = "TTbar";

    const int nbins = 34;

    // drawing error from decorrelation
    TFile *decorrelationUncertaintyFile = new TFile("/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/decorrelationComparison_MC."+sample+".root");
    TH1* decorrelationUncertainty_ori = (TH1*)decorrelationUncertaintyFile->Get("decorrelation_uncertainty");
    TH1* decorrelationUncertainty = (TH1*)decorrelationUncertainty_ori->Clone();
    for (int bin = -1; bin <= nbins; bin++) {
        decorrelationUncertainty->SetBinContent( bin, abs( decorrelationUncertainty->GetBinContent(bin)-1) );
    }

    // smoothening the histogram
    TH1* decorrelationUncertainty_preSmooth = (TH1*)decorrelationUncertainty->Clone();
    for (int bin = 0; bin <= nbins; bin++) {

        double previousBin = decorrelationUncertainty_preSmooth->GetBinContent(bin-1);
        double currentBin = decorrelationUncertainty_preSmooth->GetBinContent(bin);
        double nextBin = decorrelationUncertainty_preSmooth->GetBinContent(bin+1);

        double average = (previousBin + currentBin + nextBin)/3;
        if(previousBin == 1) average = (currentBin + nextBin)/2;
        else if(nextBin == 1) average = (previousBin + currentBin)/2;

        if(average == 1) average = 0;
        decorrelationUncertainty->SetBinContent(bin, average);
    }

    TFile *outFile = new TFile("/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/decorrelationComparison_"+sample+"_smooth.root", "RECREATE");
    outFile->cd();
    decorrelationUncertainty->GetYaxis()->SetRangeUser(0, 0.5);
    decorrelationUncertainty->Write();
}
