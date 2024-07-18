
void btagYieldUnc() {

    TString year = "UL16preVFP";
    TString channel = "mu";
    TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/"+year+"/hadded/";
    TString filename_base = "uhh2.AnalysisModuleRunner.MC.";

    std::vector<TString> samples = {"TTbar", "WJets", "ST", "QCD", "VV", "DYJets", "TstarTstar", "TstarTstar_Spin32"};
    //std::vector<TString> samples = {"TstarTstar", "TstarTstar_Spin32"};

    TString histname = "pt_"+channel;

    TString folder_before = "HOTVRcut_"+channel;
    //TString folder_after = "bcorrections_"+channel;
    TString folder_crosscheck = "byield_"+channel;

    bool writeSFsToFile = true;

    std::vector<TH1D*> histograms_to_store;
    for (const auto sample : samples) {
        std::cout << "Processing " << sample << "..." << std::endl;

        TFile *input = TFile::Open(path+filename_base+sample+".root");
        if(!input) cout << "Empty file" << endl;
        TH1D *hist_before = (TH1D*)input->Get(folder_before+"/"+histname);
        if(!hist_before) cout << "Empty hist before" << endl;
        TH1D *hist_crosscheck = (TH1D*)input->Get(folder_crosscheck+"/"+histname);
        if(!hist_crosscheck) cout << "Empty hist crosscheck" << endl;

        // clone it, divide, then close the file
        TH1D *hist_ratio = (TH1D*)hist_crosscheck->Clone();
        hist_ratio->Divide(hist_before);

        hist_ratio->SetName(sample);

        // this is the scale factor, now lets store it!
        histograms_to_store.push_back(hist_ratio);

    }

    if(writeSFsToFile) {
        // outputting the 2D histograms for scaling
        TFile *output = TFile::Open("files/btagyield/btagYieldUncs_"+year+"_"+channel+".root", "RECREATE");
        for (auto histogram : histograms_to_store) {
            histogram->Write();
        }
    }

}