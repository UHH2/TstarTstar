// helper to extract sample norms for all PDF sets after the preselection

void sampleNorm(){

    TString file_dir = "/nfs/dust/cms/user/flabe/TstarTstar/data/Preselection";
    vector<TString> year = {
        "UL16preVFP",
        "UL16postVFP",
        "UL17",
        "UL18"
    };

    vector<TString> masspoints = {"700", "800", "900", "1000", "1100", "1200", "1300", "1400", "1500", "1600", "1700", "1800", "1900", "2000", "2250", "2500", "2750"};

    vector<TString> sample_names;
    for (auto masspoint : masspoints) {
        sample_names.push_back("TstarTstar_M-" + masspoint);
        sample_names.push_back("TstarTstar_Spin32_M-" + masspoint);
    }

    for(unsigned int j=0; j<year.size(); j++){
        cout << "year: " << year.at(j) << endl;

        for(unsigned int i=0; i<sample_names.size(); i++){
            cout << sample_names.at(i) << endl;

            TString filename;
            filename = file_dir + "/" + year.at(j) + "/hadded/uhh2.AnalysisModuleRunner.MC." + sample_names.at(i) + ".root";

            TFile* file_in = new TFile(filename, "READ");

            float nominal = ((TH1F*)(file_in->Get("PDFNorm/nominal")))->GetBinContent(1);
            
            std::ofstream file_out;
            std::ofstream file_out_copy;
            
            file_out.open("files/SignalNorm_" + year.at(j) + "_" + sample_names.at(i) + ".txt");

            for(int i = 1; i<101; i++){
                stringstream ss_name;
                ss_name << "PDFNorm/PDF_" << i;
                string s_name = ss_name.str();
                const char* char_name = s_name.c_str();
                float sum_pdf =  ((TH1F*)(file_in->Get(char_name)))->GetBinContent(1);
                float norm_pdf = nominal / sum_pdf;
                file_out << "pdf" << i << " " << norm_pdf << std::endl;
                file_out_copy << "pdf" << i << " " << norm_pdf << std::endl;
            }

            delete file_in;
        }
    }
}
