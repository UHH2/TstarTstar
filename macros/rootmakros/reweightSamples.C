// A short script to prepare the Root files for the DNN training
// author: F.Labe
// date: 30.06.2010
// Run it with following command:
// root -l -b -q reweightSamples.C

void reweightSamples(){

  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/data/Selection/UL16postVFP/";
  TString backgroundName = "hadded/uhh2.AnalysisModuleRunner.MC.TTbar.root";
  TString signalName = "hadded/uhh2.AnalysisModuleRunner.MC.TstarTstar_M-";
  std::vector<TString> masspoints = {"700", "800", "900", "1000", "1100", "1300", "1400", "1500", "1600", "1700", "1800", "1900", "2000", "2250", "2500", "2750"};
  TString add_masspoint = "1200";

  TCanvas *canvas = new TCanvas("chist", "c", 800, 600);

  TString plotdir = "plots/reweighting/";
  TString datadir = "files/";

  // ##########################################
  // ## Calculating ST ratio for reweighting ##
  // ##########################################

  std::cout << "Start getting ST ratio histogram." << endl;

  // get Signal ST hist
  TH1D *ST_sig;
  {
    TFile *file = new TFile(path+signalName+add_masspoint+".root", "READ");
    ST_sig = (TH1D*)file->Get("Aftertopptreweighting/pt_ST");
  }
  for (const auto & masspoint : masspoints){
    TFile *file = new TFile(path+signalName+masspoint+".root", "READ");
    ST_sig->Add((TH1D*)file->Get("Aftertopptreweighting/pt_ST"));
  }
  ST_sig->SetLineWidth(4);
  ST_sig->Draw();
  canvas->SaveAs(plotdir+"ST_sig.pdf");

  // get BG ST hist
  TH1D *ST_bkg;
  {
    TFile *file = new TFile(path+backgroundName, "READ");
    ST_bkg = (TH1D*)file->Get("Aftertopptreweighting/pt_ST");
  }
  ST_bkg->SetLineWidth(4);
  ST_bkg->Draw();
  canvas->SaveAs(plotdir+"ST_bkg.pdf");

  // calculate ratio hist
  // Now have to divide Background weights by this value
  TH1D *ST_ratio = ST_bkg;
  ST_ratio->Divide(ST_sig);

  // Printing Hist
  canvas->SetLogy();
  ST_ratio->SetLineWidth(4);
  ST_ratio->Draw();
  canvas->SaveAs(plotdir+"ST_ratio.pdf");

  //canvas->Close();
  //gSystem->ProcessEvents();

  double sig = ST_sig->Integral();
  double bkg = ST_bkg->Integral();
  std::cout << "Ratio sig/bkg: " << sig/bkg << std::endl;

  TH1D *ST_sig_ori = new TH1D(*ST_sig);
  TH1D *ST_sig_ori_2 = new TH1D(*ST_sig);
  TH1D *ST_bkg_ori = new TH1D(*ST_bkg);
  TH1D *ST_bkg_ori_2 = new TH1D(*ST_bkg);

  // Saving the hists.
  TFile *file = new TFile(datadir+"ST_weights.root", "RECREATE");
  ST_ratio->SetName("ST_ratio");
  ST_ratio->Write();

  // split into parts
  std::vector<TH1D*> ST_sig_vec;
  std::vector<TH1D*> ST_bkg_vec;
  for (uint i = 0; i < 6; i++){
    TString numberstring = std::to_string(i);
    ST_sig_vec.push_back(new TH1D("ST_split_sig_"+numberstring, "ST_split_sig_"+numberstring, 30, 0, 3000));
    ST_bkg_vec.push_back(new TH1D("ST_split_bkg_"+numberstring, "ST_split_bkg_"+numberstring, 30, 0, 3000));
    for(uint j = 1; j <= 5; j++){
      ST_sig_vec.at(i)->SetBinContent((i*5)+j, ST_sig->GetBinContent((i*5)+j));
      ST_bkg_vec.at(i)->SetBinContent((i*5)+j, ST_bkg->GetBinContent((i*5)+j));
    }
  }

  // plot the splitting for check
  ST_sig_vec.at(0)->SetLineColor(1);
  ST_sig_vec.at(0)->SetMaximum(1e4);
  ST_sig_vec.at(0)->Draw();
  for(uint i = 1; i < 6; i++){
    ST_sig_vec.at(i)->SetLineColor(i+1);
    ST_sig_vec.at(i)->Draw("same");
  }
  canvas->SaveAs(plotdir+"ST_sig_split.pdf");
  ST_bkg_vec.at(0)->SetLineColor(1);
  ST_bkg_vec.at(0)->SetMaximum(1e4);
  ST_bkg_vec.at(0)->SetMinimum(1e-3);
  ST_bkg_vec.at(0)->Draw();
  for(uint i = 1; i < 6; i++){
    ST_bkg_vec.at(i)->SetLineColor(i+1);
    ST_bkg_vec.at(i)->Draw("same");
  }
  canvas->SaveAs(plotdir+"ST_bkg_split.pdf");

  // scale split stuff
  for (uint i = 0; i < 6; i++){
    ST_sig_vec.at(i)->Scale(1/(ST_sig_vec.at(i)->Integral()));
    ST_bkg_vec.at(i)->Scale(1/(ST_bkg_vec.at(i)->Integral()));
  }

  TH1D *ST_sig_split = new TH1D("ST_sig_split", "ST_sig_split", 30, 0, 3000);
  TH1D *ST_bkg_split = new TH1D("ST_bkg_split", "ST_bkg_split", 30, 0, 3000);
  for (uint i = 0; i < 6; i++){
    for(uint j = 1; j <= 5; j++){
      ST_sig_split->SetBinContent((i*5)+j, ST_sig_vec.at(i)->GetBinContent((i*5)+j));
      ST_bkg_split->SetBinContent((i*5)+j, ST_bkg_vec.at(i)->GetBinContent((i*5)+j));
    }
  }

  ST_sig_split->Write();
  ST_bkg_split->Write();

  // scale and save total hists
  ST_sig->Scale(1/(ST_sig->Integral()));
  ST_sig->SetName("ST_sig");
  ST_sig->Write();
  ST_bkg->Scale(1/(ST_bkg->Integral()));
  ST_bkg->SetName("ST_bkg");
  ST_bkg->Write();

  ST_sig_ori->Divide(ST_sig);
  ST_bkg_ori->Divide(ST_bkg);

  ST_sig_ori->SetLineWidth(4);
  ST_sig_ori->Draw();
  canvas->SaveAs(plotdir+"ST_sig_scaled.pdf");

  ST_bkg_ori->SetLineWidth(4);
  ST_bkg_ori->Draw();
  canvas->SaveAs(plotdir+"ST_bkg_scaled.pdf");

  TH1D *ST_sig_rebinned = (TH1D*) ST_sig->Clone();
  TH1D *ST_bkg_rebinned = (TH1D*) ST_bkg->Clone();
  for(uint i = 0; i < 6; i++){
    double sig_val = ST_sig_rebinned->GetBinContent((i*5)+3);
    double bkg_val = ST_bkg_rebinned->GetBinContent((i*5)+3);
    for(uint j = 1; j <= 5; j++){
      ST_sig_rebinned->SetBinContent((i*5)+j, sig_val);
      ST_bkg_rebinned->SetBinContent((i*5)+j, bkg_val);
    }
  }

  ST_sig_rebinned->SetName("ST_sig_rebinned");
  ST_bkg_rebinned->SetName("ST_bkg_rebinned");

  ST_sig_rebinned->SetLineWidth(4);
  ST_sig_rebinned->Draw();
  canvas->SaveAs(plotdir+"ST_sig_rebinned.pdf");
  ST_bkg_rebinned->SetLineWidth(4);
  ST_bkg_rebinned->Draw();
  canvas->SaveAs(plotdir+"ST_bkg_rebinned.pdf");

  ST_sig_rebinned->Scale(1/ST_sig_rebinned->Integral());
  ST_bkg_rebinned->Scale(1/ST_bkg_rebinned->Integral());

  ST_sig_rebinned->Write();
  ST_bkg_rebinned->Write();

  file->Close();
  delete file;

}
