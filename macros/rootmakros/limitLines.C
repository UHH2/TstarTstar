// Just a basic macro to create simple plots from given hists.
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2019
// Run it with following command:
// root -l -b -q FitModelLine.C

void FitModelLine(){
  
  // define masspoints
  std::vector<int> masspoints;
  
  // TODO get input values from some file
  std::vector<std::vector<double>> limitValues;
  std::vector<TString> limitLabels;
  
  assert(limitLabels.size() == limitValues.size());

  for(uint i = 0; i < limitValues.size(); i++) {
    // create graph from values

    // plot graph
    // 1st is theory line
    // 2nd is "main" limit


  }

  TCanvas *c1 = new TCanvas("chist", "c", 800, 600);
  TGraph *graph = new TGraph(10, x, y_ratio);

  TF1 *f1 = new TF1("f1", "pol 1", 700, 1600);
  graph->Fit("f1", "N");

  graph->SetTitle("Sample / Theory ratio");
  graph->SetMarkerStyle(2);
  graph->Draw("AP");  
  f1->Draw("same");

  c1->SaveAs("Plot_Fit_Ratio.pdf");

  //calculate new values from ratio
  Double_t y_theory_new[14] = {0};
  for (uint i = 0; i < 10; i++){
    y_theory_new[i] = y_samples[i]/f1->Eval((i+6)*100);
  }
  for (uint i = 10; i < 14; i++){
    y_theory_new[i] = y_samples_more[i]/f1->Eval((i+6)*100);
  }

  // plotting stuff
  TCanvas *c = new TCanvas("c", "c", 800, 600);
  TLegend *legend = new TLegend(0.2,0.2);
  TGraph *graph_theory = new TGraph(10, x, y_theory);
  TGraph *graph_samples = new TGraph(14, x_2, y_samples_more);
  TGraph *graph_theory_new = new TGraph(14, x_2, y_theory_new);
  
  // fit new line
  TH1D *hist = new TH1D("hist", "hist", 14, 650, 2050);
  for (uint i = 0; i < 14; i++){
    hist->Fill(x_2[i], y_theory_new[i]);
  }

  TF1 *f2 = new TF1("exp", "expo", 650, 2050);
  hist->Fit("exp", "LMI", "", 900, 2000);
  
  c->SetLogy();

  graph_theory->SetMarkerStyle(2);
  graph_samples->SetMarkerStyle(2);
  graph_theory_new->SetMarkerStyle(2);

  graph_theory->SetMarkerColor(2);
  graph_samples->SetMarkerColor(3);
  graph_theory_new->SetMarkerColor(4);

  graph_theory_new->SetMaximum(10);
  graph_theory_new->SetMinimum(1e-7);

  graph_theory_new->SetTitle("Model lines & fit");

  graph_theory_new->Draw("AP");
  graph_theory->Draw("P SAME");
  graph_samples->Draw("P SAME");
  //hist->Draw("same");
  f2->SetLineColor(1);
  f2->Draw("same");
  
  legend->AddEntry(graph_theory, "theory");
  legend->AddEntry(graph_samples, "samples");
  legend->AddEntry(graph_theory_new, "theory fit");
  legend->Draw();  

  c->SaveAs("Plot_Fit_Models.pdf");

  TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
  TLegend *legend2 = new TLegend(0.2,0.2);

  c3->SetLogy();

  f2->SetMaximum(10);
  f2->SetMinimum(1e-4);

  f2->SetLineColor(1);
  graph_theory->SetLineColor(2);

  f2->SetTitle("Comparison of old and new model line");

  f2->Draw();
  graph_theory->Draw("L SAME");

  legend2->AddEntry(f2, "new");
  legend2->AddEntry(graph_theory, "previous");
  legend2->Draw();

  c3->SaveAs("Plot_Models.pdf");

  // output values into a text file
  ofstream myfile;
  myfile.open("model_values.txt");
  myfile << "[";
  for(uint i = 0; i < 13; i++){
    myfile << graph_theory_new->Eval((i+7)*100);
    myfile << ", ";
  }
  myfile << f2->Eval(2000) << "]" << endl;
  
}
