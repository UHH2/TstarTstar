// Just a basic macro to create simple plots from given hists.
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2010
// Run it with following command:
// root -l -b -q FitModelLine.C

void FitModelLine(){
  
  Double_t x[] = {700, 800, 900, 1000, 1200, 1300, 1400, 1500, 1600};
  Double_t x_2[] = {700, 800, 900, 1000, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
  Double_t y_theory[] = {4.92, 1.682, 6.36e-01, 2.62e-01, 5.37e-02, 2.61e-02, 1.31e-02, 6.77e-03, 3.59e-03};  
  
  Double_t y_samples[] = {0.002122, 0.0006901, 0.0002398, 0.00006587, 0.00001477, 6.329e-06, 2.813e-06, 1.271e-06, 5.863e-07, 2.755e-07};
  Double_t y_samples_more[] = {0.002122, 0.0006901, 0.0002398, 0.00006587, 0.00001477, 6.329e-06, 2.813e-06, 1.271e-06, 5.863e-07, 2.755e-07, 1.259e-07, 5.858e-08, 2.499e-08};

  assert(x.size() == y_sampels.size());
  assert(x_2.size() == y_samples_more.size());

  // calculate ratio
  Double_t y_ratio[9] = {0};
  for(uint i = 0; i < 9; i++){
    y_ratio[i] = y_samples[i]/y_theory[i];
  }

  TCanvas *c1 = new TCanvas("chist", "c", 800, 600);
  TGraph *graph = new TGraph(10, x, y_ratio);

  TF1 *f1 = new TF1("f1", "pol 1", 700, 1600);
  graph->Fit("f1", "N");

  graph->SetTitle("");
  graph->SetMarkerStyle(104);
  graph->Draw("AP");  
  f1->Draw("same");

  //c1->SaveAs("Plot_Fit_Ratio.pdf");

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
  TLegend *legend = new TLegend(0.72,0.72,0.92,0.92);
  legend->SetBorderSize(0);
  TGraph *graph_samples = new TGraph(14, x_2, y_samples_more);  
  TGraph *graph_theory = new TGraph(10, x, y_theory);
  TGraph *graph_theory_new = new TGraph(14, x_2, y_theory_new);
  
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.16);
  gPad->SetLeftMargin(0.16);

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

  graph_theory->SetMarkerSize(2);
  graph_samples->SetMarkerSize(2);
  graph_theory_new->SetMarkerSize(2);

  graph_theory->SetMarkerColor(2);
  graph_samples->SetMarkerColor(kGreen+2);
  graph_theory_new->SetMarkerColor(4);

  graph_theory_new->SetMaximum(20);
  graph_theory_new->SetMinimum(1e-8);

  graph_theory_new->SetTitle("");
  graph_theory_new->GetXaxis()->SetTitle("m_{T*} [GeV]");
  graph_theory_new->GetYaxis()->SetTitle("#sigma [pb]");

  graph_theory_new->GetXaxis()->SetTitleFont(42);   
  graph_theory_new->GetXaxis()->SetLabelFont(42);           
  graph_theory_new->GetYaxis()->SetTitleFont(42);            
  graph_theory_new->GetYaxis()->SetLabelFont(42);                          
                                                                                                     
  graph_theory_new->GetXaxis()->SetLabelSize(0.055); // 0.045                                     
  graph_theory_new->GetXaxis()->SetLabelOffset(0.01);                                  
  graph_theory_new->GetXaxis()->SetTickLength(0.03);                                                     
  graph_theory_new->GetXaxis()->SetTitleSize(0.05);                                                             
  graph_theory_new->GetXaxis()->SetTitleOffset(1.2);                                                    
                                                                                                             
  graph_theory_new->GetYaxis()->SetTitleOffset(1.2); // 1.8                                       
  graph_theory_new->GetYaxis()->SetTitleSize(0.06); // 0.05            
  graph_theory_new->GetYaxis()->SetLabelSize(0.05); // 0.045             
  graph_theory_new->GetYaxis()->SetTickLength(0.02); 
  graph_theory_new->GetYaxis()->SetLabelOffset(0.011);    

  graph_theory_new->Draw("AP");
  graph_theory_new->GetXaxis()->SetRangeUser(600,2100);
  graph_theory_new->Draw("AP");
  c->Update();
  graph_theory->Draw("P SAME");
  graph_samples->Draw("P SAME");
  //hist->Draw("same");
  f2->SetLineColor(1);
  //f2->Draw("same");

  graph_theory_new->GetXaxis()->SetRangeUser(600,2100);
  
  legend->AddEntry(graph_theory, "theory", "P");
  legend->AddEntry(graph_samples, "samples", "P");
  legend->AddEntry(graph_theory_new, "theory fit", "P");
  legend->Draw();  

  c->SaveAs("Plot_Fit_Models.pdf");

  // output values into a text file
  ofstream myfile;
  myfile.open("model_values.txt");
  myfile << "[";
  for(uint i = 0; i < 13; i++){
    if(i+7==11) continue;
    myfile << graph_theory_new->Eval((i+7)*100);
    myfile << ", ";
  }
  myfile << f2->Eval(2000) << "]" << endl;
  
}
