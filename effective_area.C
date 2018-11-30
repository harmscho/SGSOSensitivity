TFile f("performance.root");
TGraph* HistToGraph(TH1F* h, double start_val = 2, double stop_val = 5.4) {
  TGraph * gr = new TGraph();
  int ib = h->GetXaxis()->FindBin(start_val);
  int ie = h->GetXaxis()->FindBin(stop_val);
  
  for (int i = ib; i <= ie;i++) {
    gr->SetPoint(gr->GetN(),h->GetXaxis()->GetBinCenter(i),h->GetBinContent(i));
  }
  return gr;
}



void effective_area() {
  f.ls();
  TCanvas* c = new TCanvas();
  c->SetGrid();
  c->SetLogy();
  c->SetTicks();
  TGraph* g_area_sgso_in = HistToGraph((TH1F*)f.Get("hArea_sgso_i_g"),1,5.5);
  TGraph* g_area_sgso_out =HistToGraph((TH1F*)f.Get("hArea_sgso_o_g"),1,5.5);
  
  TH1F* hSum_g = (TH1F*)f.Get("hArea_sgso_i_g");
  hSum_g->Add((TH1F*)f.Get("hArea_sgso_o_g"));
  TGraph* g_area_sgso = HistToGraph(hSum_g,1,5.5);
  g_area_sgso->SetLineColor(kRed);
  g_area_sgso_in->SetLineColor(kRed);
  g_area_sgso_out->SetLineColor(kRed);
  g_area_sgso_in->SetLineWidth(3);
  g_area_sgso_out->SetLineWidth(3);
  g_area_sgso->SetLineWidth(3);

  g_area_sgso->Draw("AL");
  g_area_sgso_out->Draw("samel");
  g_area_sgso_in->Draw("SAMEL");
  
  
  TGraph* g_area_sgso_in_p = HistToGraph((TH1F*)f.Get("hArea_sgso_i_p"),1,5.5);
  TGraph* g_area_sgso_out_p =HistToGraph((TH1F*)f.Get("hArea_sgso_o_p"),1,5.5);
  TH1F* hSum_p = (TH1F*)f.Get("hArea_sgso_i_p");
  hSum_p->Add((TH1F*)f.Get("hArea_sgso_o_p"));
  TGraph* g_area_sgso_p = HistToGraph(hSum_p,1,5.5);
  
  
  g_area_sgso_in_p->SetLineWidth(3);
  g_area_sgso_out_p->SetLineWidth(3);
  g_area_sgso_p->SetLineWidth(3);
  
  g_area_sgso_in_p->SetLineColor(kBlue);
  g_area_sgso_out_p->SetLineColor(kBlue);
  g_area_sgso_p->SetLineColor(kBlue);
  g_area_sgso_p->Draw("samel");
  g_area_sgso_in_p->Draw("SAMEL");
  g_area_sgso_out_p->Draw("samel");


  
}

