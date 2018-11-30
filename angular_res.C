TFile f("performance.root");
TFile fcta("CTA-Performance-South-50h_20150511.root");
TGraph* HistToGraph(TH1F* h, double start_val = 2, double stop_val = 5.4) {
  TGraph * gr = new TGraph();
  int ib = h->GetXaxis()->FindBin(start_val);
  int ie = h->GetXaxis()->FindBin(stop_val);
  
  for (int i = ib; i <= ie;i++) {
    gr->SetPoint(gr->GetN(),h->GetXaxis()->GetBinCenter(i),h->GetBinContent(i));
  }
  return gr;
}



void angular_res() {
  f.ls();
  TCanvas* c = new TCanvas();
  c->SetGrid();
  c->SetLogy();
  c->SetTicks();
  TGraph* g_res_sgso_in = HistToGraph((TH1F*)f.Get("hRes_sgso_i_g"),2.5);
  TGraph* g_res_sgso_out =HistToGraph((TH1F*)f.Get("hRes_sgso_o_g"),2.9);
  TGraph* g_res_hawc =HistToGraph((TH1F*)f.Get("hRes_hawc_g"),2.9);
  TGraph* g_res_cta =HistToGraph((TH1F*)fcta.Get("AngRes"),-1.8,2.2);
  shiftGraph(*g_res_cta,3);
  
  g_res_sgso_in->SetLineColor(kRed);
  g_res_sgso_out->SetLineStyle(7);
  g_res_sgso_out->SetLineColor(kRed);
  g_res_cta->SetLineColor(kBlue+2);
  g_res_sgso_out->SetLineWidth(3);
  g_res_sgso_in->SetLineWidth(3);
  g_res_hawc->SetLineWidth(3);
  g_res_hawc->SetLineWidth(3);
  g_res_cta->SetLineWidth(3);
  g_res_hawc->SetLineColor(kMagenta+1);

  g_res_cta->SetTitle(";log_{10}(E_{R} / GeV);#sigma_{68%} [#circ]");
  g_res_cta->Draw("al");
  g_res_cta->SetMaximum(1.8);
  g_res_cta->SetMinimum(0.02);
  g_res_sgso_in->Draw("samel");
  g_res_sgso_out->Draw("same l");
  g_res_hawc->Draw("same l");
  
  
  
  
}
