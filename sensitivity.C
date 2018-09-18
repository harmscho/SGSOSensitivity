//
// Sensitivity Calculations for SGSO
//
//This program estimates the sensitivity of a strawmans design of the SGSO.
//It loads performance figures and calculates integral (norm_integral_flux) or differetial sensitvity (diff_sens) for 5-sigma detections.
// A given energy spectrum  will be transformed into event-rate in "reconstructed" gamma-ray energy.
//
// source size and scaling factors  for background, point-spread-function, and area are optional arguments for the calculations
//
// harmscho@mpi-hd.mpg.de
//
TH2F* hEtrueErec_sgso_o_g= NULL;
TH2F* hEtrueErec_sgso_i_g= NULL;
TH2F* hEtrueErec_hawc_g= NULL;
TH2F* hEtrueErec_sgso_o_p= NULL;
TH2F* hEtrueErec_sgso_i_p= NULL;
TH2F* hEtrueErec_hawc_p= NULL;
TH1F* hRes_sgso_i_g= NULL;
TH1F* hRes_sgso_o_g= NULL;
TH1F* hRes_hawc_g= NULL;

TH1F* hArea_hawc_g= NULL;
TH1F* hArea_hawc_p= NULL;
TH1F* hArea_sgso_i_g= NULL;
TH1F* hArea_sgso_o_g= NULL;
TH1F* hArea_sgso_i_p= NULL;
TH1F* hArea_sgso_o_p= NULL;

const double TeV2Ergs = 1.60218;
const double MIN = 60;
const double HOUR = MIN *60;
const double DAY = 24 *HOUR;
const double YEAR = 365.25 * DAY;

//Open up the file with performance figures.
TFile* gPerFile = NULL;
void open_file() {
  gPerFile = new TFile("performance.root");
  hEtrueErec_sgso_o_g = (TH2F*)gPerFile->Get("hEtrueErec_sgso_o_g");
  hEtrueErec_sgso_i_g = (TH2F*)gPerFile->Get("hEtrueErec_sgso_i_g");
  hEtrueErec_hawc_g = (TH2F*)gPerFile->Get("hEtrueErec_hawc_g");
  hEtrueErec_sgso_o_p = (TH2F*)gPerFile->Get("hEtrueErec_sgso_o_p");
  hEtrueErec_sgso_i_p = (TH2F*)gPerFile->Get("hEtrueErec_sgso_i_p");
  hEtrueErec_hawc_p = (TH2F*)gPerFile->Get("hEtrueErec_hawc_p");
  hRes_sgso_i_g = (TH1F*)gPerFile->Get("hRes_sgso_i_g");
  hRes_sgso_o_g = (TH1F*)gPerFile->Get("hRes_sgso_o_g");
  hRes_hawc_g = (TH1F*)gPerFile->Get("hRes_hawc_g");
  hArea_hawc_g = (TH1F*)gPerFile->Get("hArea_hawc_g");
  hArea_hawc_p = (TH1F*)gPerFile->Get("hArea_hawc_p");
  hArea_sgso_i_g = (TH1F*)gPerFile->Get("hArea_sgso_i_g");
  hArea_sgso_o_g = (TH1F*)gPerFile->Get("hArea_sgso_o_g");
  hArea_sgso_i_p = (TH1F*)gPerFile->Get("hArea_sgso_i_p");
  hArea_sgso_o_p = (TH1F*)gPerFile->Get("hArea_sgso_o_p");
}
//
//Calculate the expected rate for a given spectrum
//Input spectrum should be in f(E) [TeV^-1 m^-2 s^-1], with E in [TeV]
//
TH1F* get_array_rate(TH2F* hEtrueEreco, TH1F* hArea,TF1* fSpec,TString sTitle = "") {
  TH1F* hRate =(TH1F*)hArea->Clone(sTitle);
  hRate->Reset();
  for (int iy = 1; iy <= hEtrueEreco->GetNbinsY(); iy++) {
    TAxis* trueAxis = hEtrueEreco->GetYaxis();//In GeV
    double Etrue = pow(10,trueAxis->GetBinCenter(iy));// GeV
    TH1D* hReco = hEtrueEreco->ProjectionX("hReco",iy,iy);
    //skip below threshold
    if (hReco->GetEntries() == 0) {
      delete hReco;
      continue;
    }
    hReco->Scale(1./hReco->GetEntries());//normalize
    double dE = pow(10,trueAxis->GetBinUpEdge(iy)) - pow(10,trueAxis->GetBinLowEdge(iy));
    //    cout << Etrue * 1e-3 << "\t" <<fSpec->Eval(Etrue * 1e-3) << endl;
    double F = fSpec->Eval(Etrue*1e-3) * dE * hArea->GetBinContent(hArea->FindBin(trueAxis->GetBinCenter(iy)));
    F *= 1e-3; // histrograms are stored in GeV.
    TAxis* recoAxis = hRate->GetXaxis();
    for (int id = 1; id <= hRate->GetNbinsX(); id++) {
      //rate histogram is in GeV
      double dEreco = pow(10,recoAxis->GetBinUpEdge(id)) - pow(10,recoAxis->GetBinLowEdge(id));
      hRate->Fill(recoAxis->GetBinCenter(id),F * hReco->GetBinContent(id)/dEreco);
    }
    delete hReco;
  }
  return hRate;
}

//
//Estimate the Background in a logE bin.
//
double get_background(double logE /*GeV*/,double dlE /*GeV*/,TH1F* hRate,TH1F* hRes ,double T = 1, double rescaleBKG  = 1, double rescalePSF = 1, double src_radius_68 = 0) {
  int iStart = hRate->FindBin(logE-dlE/2);
  int iStop = hRate->FindBin(logE+dlE/2);
  double dEtrue = pow(10,logE + dlE/2) - pow(10,logE - dlE/2);
  double dEbin = pow(10,hRate->GetXaxis()->GetBinUpEdge(iStop)) - pow(10,hRate->GetXaxis()->GetBinLowEdge(iStart));
  double corrFactor = dEtrue/dEbin; //few percent  correction factor to account for mismatch in binning
  double sum = 0;
  for (int i = iStart; i <=iStop; i++) {
    double dE =  pow(10,hRate->GetXaxis()->GetBinUpEdge(i)) - pow(10,hRate->GetXaxis()->GetBinLowEdge(i));
    double rate = hRate->GetBinContent(i);
    double psf = hRes->GetBinContent(i) * rescalePSF ;
    
    double radius = sqrt(psf*psf + src_radius_68 * src_radius_68);
    double solid_angle = 2 * TMath::Pi() * (1. - cos(radius * TMath::DegToRad()));
    sum += rate * dE * T * solid_angle * rescaleBKG;
  }
  return sum * corrFactor;
}

//
//Estimate a signal in a logE bin
//
double get_signal(double logE /*GeV*/,double dlE /*GeV*/, TH1F* hRate, double T = 1, bool correctForPSF = true) {
  int iE = hRate->FindBin(logE);
  int iStart = hRate->FindBin(logE-(dlE/2));
  int iStop = hRate->FindBin(logE + (dlE/2));
  double dEtrue = pow(10,logE + (dlE/2)) - pow(10,logE - (dlE/2));
  double dEbin = pow(10,hRate->GetXaxis()->GetBinUpEdge(iStop)) - pow(10,hRate->GetXaxis()->GetBinLowEdge(iStart));
  double corrFactor = dEtrue/dEbin; //few percent correction factor to account for mismatch in binning
  
  if (correctForPSF) corrFactor *= 0.68;
  
  double sum = 0;
  for (int i = iStart; i <=iStop; i++) {
    double dE =  pow(10,hRate->GetXaxis()->GetBinUpEdge(i)) - pow(10,hRate->GetXaxis()->GetBinLowEdge(i));
    double rate = hRate->GetBinContent(i);
    
    sum += rate  * dE * T; // 0.68 for using the PSF
  }
  return sum*corrFactor;
}

//
// Helper function to combine two diff sens curves
//
TH1F* combine_diff_sens(TH1F* h1, TH1F* h2) {
  TH1F* hCom = (TH1F*)h1->Clone();
  for (int i = 1; i <= h1->GetNbinsX(); i++) {
    double val1 = h1->GetBinContent(i);
    double val2 = h2->GetBinContent(i);
    if (val1 != 0 && val2 != 0) {
      if (val1 < val2) {
        hCom->SetBinContent(i,val1);
      } else hCom->SetBinContent(i,val2);
    } else if (val1 == 0 && val2 != 0) {
      hCom->SetBinContent(i,val2);
    } else if (val2 == 0 && val1 != 0) {
      hCom->SetBinContent(i,val1);
    }
  }
  return hCom;
}
//
//cta-south diff sense https://www.cta-observatory.org/science/cta-performance/
//
TGraph* diff_sens_cta() {
  TFile* fCTA = new TFile("CTA-Performance-South-50h_20150511.root");
  TH1F* h = (TH1F*)fCTA->Get("DiffSens");
  TH1F* hGev = new TH1F("hDiff",";log_{10} (E_{E_{#gamma '} / GeV);E_{#gamma '}^{2}#times F_{ref} (ergs cm^{-2} s^{-1})",5*5,1,6);
  hGev->Reset();
  TGraph* g = new TGraph();
  for (int i = 1; i <= hGev->GetNbinsX(); i++) {
    double gev = hGev->GetBinCenter(i);
    double bin = h->FindBin(gev-3);
    if (h->GetBinContent(bin) != 0) {
      hGev->SetBinContent(i,h->GetBinContent(bin));
      g->SetPoint(g->GetN(),gev,h->GetBinContent(bin));
    }
  }
  g->SetLineColor(kBlue+2);
  g->SetLineWidth(2);
  return g;
}

//
//get the hawc sensitivy http://iopscience.iop.org/article/10.3847/1538-4357/aa7555/meta
//
TGraph* diff_sens_hawc() {
  ifstream ifs("hawc_crab.txt");
  double Phi,E;
  TGraph* g = new TGraph();
  while (ifs >> E >> Phi) {
    g->SetPoint(g->GetN(), log10(E)+3, TeV2Ergs * E * E* Phi);
  }
  g->SetLineWidth(3);
  g->SetLineColor(kRed);
  g->SetLineColor(kMagenta+1);
  return g;
}

//
// simple helper function to get a graph from a histogram
//
TGraph* hist_to_graph(TH1F* h) {
  TGraph* g = new TGraph();
  for (int i = 1; i <= h->GetNbinsX();i++) {
    if (h->GetBinContent(i) != 0)
      g->SetPoint(g->GetN(),h->GetBinCenter(i),h->GetBinContent(i));
  }
  g->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  g->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
  g->SetLineWidth(3);
  return g;
}

//
//Calculate the differential sensitivity for different arrays.
//
TH1F* diff_sens(TF1* fRefSpec, TString sDet = "hawc", double period = YEAR, double frac_per_day = 0.25, double scalePSF = 1, double scaleBKG = 1, double scaleArea = 1, double sourceSize = 0, double nsigma = 5){
  if (gPerFile == NULL) open_file();
  
  TH1F* hDiff = new TH1F("",";log_{10} (E_{#gamma '} / GeV);E_{#gamma '}^{2} #times #Phi_{ref} (ergs cm^{-2} s^{-1})",5*5-3,1,5.4);
  double T = period * frac_per_day; //six hours per day
  TH2F* hDisp_g = NULL;
  TH2F* hDisp_p = NULL;
  TH1F* hArea_g = NULL;
  TH1F* hArea_p = NULL;
  TH1F* hRes = NULL;
  if (sDet == "hawc") {
    hDisp_g = hEtrueErec_hawc_g;
    hDisp_p = hEtrueErec_hawc_p;
    hArea_g = (TH1F*)hArea_hawc_g->Clone();
    hArea_p = (TH1F*)hArea_hawc_p->Clone();
    hRes = hRes_hawc_g;
  } else if (sDet == "sgso_i"){
    hDisp_g = hEtrueErec_sgso_i_g;
    hDisp_p = hEtrueErec_sgso_i_p;
    hArea_g = (TH1F*)hArea_sgso_i_g->Clone();
    hArea_p = (TH1F*)hArea_sgso_i_p->Clone();
    hRes = hRes_sgso_i_g;
  } else if (sDet == "sgso_o"){
    hDisp_g = hEtrueErec_sgso_o_g;
    hDisp_p = hEtrueErec_sgso_o_p;
    hArea_g = (TH1F*)hArea_sgso_o_g->Clone();
    hArea_p = (TH1F*)hArea_sgso_o_p->Clone();
    hRes = hRes_sgso_o_g;
  } else if (sDet == "sgso") {
    TH1F* hDiffSensInner = diff_sens(fRefSpec,"sgso_i",period,frac_per_day,scalePSF,scaleBKG,scaleArea,sourceSize);
    TH1F* hDiffSensOuter = diff_sens(fRefSpec,"sgso_o",period,frac_per_day,scalePSF,scaleBKG,scaleArea,sourceSize);
    return combine_diff_sens(hDiffSensInner,hDiffSensOuter);
  } else {
    cout << "Request non-existin detector: "<< sDet << endl;
    return NULL;
  }
  hArea_g->Scale(scaleArea);
  hArea_p->Scale(scaleArea);
  
  ///CR rate: 0.096 (E/TeV)^{-2.7} TeV^{-1} s^{-1} m^{-2} sr^-1
  TF1* fCR = new TF1("fCR","[0]*pow(x/[2],-[1])");
  fCR->SetParameter(0,0.096);
  fCR->SetParameter(1,2.7);
  fCR->SetParameter(2,1);//working in TeV
  TH1F* hRateCR = get_array_rate(hDisp_p,hArea_p,fCR,"hRateCR");
  TH1F* hRateGam = get_array_rate(hDisp_g,hArea_g,fRefSpec,"hRefRate");
  
  for (int i = 1; i <= hDiff->GetNbinsX();i++) {
    double logE = hDiff->GetBinCenter(i);
    double dE = hDiff->GetXaxis()->GetBinWidth(2);
    double bkg = get_background(logE,dE,hRateCR,hRes,T,scaleBKG,scalePSF,sourceSize);
    bool psfCor = !(sourceSize > 0);
    double sig = get_signal(logE,dE,hRateGam,T,psfCor);
    double Eergs = pow(10,logE-3)*TeV2Ergs;
    //rescale ref flux
    double nsig = sqrt(bkg) * nsigma;//5 sigma above background
    if (nsig < 10) nsig = 10;//minimal 10 photons per bin
    double ratio = nsig/sig;
    double F = fRefSpec->Eval(pow(10,logE-3)) * ratio;
    F /= TeV2Ergs;//in ergs
    F /= 1e4;//to cm2
    
    double ds = Eergs * Eergs * F;
    if (bkg != 0) {
      //      cout  <<  std::scientific << setprecision(3) <<
      //      std::setw(10)  << logE <<
      //      std::setw(10) << " bkg: " << bkg <<
      //      std::setw(10) << "nsig " <<nsig <<
      //      std::setw(10) << " sig " << sig <<
      //      std::setw(10) << " F " << F <<  endl;
      hDiff->SetBinContent(i,ds);
    }
  }
  return hDiff;
}


//likelihood only works for SGSO. Combines innner and outer array
TH1F* diff_sens_lh(TF1* fRefSpec, TString sDet = "sgso", double period = YEAR, double frac_per_day = 0.25, double scalePSF = 1, double scaleBKG = 1, double scaleArea = 1, double sourceSize = 0, double nsigma = 5){
  if (gPerFile == NULL) open_file();
  
  TH1F* hDiff = new TH1F("",";log_{10} (E_{#gamma '} / GeV);E_{#gamma '}^{2} #times #Phi_{ref} (ergs cm^{-2} s^{-1})",5*5-3,1,5.4);
  double T = period * frac_per_day; //six hours per day
  TH2F* hDisp_g = NULL;
  TH2F* hDisp_p = NULL;
  TH1F* hArea_g = NULL;
  TH1F* hArea_p = NULL;
  TH1F* hRes = NULL;

  if (sDet != "sgso"){
    cout << "Only works for comnbination of inner and outer (sgso), not for detector: "<< sDet << endl;
    return NULL;
  }
  
  ///CR rate: 0.096 (E/TeV)^{-2.7} TeV^{-1} s^{-1} m^{-2} sr^-1
  TF1* fCR = new TF1("fCR","[0]*pow(x/[2],-[1])");
  fCR->SetParameter(0,0.096);
  fCR->SetParameter(1,2.7);
  fCR->SetParameter(2,1);//working in TeV
  TH1F* hRateCR = NULL;
  TH1F* hRateGam = NULL;
  
  for (int i = 1; i <= hDiff->GetNbinsX();i++) {
    double logE = hDiff->GetBinCenter(i);
    double dE = hDiff->GetXaxis()->GetBinWidth(2);
    bool psfCor = !(sourceSize > 0);
    double Eergs = pow(10,logE-3)*TeV2Ergs;
    //rescale ref flux
    
    double norm = 1000;
    double sumSig = 0;
    //combined sensitivity
    
    //signal outer
    hDisp_g = hEtrueErec_sgso_o_g;
    hDisp_p = hEtrueErec_sgso_o_p;
    hArea_g = (TH1F*)hArea_sgso_o_g->Clone();
    hArea_p = (TH1F*)hArea_sgso_o_p->Clone();
    hRes = hRes_sgso_o_g;
    hRateCR = get_array_rate(hDisp_p,hArea_p,fCR,"hRateCR");
    hRateGam = get_array_rate(hDisp_g,hArea_g,fRefSpec,"hRefRate");
    double sig_out = get_signal(logE,dE,hRateGam,T,psfCor);
    double bkg_out = get_background(logE,dE,hRateCR,hRes,T,scaleBKG,scalePSF,sourceSize);
      
    //signal inner
    hDisp_g = hEtrueErec_sgso_i_g;
    hDisp_p = hEtrueErec_sgso_i_p;
    hArea_g = (TH1F*)hArea_sgso_i_g->Clone();
    hArea_p = (TH1F*)hArea_sgso_i_p->Clone();
    hRes = hRes_sgso_i_g;
    hRateCR = get_array_rate(hDisp_p,hArea_p,fCR,"hRateCR");
    hRateGam = get_array_rate(hDisp_g,hArea_g,fRefSpec,"hRefRate");
    double sig_in = get_signal(logE,dE,hRateGam,T,psfCor);
    double bkg_in = get_background(logE,dE,hRateCR,hRes,T,scaleBKG,scalePSF,sourceSize);
    
    //calculate the normalisation needed of reference spectrum
    //that results in a 5 sigma detection
    norm = 1000;
    double prev_norm = 10 * norm;
    double TS = 1e3;
    double scaleDown = 0.5;
    double prevScaleDown = scaleDown;
    // iteration counter, in case we do not converge...
    int iterations = 0;
    sumSig = 0;
    while ( TMath::Abs(TS - nsigma*nsigma) > 1 && iterations < 500) {
      //Setting stepsizes in flux norm
      if (TS < nsigma*nsigma) {
        norm = prev_norm;
        scaleDown = prevScaleDown/2;
      }
      prev_norm = norm;
      prevScaleDown = scaleDown;
      norm *= (1.-scaleDown);
      // likelihood calculation
      double Lbkg = 0;
      double Lsig = 0;
      
      double Pbkg_in = TMath::Poisson(bkg_in,bkg_in);
      double Psig_in = TMath::Poisson(norm*sig_in + bkg_in,bkg_in);
      if (norm*sig_in != 0 && bkg_in != 0) {
        if (Pbkg_in < 1e-20) Pbkg_in = 1e-20;
        if (Psig_in < 1e-20) Psig_in = 1e-20;
        Lbkg += -log(Pbkg_in);
        Lsig += -log(Psig_in);
      }
      
      double Pbkg_out = TMath::Poisson(bkg_out,bkg_out);
      double Psig_out = TMath::Poisson(norm*sig_out + bkg_out,bkg_out);
      if (norm*sig_out != 0 && bkg_out != 0) {
        if (Pbkg_out < 1e-20) Pbkg_out = 1e-20;
        if (Psig_out < 1e-20) Psig_out = 1e-20;
        Lbkg += -log(Pbkg_out);
        Lsig += -log(Psig_out);
      }
      sumSig = norm*sig_out + norm*sig_in;
      TS = 2 * (Lsig - Lbkg);
      iterations++;
    }//end loop
    
    if (sumSig < 10) norm *= 10./sumSig;
    
    double F = fRefSpec->Eval(pow(10,logE-3)) * norm;
    F /= TeV2Ergs;//in ergs
    F /= 1e4;//to cm2
    
    double ds = Eergs * Eergs * F;
    if (bkg_in+bkg_out != 0) {
      //      cout  <<  std::scientific << setprecision(3) <<
      //      std::setw(10)  << logE <<
      //      std::setw(10) << " bkg: " << bkg <<
      //      std::setw(10) << "nsig " <<nsig <<
      //      std::setw(10) << " sig " << sig <<
      //      std::setw(10) << " F " << F <<  endl;
      hDiff->SetBinContent(i,ds);
    }
  }
  return hDiff;
}




//
// function declaration
//
double norm_integral_flux(TF1* fRefSpec, TString sDet = "hawc", double period = YEAR, double frac_per_day = 0.25, double scalePSF = 1, double scaleBKG = 1, double scaleArea = 1, double sourceSize = 0);

//
//small wrapper to combine inner and outer part of sgso.
//
double comb_norm_integral_flux(TF1* fRefSpec, TString sDet = "hawc", double period = YEAR, double frac_per_day = 0.25, double scalePSF = 1, double scaleBKG = 1, double scaleArea = 1, double sourceSize = 0) {
  double normIn = norm_integral_flux(fRefSpec,"sgso_i",period,frac_per_day,scalePSF,scaleBKG,scaleArea,sourceSize);
  double normOut = norm_integral_flux(fRefSpec,"sgso_o",period,frac_per_day,scalePSF,scaleBKG,scaleArea,sourceSize);
  if (normIn < normOut) {
    return normIn;
  } else return normOut;
}

//
// Takes as reference spectrum, and returns the scaling factor that would lead to a 5 -sigma detection.
// Reference spectrum f(E) should be in units of [TeV^-1 m^-2 s^-1], with E in [TeV]
//
double norm_integral_flux(TF1* fRefSpec, TString sDet, double period, double frac_per_day, double scalePSF, double scaleBKG, double scaleArea, double sourceSize) {
  
  //histogram to use binning from
  TH1F* hSens = new TH1F("",";log_{10} (E_{#gamma '} / GeV);E_{#gamma '}^{2} #times #Phi_{ref} (ergs cm^{-2} s^{-1})",5*5,1,6);
  
  //getting performance figures
  TH2F* hDisp_g = NULL;
  TH2F* hDisp_p = NULL;
  TH1F* hArea_g = NULL;
  TH1F* hArea_p = NULL;
  TH1F* hRes = NULL;
  if (sDet == "hawc") {
    hDisp_g = hEtrueErec_hawc_g;
    hDisp_p = hEtrueErec_hawc_p;
    hArea_g = (TH1F*)hArea_hawc_g->Clone();
    hArea_p = (TH1F*)hArea_hawc_p->Clone();
    hRes = hRes_hawc_g;
  } else if (sDet == "sgso_i"){
    hDisp_g = hEtrueErec_sgso_i_g;
    hDisp_p = hEtrueErec_sgso_i_p;
    hArea_g = (TH1F*)hArea_sgso_i_g->Clone();
    hArea_p = (TH1F*)hArea_sgso_i_p->Clone();
    hRes = hRes_sgso_i_g;
  } else if (sDet == "sgso_o"){
    hDisp_g = hEtrueErec_sgso_o_g;
    hDisp_p = hEtrueErec_sgso_o_p;
    hArea_g = (TH1F*)hArea_sgso_o_g->Clone();
    hArea_p = (TH1F*)hArea_sgso_o_p->Clone();
    hRes = hRes_sgso_o_g;
  } else if (sDet == "sgso") {
    return comb_norm_integral_flux(fRefSpec,"",period,frac_per_day,scalePSF,scaleBKG,scaleArea,sourceSize);
  } else {
    cout << "Request non-existin detector: "<< sDet << endl;
    return 0;
  }
  hArea_g->Scale(scaleArea);
  hArea_p->Scale(scaleArea);
  double T = period * frac_per_day; // total observation time
  
  ///CR rate: 0.096 (E/TeV)^{-2.7} TeV^{-1} s^{-1} m^{-2} sr^-1
  TF1* fCR = new TF1("fCR","[0]*pow(x/[2],-[1])");
  fCR->SetParameter(0,0.096);
  fCR->SetParameter(1,2.7);
  fCR->SetParameter(2,1);
  //Getting Array background and gamma-ray rates
  TH1F* hRateCR = get_array_rate(hDisp_p,hArea_p,fCR,"hRateCR");
  TH1F* hRateGam = get_array_rate(hDisp_g,hArea_g,fRefSpec,"hRefRate");
  
  //calculate the normalisation needed of reference spectrum
  //that results in a 5 sigma detection
  double norm = 1000;
  double prev_norm = 10 * norm;
  double TS = 1e3;
  double scaleDown = 0.5;
  double prevScaleDown = scaleDown;
  // iteration counter, in case we do not converge...
  int iterations = 0;
  while ( TMath::Abs(TS - 25) > 1 && iterations < 500) {
    //Setting stepsizes in flux norm
    if (TS < 25) {
      norm = prev_norm;
      scaleDown = prevScaleDown/2;
    }
    prev_norm = norm;
    prevScaleDown = scaleDown;
    norm *= (1.-scaleDown);
    
    // likelihood calculation
    double Lbkg = 0;
    double Lsig = 0;
    double sumSig = 0;
    double sumBkg= 0;
    for (int i = 1; i <= hSens->GetNbinsX(); i++) {
      double logE = hSens->GetBinCenter(i);
      double dE = hSens->GetXaxis()->GetBinWidth(2);
      double bkg = get_background(logE,dE,hRateCR,hRes,T,scaleBKG,scalePSF,sourceSize);
      sumBkg += bkg;
      bool psfCor = !(sourceSize > 0);
      double sig = norm * get_signal(logE,dE,hRateGam,T);
      sumSig += sig;
      double Pbkg = TMath::Poisson(bkg,bkg);
      double Psig = TMath::Poisson(sig + bkg,bkg);
      if (sig != 0 && bkg != 0) {
        if (Pbkg < 1e-20) Pbkg = 1e-20;
        if (Psig < 1e-20) Psig = 1e-20;
        Lbkg += -log(Pbkg);
        Lsig += -log(Psig);
      }
    }
    TS = 2 * (Lsig - Lbkg);
    iterations++;
  }//end loop
  
  if (iterations == 500) {
    cout << "WARNING: couldn't find 5-sigma limit" << endl;
  }
  delete hSens;
  return norm ;
}

//
//Make Differential Sensitivity Figure
//
//
//Make Differential Sensitivity Figure
//
void diff_sens_figure() {
  
  //Taking a simple power-law crab spectrum as a reference source
  ///crab = 3.2e-7 * (E/TeV)^2.49 m^2 s^-1 TeV^-1 http://iopscience.iop.org/article/10.1086/306005/meta
  TF1* fCrab = new TF1("fCrab","[0]*pow(x/[2],-[1])");
  fCrab->SetParameter(0,3.2e-7); //m^2 s^-1 TeV^-1
  fCrab->SetParameter(1,2.5);// m^2 s^-1 TeV^-1
  fCrab->SetParameter(2,1); // in TeV
  
  //scaling factors for background, point spread functions, and instrumented area
  double rPSF = 1;
  double rBKG = 1;
  double rA = 1;
  //
  TH1F* hSimHAWC = diff_sens(fCrab,"hawc",507*DAY,6.*HOUR/DAY,rPSF,rBKG,rA);
  
  //setting up the canvas
  TCanvas* cDiff = new TCanvas();
  cDiff->SetLogy();
  cDiff->SetGrid();
  cDiff->SetTicks();
  
  TGraph* gCTA = diff_sens_cta();
  gCTA->GetXaxis()->SetLabelSize(0.05);
  gCTA->GetYaxis()->SetLabelSize(0.05);
  gCTA->GetXaxis()->SetTitleSize(0.05);
  gCTA->GetYaxis()->SetTitleSize(0.05);
  gCTA->GetYaxis()->SetTitleOffset(1.20);
  gCTA->GetXaxis()->SetTitle("log_{10} (E_{ r} / GeV)");
  gCTA->GetYaxis()->SetTitle("E^{2}#times Flux Sens. [ergs cm^{-2} s^{-1}]");
  gCTA->SetMinimum(1e-14);
  gCTA->SetMaximum(1e-10);
  gCTA->Draw("AL");
  
  //Asume a factor 0.5 less hadrons, and 0.75 in point spread function for SGSO
  rA = 1;
  rPSF = 1;
  rBKG = 0.5;
  
  //1 year of SGSO
  TH1F* hSGSO = diff_sens_lh(fCrab,"sgso",YEAR,6.*HOUR/DAY,1,rBKG,rA);
  TGraph* gSGSO = hist_to_graph(hSGSO);
  gSGSO->SetLineColor(kRed+2);
  //  gSGSO->Draw("SAMEL");
  
  //  TH1F* hSGSO_lh = diff_sens_lh(fCrab,"sgso_c",YEAR,6.*HOUR/DAY,1,rBKG,rA);
  //  TGraph* gSGSO_lh = hist_to_graph(hSGSO_lh);
  
  
  //convervative
  TH1F* hSGSO_con= diff_sens_lh(fCrab,"sgso",YEAR,6.*HOUR/DAY,1,1,rA);
  TGraph* gSGSO_con = hist_to_graph(hSGSO_con);
  gSGSO_con->SetLineColor(kRed+3);
  //  gSGSO_con->Draw("SAMEL");
  
  //optimistic
  TH1F* hSGSO_opt= diff_sens_lh(fCrab,"sgso",YEAR,6.*HOUR/DAY,0.8,0.5,rA);
  TGraph* gSGSO_opt = hist_to_graph(hSGSO_opt);
  gSGSO_opt->SetLineColor(kRed);
  
  
  TGraphAsymmErrors* gr_1year = new TGraphAsymmErrors();
  for (int i = 0; i < gSGSO_opt->GetN(); i++) {
    gr_1year->SetPoint(i,gSGSO->GetX()[i],gSGSO->GetY()[i]);
    gr_1year->SetPointEYhigh(i,gSGSO_con->GetY()[i]-gSGSO->GetY()[i]);
    gr_1year->SetPointEYlow(i,gSGSO->GetY()[i] - gSGSO_opt->GetY()[i]);
  }
  gr_1year->SetFillColorAlpha(2,0.4);
  gr_1year->Draw("same 4");
  gSGSO_opt->Draw("SAMEL");
  gSGSO_con->SetLineColor(2);
  gSGSO_con->Draw("SAMEL");
  //  gSGSO_lh->Draw("samel");
  
  cout << "\n\n...5 year scenario... " << endl;
  cout << "...Middle... " << endl;
  //5 years of SGSO
  TH1F* hSGSO_5yr = diff_sens_lh(fCrab,"sgso",5 * YEAR,6.*HOUR/DAY,rPSF,rBKG,rA);
  TGraph* gSGSO_5yr = hist_to_graph(hSGSO_5yr);
  gSGSO_5yr->SetLineColor(kRed+2);
  gSGSO_5yr->SetLineStyle(2);
  //  gSGSO_5yr->Draw("same L");
  
  //convervative
  cout << "...convervative... " << endl;
  TH1F* hSGSO_5yr_cons=  diff_sens_lh(fCrab,"sgso",5 * YEAR,6.*HOUR/DAY,1,1,rA);
  TGraph* gSGSO_5yr_cons = hist_to_graph(hSGSO_5yr_cons);
  gSGSO_5yr_cons->SetLineColor(kRed+3);
  //  gSGSO_5yr_cons->SetLineStyle(2);
  //  gSGSO_5yr_cons->Draw("same L");
  cout << "...best ... " << endl;
  TH1F* hSGSO_5yr_opt= diff_sens_lh(fCrab,"sgso",5*YEAR,6.*HOUR/DAY,0.8,0.5,rA);
  TGraph* gSGSO_5yr_opt = hist_to_graph(hSGSO_5yr_opt);
  //  gSGSO_5yr_opt->SetLineStyle(2);
  //  gSGSO_5yr_opt->SetLineColor(kRed);
  //  gSGSO_5yr_opt->Draw("SAMEL");
  
  
  TGraphAsymmErrors* gr_5year = new TGraphAsymmErrors();
  for (int i = 0; i < gSGSO_5yr->GetN(); i++) {
    gr_5year->SetPoint(i,gSGSO_5yr->GetX()[i],gSGSO_5yr->GetY()[i]);
    gr_5year->SetPointEYhigh(i,gSGSO_5yr_cons->GetY()[i]-gSGSO_5yr->GetY()[i]);
    gr_5year->SetPointEYlow(i,gSGSO_5yr->GetY()[i] - gSGSO_5yr_opt->GetY()[i]);
  }
  gr_5year->SetFillColorAlpha(kRed+3,0.4);
  gr_5year->Draw("same 4");
  gSGSO_5yr_opt->SetLineColor(kRed+3);
  gSGSO_5yr_opt->Draw("SAMEL");
  gSGSO_5yr_cons->SetLineColor(kRed+3);
  gSGSO_5yr_cons->Draw("SAMEL");
  
  
  // the published sensitivity from HAWC-crab paper
  TGraph* gHAWCCrab = diff_sens_hawc();
  gHAWCCrab->Draw("same L");
}

//
//integral power law sensitivity figure
//
void pl_detect_figure() {
  TF1* fPowerLaw = new TF1("fPowerLaw","[0]*pow(x/[2],-[1])",200,1e6);
  double crabNorm = 3.2e-7;
  fPowerLaw->SetParameter(0,3.2e-7);//m^2 s^-1 TeV^-1
  fPowerLaw->SetParameter(1,2.5);
  fPowerLaw->SetParameter(2,1);//TeV
  
  double rA = 1;
  double rPSF = 0.75;
  double rBKG = 0.5;
  
  TCanvas* can1 = new TCanvas("can1","can1",500,450);
  can1->SetGrid();
  can1->SetTicks();
  can1->SetLogy();
  can1->SetRightMargin(0.1);
  
  TGraph* g;
  double color = 51;
  double dc = 50./40;
  double specIndex = 0;
  for (int i =0; i < 40; i++) {
    fPowerLaw->SetParameter(0,crabNorm);
    fPowerLaw->SetParameter(1,specIndex);
    double normI =norm_integral_flux(fPowerLaw,"sgso",1*YEAR,6.*HOUR/DAY,rPSF,rBKG,rA);
    //    cout << i << "\t" << specIndex  << "\t" <<  normI * crabNorm << "\t" <<  normI << "\t" << color <<endl;
    fPowerLaw->SetParameter(0,crabNorm * normI);
    g = new TGraph();
    double dx = 0.05;
    for (int j = 0; j < 80; j++) {
      double E = pow(10,-1+j*dx);//E in TeV
                                 //Xaxis in GeV / Yaxis in ergs cm-2 s-1.
      g->SetPoint(j,log10(E)+3,TeV2Ergs * pow(E,2)/1e4 * fPowerLaw->Eval(E));
    }
    g->SetLineColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(4);
    if (i==0) {
      g->SetMinimum(1e-15);
      g->SetMaximum(1e-10);
      g->GetXaxis()->SetLabelSize(0.05);
      g->GetYaxis()->SetLabelSize(0.05);
      g->GetXaxis()->SetTitleSize(0.05);
      g->GetYaxis()->SetTitleSize(0.05);
      g->GetYaxis()->SetTitleOffset(1.35);
      g->GetYaxis()->SetTitle("E^{2} dN/dE [ergs cm^{-2} s^{-1}] ");
      g->GetXaxis()->SetTitle("log_{10}( E / GeV)");
      g->DrawClone("AL");
    } else {
      g->DrawClone("SAMEL");
    }
    specIndex += 0.1;
    color += dc;
  }
}

//
// Estimates the surface brightness needed for a 3 sigma
// differential detection.
// for a dirty estimation of fermi-bubbles
void large_source_figure(int ndeg = 10) {
  cout << "\n Large source figure (approximation for FermiBubble)..." << endl;
  TCanvas* cLarge = new TCanvas();
  cLarge->SetLogy();
  cLarge->SetGrid();
  double rA = 1;
  double rPSF = 0.75;
  double rBKG = 0.5;
  TF1* fCrab = new TF1("fCrab","[0]*pow(x/[2],-[1])");
  fCrab->SetParameter(0,3.2e-7); //m^2 s^-1 TeV^-1
  fCrab->SetParameter(1,2.5);// m^2 s^-1 TeV^-1
  fCrab->SetParameter(2,1); // in TeV
  
  double color = 51;
  double dc = 49./ndeg;
  for (int i = 1; i <= ndeg; i++ ) {
    double sourceRadius = 2*i;
    
    //Only using the inner part
    TH1F* hSGSO10deg = diff_sens(fCrab,"sgso_i",5*YEAR,6.*HOUR/DAY,rPSF,rBKG,rA,sourceRadius,3.);
    TGraph* gSGSO10deg = hist_to_graph(hSGSO10deg);
    
    double solid_angle = 2 * TMath::Pi() * (1 - cos(sourceRadius * TMath::DegToRad()));
    gSGSO10deg->GetYaxis()->SetTitle("[GeV cm^{-2} s^{-1} sr^{-1}]");
    for (int j = 0; j < gSGSO10deg->GetN(); j++ ) {
      gSGSO10deg->GetY()[j] *= 1.e3/(solid_angle * TeV2Ergs);
      //fermi bubbles have roughly 20 degrees radius
      if (sourceRadius == 20)
        cout << pow(10,gSGSO10deg->GetX()[j])  << ", " << gSGSO10deg->GetY()[j] << endl;
    }
    gSGSO10deg->SetLineColor(color);
    gSGSO10deg->SetMinimum(1e-9);
    gSGSO10deg->SetMaximum(2e-6);
    if (i == 1) {
      gSGSO10deg->DrawClone("AL");
    } else {
      gSGSO10deg->DrawClone("SAMEL");
    }
    color += dc;
  }
}

// TeV emission from CR interaction with nearby clouds.
// Print the rescaling needed of the typical flux that would lead to
// a 5 sigma detection.
//clouds spectra from https://arxiv.org/pdf/1112.5541.pdf
void clouds() {
  cout << "\nCloud estimation..." << endl;
  TF1* fPowerLaw = new TF1("fPowerLaw","[0]*pow(x/[2],-[1])",200,1e6);
  
  //norm at 10 GeV
  double norm = 1e-3; // in (TeV m2 s-1)
  double refEnergy = 1e-2; // in GeVs
  
  //index
  double index = 2.7;
  //
  fPowerLaw->SetParameter(0,norm);//m^2 s^-1 TeV^-1
  fPowerLaw->SetParameter(1,index);
  fPowerLaw->SetParameter(2,refEnergy);//TeV
  
  //typical source radius of the clouds
  double sourceRadius = 5;
  double rPSF = 1;
  double rBKG = 0.5;
  
  double normI = norm_integral_flux(fPowerLaw,"sgso",5*YEAR,6.*HOUR/DAY,rPSF,rBKG,1,sourceRadius);
  cout << "Flux needs to be  " << normI
  << " times higher than typical flux for 5sigma detection" << endl;
  
}


void duration_detectability() {
  //
  //integral power law sensitivity figure
  //
  
  TF1* fPowerLaw = new TF1("fPowerLaw","[0]*pow(x/[2],-[1])",200,1e6);
  double crabNorm = 3.2e-7;
  fPowerLaw->SetParameter(0,3.2e-7);//m^2 s^-1 TeV^-1
  fPowerLaw->SetParameter(1,2.5);
  fPowerLaw->SetParameter(2,1);//TeV
  
  double rA = 1;
  double rPSF = 0.75;
  double rBKG = 0.5;
  
  TCanvas* can1 = new TCanvas("can1","can1",500,450);
  can1->SetGrid();
  can1->SetTicks();
  can1->SetLogy();
  can1->SetLogx();
  can1->SetRightMargin(0.1);
  
  const int nDur = 20;
  double dDur = 6./nDur;
  
  TGraph* g;
  double color = 51;
  double dc = 48./4;
  double specIndex = 2.0;
  for (int iSpec = 0; iSpec < 5; iSpec++) {
    double duration = 0;
    fPowerLaw->SetParameter(0,crabNorm);
    fPowerLaw->SetParameter(1,specIndex);
    g = new TGraph();
    for (int i =0; i < nDur; i++) {
      
      double normI =norm_integral_flux(fPowerLaw,"sgso",pow(10,duration),1,rPSF,rBKG,rA);
      
      g->SetPoint(i,pow(10,duration),normI*fPowerLaw->Eval(0.1));
      duration += dDur;
    }
    
    g->SetLineColorAlpha(color,0.5);
    g->SetLineWidth(4);
    if (iSpec==0) {
      g->SetMinimum(1e-6);
      g->SetMaximum(1e-2);
      g->GetXaxis()->SetLabelSize(0.05);
      g->GetYaxis()->SetLabelSize(0.05);
      g->GetXaxis()->SetTitleSize(0.05);
      g->GetXaxis()->SetTitleOffset(1.35);
      g->GetYaxis()->SetTitleSize(0.05);
      g->GetYaxis()->SetTitleOffset(1.35);
      g->GetYaxis()->SetTitle("#phi(100 GeV) [m^{2} s^{-1} TeV^{-1}]");
      g->GetXaxis()->SetTitle("duration [s]");
      g->DrawClone("AL");
    } else {
      g->DrawClone("SAMEL");
    }
    specIndex += 0.5;
    color += dc;
  }
  
  
  //PKS 2155-304
  cout << "PKS2155 detected in ";
  TF1* fSpectrumCut = new TF1("PKS2155Cut","2.06e-6*TMath::Min(pow(x,-2.71),pow(0.43,0.82)*pow(x,-3.53))*(0.5+TMath::Sign(0.5,5.0-x))",0.2,1.0e3);
  for (int i = 10; i < 500; i++) {
    double norm = norm_integral_flux(fSpectrumCut,"sgso",i,1,rPSF,rBKG,rA);
    if (norm < 1) {
      cout << i << " seconds,  flux at 100 GeV is " << fSpectrumCut->Eval(0.1) << " [m2 TeV s]^1"  << endl;
      break;
    }
  }

  //fCrab
  fPowerLaw->SetParameter(0,crabNorm);//m^2 s^-1 TeV^-1
  fPowerLaw->SetParameter(1,2.5);
  fPowerLaw->SetParameter(2,1);//TeV
  cout << "Crab detected in ";
  for (int i = 540; i < 1000; i++) {
    double norm = norm_integral_flux(fPowerLaw,"sgso",i,1,rPSF,rBKG,rA);
    if (norm < 1) {
      cout << i << " seconds,  flux at 100 GeV is " << fPowerLaw->Eval(0.1) << " [m2 TeV s]^1"  << endl;
      break;
    }
  }

  
}

//
// Function to get figures
//
void sensitivity() {
  open_file();
  diff_sens_figure();
  /*pl_detect_figure();
  large_source_figure();
  clouds();
  duration_detectability();*/
  
}








