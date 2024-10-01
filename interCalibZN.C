#if !defined(__CINT__) || defined(__MAKECINT__)

#include <stdio.h>
#include <stdlib.h>
#include <TROOT.h>
#include <Riostream.h>
#include <TClassTable.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TTree.h>
#include <TLegend.h>
#include <TString.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMinuit.h>
#include <TDirectoryFile.h>

#endif

void setCosm();
void cumulate(int ih, double tc, double t1, double t2, double t3, double t4, double w);
static void fcn(int &np, double *gin, double &f, double *par, int iflag);
int minimization(int ih); 
TCanvas* plotTProfiles(std::vector<TProfile*> profiles, const char* name);

constexpr double kVeryNegative = -1.e12; // assign this when ZDC when energy is negative 
//
constexpr int ndet = 2; 
constexpr int npar = 6; 
//
static double mSum[ndet][npar][npar]; 
static double mAdd[npar][npar];
//
const float xmin = -0.5;
const float xmax = 8000;
//
const float lowbnd = 0.1;
const float uppbnd = 10.;
const float start = 0.8; 
const float step = 0.01;

bool deming = true;
//
double coeff[ndet][npar], coefferr[ndet][npar]; 

TH1F *ZNCcalibcoeff = new TH1F("ZNCcalibcoeff","Calib Coeff ZNC", 4, 0, 4);
ZNCcalibcoeff->GetXaxis()->SetBinLabel(1,"pmq1");
ZNCcalibcoeff->GetXaxis()->SetBinLabel(2,"pmq2");
ZNCcalibcoeff->GetXaxis()->SetBinLabel(3,"pmq3");
ZNCcalibcoeff->GetXaxis()->SetBinLabel(4,"pmq4");

TH1F *ZNAcalibcoeff = new TH1F("ZNAcalibcoeff","Calib Coeff ZNA", 4, 0, 4);
ZNAcalibcoeff->GetXaxis()->SetBinLabel(1,"pmq1");
ZNAcalibcoeff->GetXaxis()->SetBinLabel(2,"pmq2");
ZNAcalibcoeff->GetXaxis()->SetBinLabel(3,"pmq3");
ZNAcalibcoeff->GetXaxis()->SetBinLabel(4,"pmq4");


void interCalibZN (TString inputfilename = "treeZNpms.root")
{
  TFile *fin = TFile::Open(inputfilename.Data(),"READ");
  if(!fin) {
    printf("  +++ FILE %s not found!!!\n\n", inputfilename.Data());
    return;
  }

  TDirectoryFile* df = (TDirectoryFile*)fin->Get("DF_2336995099171200");
  if(!df){
	  printf("Directory missing -> PLEASE CHECK FILE NAME!!! \n");
	  return;
  }

  TTree* tree = (TTree*)df->Get("O2zdcic");
  if(!tree){
	  printf("  + tree O2zdcic is not there!!!\n");
	  return;
  }  

  tree->Print("");
  Int_t nentries = (Int_t) tree->GetEntries();
  printf("      Tree has %d entries\n\n", nentries);


  float pmcZNA;
  float pm1ZNA;
  float pm2ZNA;
  float pm3ZNA;
  float pm4ZNA;
  float pmcZNC;
  float pm1ZNC;
  float pm2ZNC;
  float pm3ZNC;
  float pm4ZNC;
  
  tree->SetBranchAddress("fpmcZNA", &pmcZNA);
  tree->SetBranchAddress("fpm1ZNA", &pm1ZNA);
  tree->SetBranchAddress("fpm2ZNA", &pm2ZNA);
  tree->SetBranchAddress("fpm3ZNA", &pm3ZNA);
  tree->SetBranchAddress("fpm4ZNA", &pm4ZNA);
  tree->SetBranchAddress("fpmcZNC", &pmcZNC);
  tree->SetBranchAddress("fpm1ZNC", &pm1ZNC);
  tree->SetBranchAddress("fpm2ZNC", &pm2ZNC);
  tree->SetBranchAddress("fpm3ZNC", &pm3ZNC);
  tree->SetBranchAddress("fpm4ZNC", &pm4ZNC);

  TH1F *hpmcZNC = new TH1F("hpmcZNC"," PMC ZNC", 40, xmin, xmax);
  TH1F *hpmcZNA = new TH1F("hpmcZNA"," PMC ZNA", 40, xmin, xmax);
  TH1F *hsumZNC = new TH1F("hsumZNC"," uncalibrated SUMq ZNC", 40, xmin, xmax);
  TH1F *hsumZNA = new TH1F("hsumZNA"," uncalibrated SUMq ZNA", 40, xmin, xmax);

  TH1F *hsumZNC_cal = new TH1F("hsumZNC_cal","calibrated SUMq ZNC", 40, xmin, xmax);
  TH1F *hsumZNA_cal = new TH1F("hsumZNA_cal","calibrated SUMq ZNA", 40, xmin, xmax);

  int nbins = 100; 
  TH2F *hpmcZNC_vs_sum = new TH2F("hpmcZNC_vs_sum", "pmcq vs. sum ZNC", nbins, xmin, xmax, nbins, xmin, xmax); 
  TH2F *hpmcZNA_vs_sum = new TH2F("hpmcZNA_vs_sum", "pmcq vs. sum ZNA", nbins, xmin, xmax, nbins, xmin, xmax); 

  TH2F *hpmcZNC_vs_sum_cal = new TH2F("hpmcZNC_vs_sum_cal", "pmcq vs. calibrated sum ZNC", nbins, xmin, xmax, nbins, xmin, xmax); 
  TH2F *hpmcZNA_vs_sum_cal = new TH2F("hpmcZNA_vs_sum_cal", "pmcq vs. calibrated sum ZNA", nbins, xmin, xmax, nbins, xmin, xmax); 

  TProfile *hpmcZNC_vs_sum_tp = new TProfile("hpmcZNC_vs_sum_tp", "Before", nbins, xmin, xmax);
  TProfile *hpmcZNA_vs_sum_tp = new TProfile("hpmcZNA_vs_sum_tp", "Before", nbins, xmin, xmax);

  TProfile *hpmcZNC_vs_sum_cal_tp = new TProfile("hpmcZNC_vs_sum_cal_tp", "After", nbins, xmin, xmax);
  TProfile *hpmcZNA_vs_sum_cal_tp = new TProfile("hpmcZNA_vs_sum_cal_tp", "After", nbins, xmin, xmax);


  for(Int_t iev=0; iev<nentries; iev++){
    // So the loop is done for each entry of the tree. 
    // Because we want to check the towers and sum for each entry. 

    tree->GetEntry(iev);
    //
    bool isZNChit = true, isZNAhit = true;

    if (pmcZNC < kVeryNegative) {
      pmcZNC = kVeryNegative;
      isZNChit = false;
    }
    if (pmcZNA < kVeryNegative) {
      pmcZNA = kVeryNegative;
      isZNAhit = false;
    }

    if(pm1ZNC==0 || pm2ZNC==0 || pm3ZNC==0 || pm4ZNC==0 ) {
      // printf("One of your ZNC towers has zero energy... Thats weird. (%i) \n", iev );
      isZNChit = false;
      }
    if(pm1ZNA==0 || pm2ZNA==0 || pm3ZNA==0 || pm4ZNA==0 ) {
      // printf("One of your ZNA towers has zero energy... Thats weird. (%i) \n", iev );
      isZNAhit = false;
      }

    //
    auto sumZNC = kVeryNegative;
    auto sumZNA = kVeryNegative;
    //
    if(isZNChit) {
      float sumZNC = pm1ZNC + pm2ZNC + pm3ZNC + pm4ZNC; 
      hpmcZNC->Fill(pmcZNC);
      hsumZNC->Fill(sumZNC);
      hpmcZNC_vs_sum->Fill(pmcZNC, sumZNC);
      hpmcZNC_vs_sum_tp->Fill(pmcZNC, sumZNC-pmcZNC);
      //
      cumulate(0, pmcZNC, pm1ZNC, pm2ZNC, pm3ZNC, pm4ZNC, 1.);
    }
    if(isZNAhit) {
      float sumZNA = pm1ZNA + pm2ZNA + pm3ZNA + pm4ZNA; 
      hpmcZNA->Fill(pmcZNA);
      hsumZNA->Fill(sumZNA);
      hpmcZNA_vs_sum->Fill(pmcZNA, sumZNA);
      hpmcZNA_vs_sum_tp->Fill(pmcZNA, sumZNA-pmcZNA);
      //
      cumulate(1, pmcZNA, pm1ZNA, pm2ZNA, pm3ZNA, pm4ZNA, 1.);
    }

  } 
  // So here we used cumulate for all entries in the tree. 

  //
  for (int ih = 0; ih < ndet; ih++) {
    int ierr = minimization(ih);
    if (mSum[ih][5][5] >= 0.) {
      if (ierr) {
        printf(" -- FAILED processing data for detector %d\n\n", ih);
      } 
    } else {
      printf(" -- FAILED processing data for det. %d, reason might be TOO FEW EVENTS! [nev = %1.0f]", ih, mSum[ih][5][5]);
    }
  }


  for(Int_t iev=0; iev<nentries; iev++){
    // Loop again with the calibration constants! 

    tree->GetEntry(iev);
  
    bool isZNChit = true, isZNAhit = true;

    if (pmcZNC < kVeryNegative) {
      isZNChit = false;
    } 
 
    if (pmcZNA < kVeryNegative) {
      isZNAhit = false;
    }
    if(pm1ZNC==0 || pm2ZNC==0 || pm3ZNC==0 || pm4ZNC==0 ) {
      // printf("One of your ZNC towers has zero energy... Thats weird. (%i) \n", iev );
      isZNChit = false;
      }

    if(pm1ZNA==0 || pm2ZNA==0 || pm3ZNA==0 || pm4ZNA==0 ) {
      // printf("One of your ZNA towers has zero energy... Thats weird. (%i) \n", iev );
      isZNAhit = false;
      }

    if(isZNChit) {
      double sumqcorrC = pm1ZNC * (coeff[0][1]) + pm2ZNC * (coeff[0][2]) + pm3ZNC * (coeff[0][3]) + pm4ZNC * (coeff[0][4]);
      hsumZNC_cal->Fill(sumqcorrC);
      hpmcZNC_vs_sum_cal->Fill(pmcZNC, sumqcorrC);
      hpmcZNC_vs_sum_cal_tp->Fill(pmcZNC, sumqcorrC-pmcZNC);

    }
    if(isZNAhit) {
      float sumqcorrA = pm1ZNA * (coeff[1][1]) + pm2ZNA * (coeff[1][2]) + pm3ZNA * (coeff[1][3]) + pm4ZNA * (coeff[1][4]);
      hsumZNA_cal->Fill(sumqcorrA);
      hpmcZNA_vs_sum_cal->Fill(pmcZNA, sumqcorrA);
      hpmcZNA_vs_sum_cal_tp->Fill(pmcZNA, sumqcorrA-pmcZNA);
    }

  } 

  TCanvas* cZNC = plotTProfiles({hpmcZNC_vs_sum_tp, hpmcZNC_vs_sum_cal_tp}, "ZNC"); 
  TCanvas* cZNA = plotTProfiles({hpmcZNA_vs_sum_tp, hpmcZNA_vs_sum_cal_tp}, "ZNA"); 

  TFile *fout = new TFile("ZNic.root","RECREATE");
  fout->cd();

  hpmcZNC->Write();
  hpmcZNA->Write();
  hsumZNC->Write();
  hsumZNA->Write();

  hsumZNC_cal->Write();
  hsumZNA_cal->Write();

  hpmcZNC_vs_sum->SetOption("COLZ");
  hpmcZNA_vs_sum->SetOption("COLZ");
  hpmcZNC_vs_sum->Write(); 
  hpmcZNA_vs_sum->Write(); 
  hpmcZNC_vs_sum_tp->Write(); 
  hpmcZNA_vs_sum_tp->Write(); 

  hpmcZNC_vs_sum_cal->SetOption("COLZ");
  hpmcZNA_vs_sum_cal->SetOption("COLZ");
  hpmcZNC_vs_sum_cal->Write(); 
  hpmcZNA_vs_sum_cal->Write(); 
  hpmcZNC_vs_sum_cal_tp->Write(); 
  hpmcZNA_vs_sum_cal_tp->Write(); 

  cZNA->Write(); 
  cZNC->Write(); 


  ZNCcalibcoeff->Write(); 
  ZNAcalibcoeff->Write(); 

  fout->Close();


}

void cumulate(int ih, double tc, double t1, double t2, double t3, double t4, double w = 1)
{
    if ((ih<0 || ih>1)){
      return;
    }
    //mSum[ndet][npar][npar] = {0};

    double val[npar] = {0, 0, 0, 0, 0, 1};
    val[0] = tc;
    val[1] = t1;
    val[2] = t2;
    val[3] = t3;
    val[4] = t4;
    printf(" det.%d \t mSum \n ", ih);
    for (int32_t i = 0; i < npar; i++) { 
      for (int32_t j = i; j < npar; j++) {
        mSum[ih][i][j] += val[i] * val[j] * w; 
        printf(" %1.0f ", mSum[ih][i][j]);
      }
      
      printf(" \n");
    }
    printf(" \n");

}

static void fcn(int &np, double *gin, double &f, double *par, int iflag)  
{
    // Calculate chisquare
    Double_t chi = 0;
    // NOTE: changed <np to <= np to include np=4 because otherwise the 4th parameter is not determined.. 
    for (int32_t i = 0; i < npar; i++) {
      for (int32_t j = 0; j < npar; j++) {
        chi += (i == 0 ? par[i] : -par[i]) * (j == 0 ? par[j] : -par[j]) * mAdd[i][j];
      }
    }
    // printf("it changed.");

    if(deming) chi = chi / (1 + par[1] * par[1] + par[2] * par[2] + par[3] * par[3] + par[4] * par[4]);

    f = chi;
    // printf("chi_sq = %e   | np = %i \n", f, np); 
}

int minimization(int ih)
{

    printf("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= \n"); 
    printf("mAdd:\n");
   
    for (int i = 0; i < npar; i++) {
      for (int j = 0; j < npar; j++) {
        if (j < i) {
          mAdd[i][j] = mSum[ih][j][i];
        } else {
          mAdd[i][j] = mSum[ih][i][j];
        }
         printf("%.0f ", mAdd[i][j]);
      }
      printf("\n");
    }
    printf("\n");
    //

    printf("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= \n"); 

    double arglist[10];
    int ierflg = 0;

    // So this is done to fit, but what exactly are we fitting? 
    TMinuit *gMinuit = new TMinuit(npar);  
    gMinuit->SetFCN(fcn);
    
    //arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist , 1, ierflg);
    
    // Set starting values and step sizes for parameters
    // The common PM is not minimized...
    // calibration coefficient is forced to 1 and step is forced to zero
    gMinuit->mnparm(0, "c0",    1.,   0.,     1.,     1., ierflg);
    gMinuit->mnparm(1, "c1", start, step, lowbnd, uppbnd, ierflg);
    gMinuit->mnparm(2, "c2", start, step, lowbnd, uppbnd, ierflg);
    gMinuit->mnparm(3, "c3", start, step, lowbnd, uppbnd, ierflg);
    gMinuit->mnparm(4, "c4", start, step, lowbnd, uppbnd, ierflg);
    gMinuit->mnparm(5, "c5",    1.,   0.,     1.,     1., ierflg);

    // Minimization step
    //arglist[0] = 500;
    //arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist, 0, ierflg);

    
    for ( int i = 0; i < npar; i++) {
      gMinuit->GetParameter(i, coeff[ih][i], coefferr[ih][i]);
      if(i>0 && i<5){
        if(ih==0)  {
          ZNCcalibcoeff->SetBinContent(i,coeff[ih][i]);
          ZNCcalibcoeff->SetBinError(i,coefferr[ih][i]);
          }
        else if (ih==1) {
          ZNAcalibcoeff->SetBinContent(i,coeff[ih][i]);
          ZNAcalibcoeff->SetBinError(i,coefferr[ih][i]);
          }
      }
    }
    //
    if(ih==0) printf("\n\t +++++++++++++++++++ \t ZNC \t +++++++++++++++++++\n");
    else if(ih==1) printf("\n\t +++++++++++++++++++ \t ZNA \t +++++++++++++++++++\n");
    printf("\n\t _____________________ CALIBRATION COEFFICIENTS _____________________\n");
    printf("\t pmq1 %f  pmq2 %f  pmq3 %f  pmq4 %f \n", coeff[ih][1], coeff[ih][2], coeff[ih][3], coeff[ih][4]);
    printf("\t\t _____________________ ERRORS _____________________\n");
    printf("\t pmq1 %f  pmq2 %f  pmq3 %f  pmq4 %f \n\n", coefferr[ih][1], coefferr[ih][2], coefferr[ih][3], coefferr[ih][4]);

    
    
    // OFFSET NON IMPLEMENTED FOR THE TIME BEING!
    // see $O2/Detectors/ZDC/calib/src/InterCalib.cxx to see how it is done on raw data

     return ierflg;
}

// **********************************************************************************
void setCosm() 
{
   // Set gStyle
   int font = 42;
   // From plain
   gStyle->SetFrameBorderMode(0);
   gStyle->SetFrameFillColor(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetPadColor(10);
   gStyle->SetCanvasColor(10);
   gStyle->SetTitleFillColor(10);
   gStyle->SetTitleBorderSize(0);
   gStyle->SetStatColor(10);
   gStyle->SetStatBorderSize(1);
   gStyle->SetLegendBorderSize(0);
   //
   gStyle->SetDrawBorder(0);
   gStyle->SetTextFont(font);
   gStyle->SetStatFont(font);
   gStyle->SetStatFontSize(0.048);
   gStyle->SetStatX(0.92);
   gStyle->SetStatY(0.92);
   gStyle->SetStatH(0.03);
   gStyle->SetStatW(0.28);
   gStyle->SetTickLength(0.02,"y");
   gStyle->SetEndErrorSize(3);
   gStyle->SetLabelSize(0.038,"xyz");
   gStyle->SetLabelFont(font,"xyz");
   gStyle->SetLabelOffset(0.01,"xyz");
   gStyle->SetTitleFont(font,"xyz");
   gStyle->SetTitleOffset(1.2,"xyz");
   gStyle->SetTitleSize(0.045,"xyz");
   gStyle->SetPadGridX(true);
   gStyle->SetPadGridY(true);
   gStyle->SetMarkerSize(1);
   gStyle->SetPalette(kPastel);
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(1111);
   gStyle->SetOptFit(111);
   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadRightMargin(0.07);
   gStyle->SetPadTopMargin(0.08);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetEndErrorSize(0);
   gStyle->SetAxisMaxDigits(2); 
   gStyle->SetLineWidth(2);
}

TCanvas* plotTProfiles(std::vector<TProfile*> profiles, const char* name){
  setCosm(); 
  // Function to plot TProfiles on canvas and ratio to linear function> 
  TCanvas* c1Data = new TCanvas(TString::Format("TProfiles Before and after calibration %s", name), TString::Format("TProfiles Before and after calibration %s", name), 800, 800);
  TLegend* leg = new TLegend(0.20, 0.72, 0.51, 0.83);
  leg->SetTextSize(28);

  c1Data->cd(); 
  int i=0; 
  for(TProfile* profile : profiles){ 
    if(i==0){ 
      profile->GetXaxis()->SetTitle("Energy pmc");
      profile->GetYaxis()->SetTitle("Sum-pmc");
      profile->GetYaxis()->SetTitleOffset(1.9);
      profile->GetYaxis()->SetMaxDigits(3);
      profile->GetXaxis()->SetMaxDigits(3);
      profile->Draw("PLC PMC"); 
    } 
    else{
      profile->Draw("PLC PMC SAME");
      }
      leg->AddEntry(profile, TString::Format("%s", profile->GetTitle()), "lp"); 
      i++;  
  }
  leg->Draw(); 

  return c1Data; 
}