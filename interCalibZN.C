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

constexpr double kVeryNegative = -1.e12;
//
constexpr int ndet = 2;
constexpr int npar = 6;
//
static double mSum[ndet][npar][npar];
static double mAdd[npar][npar];
//
const float xmin = -0.5;
const float xmax = 399.5;
//
const float lowbnd = 0.1;
const float uppbnd = 10.;
const float start = 0.8; 
const float step = 0.1;

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

  TH1F *hpmcZNC = new TH1F("hpmcZNC"," PMC ZNC", 400, xmin, xmax);
  TH1F *hpmcZNA = new TH1F("hpmcZNA"," PMC ZNA", 400, xmin, xmax);
  TH1F *hsumZNC = new TH1F("hsumZNC"," uncalibrated SUMq ZNC", 400, xmin, xmax);
  TH1F *hsumZNA = new TH1F("hsumZNA"," uncalibrated SUMq ZNA", 400, xmin, xmax);

  for(Int_t iev=0; iev<nentries; iev++){
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
    //
    auto sumZNC = kVeryNegative;
    auto sumZNA = kVeryNegative;
    //
    if(isZNChit) {
      float sumZNC = pm1ZNC + pm2ZNC + pm3ZNC + pm4ZNC; 
      hpmcZNC->Fill(pmcZNC);
      hsumZNC->Fill(sumZNC);
      //
      cumulate(0, pmcZNC, pm1ZNC, pm2ZNC, pm3ZNC, pm4ZNC, 1.);
    }
    if(isZNAhit) {
      float sumZNA = pm1ZNA + pm2ZNA + pm3ZNA + pm4ZNA; 
      hpmcZNA->Fill(pmcZNA);
      hsumZNA->Fill(sumZNA);
      //
      cumulate(1, pmcZNA, pm1ZNA, pm2ZNA, pm3ZNA, pm4ZNA, 1.);
    }

  }
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

  TFile *fout = new TFile("ZNic.root","RECREATE");
  fout->cd();

  hpmcZNC->Write();
  hpmcZNA->Write();
  hsumZNC->Write();
  hsumZNA->Write();

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
    for (int32_t i = 0; i < np; i++) {
      for (int32_t j = 0; j < np; j++) {
        chi += (i == 0 ? par[i] : -par[i]) * (j == 0 ? par[j] : -par[j]) * mAdd[i][j];
      }
    }
    f = chi;
}

int minimization(int ih)
{
    for (int i = 0; i < npar; i++) {
      for (int j = 0; j < npar; j++) {
        if (j < i) {
          mAdd[i][j] = mSum[ih][j][i];
        } else {
          mAdd[i][j] = mSum[ih][i][j];
        }
      }
    }
    //
    double arglist[10];
    int ierflg = 0;

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
    gMinuit->mnparm(0, "c5",    1.,   0.,     1.,     1., ierflg);

    // Minimization step
    //arglist[0] = 500;
    //arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist, 0, ierflg);

    double coeff[npar], coefferr[npar];
    for ( int i = 0; i < npar; i++) {
      gMinuit->GetParameter(i, coeff[i], coefferr[i]);
      //if(ih==0)  registry.get<TH1>(HIST("ZNCcalibcoeff"))->Fill(coeff[i]);
      //else if (ih==1) registry.get<TH1>(HIST("ZNAcalibcoeff"))->Fill(coeff[i]);
    }
    //
    if(ih==0) printf("\n\t +++++++++++++++++++ \t ZNC \t +++++++++++++++++++\n");
    else if(ih==1) printf("\n\t +++++++++++++++++++ \t ZNA \t +++++++++++++++++++\n");
    printf("\n\t _____________________ CALIBRATION COEFFICIENTS _____________________\n");
    printf("\t pmq1 %f  pmq2 %f  pmq3 %f  pmq4 %f \n", coeff[1], coeff[2], coeff[3], coeff[4]);
    printf("\t\t _____________________ ERRORS _____________________\n");
    printf("\t pmq1 %f  pmq2 %f  pmq3 %f  pmq4 %f \n\n", coefferr[1], coefferr[2], coefferr[3], coefferr[4]);
    
    //double sumqcorr = c1 * x[1] + c2 * x[2] + c3 * x[3] + c4 * x[4];
    
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
   gStyle->SetLegendBorderSize(1);
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
   gStyle->SetMarkerSize(1);
   gStyle->SetPalette(1,0);
   gStyle->SetOptTitle(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(111);
   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadRightMargin(0.02);
   gStyle->SetPadTopMargin(0.08);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetEndErrorSize(0);
}