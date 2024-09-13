// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief Task for ZDC tower inter-calibration
/// \author chiara.oppedisano@cern.ch
// Minimal example to run this task:
// export OPTIONS="-b --configuration json://config.json --aod-file AO2D.root"
//  o2-analysis-timestamp ${OPTIONS} | 
//  o2-analysis-event-selection ${OPTIONS} |  
//  o2-analysis-track-propagation ${OPTIONS} |  
//  o2-analysis-mm-zdc-task-intercalib ${OPTIONS} 

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "ZDCInterCalib.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TMinuit.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::evsel;

constexpr double kVeryNegative = -1.e12;
constexpr int ndet = 2;
constexpr int npar = 6;
static double mSum[ndet][npar][npar];
static double mAdd[npar][npar];

using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using ColEvSels = soa::Join<aod::Collisions, aod::EvSels>;

struct zdcInterCalib {
  
  Produces<aod::ZDCInterCalib> zTab;

  // Configurable parameters
  Configurable<double> lowbnd{"lowbnd", 0.1, "lower bound for minimization"};
  Configurable<double> uppbnd{"uppbnd", 10., "upper bound for minimization"};
  Configurable<double> start{"start", 0.8, "start value for minimization"};
  Configurable<double> step{"step", 0.1, "step value for minimization"};
  //
  Configurable<int> nBins{"nBins", 400, "n bins"};
  Configurable<float> MaxZN{"MaxZN", 399.5, "Max ZN signal"};
  //
  HistogramRegistry registry{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    registry.add("ZNApmc", "ZNApmc; ZNA PMC; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNCpmc", "ZNCpmc; ZNC PMC; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNApm1", "ZNApm1; ZNA PM1; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNApm2", "ZNApm2; ZNA PM2; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNApm3", "ZNApm3; ZNA PM3; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNApm4", "ZNApm4; ZNA PM4; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNCpm1", "ZNCpm1; ZNC PM1; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNCpm2", "ZNCpm2; ZNC PM2; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNCpm3", "ZNCpm3; ZNC PM3; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNCpm4", "ZNCpm4; ZNC PM4; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNAsumq", "ZNAsumq; ZNA uncalib. sum PMQ; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
    registry.add("ZNCsumq", "ZNCsumq; ZNC uncalib. sum PMQ; Entries", {HistType::kTH1F, {{nBins, -0.5, MaxZN}}});
  }

  void process(ColEvSels const& cols, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcs*/)
  {
    // collision-based event selection
    for (auto& collision : cols) {
      const auto& foundBC = collision.foundBC_as<BCsRun3>();
      if (foundBC.has_zdc()) {
        const auto& zdc = foundBC.zdc();

        // To assure that ZN have a genuine signal (tagged by the relative TDC)
        // we can check that the amplitude is >0 or that ADC is NOT very negative (-inf)
        // If this is not the case, signals are set to kVeryNegative values
        double pmcZNC = zdc.energyCommonZNC();
        double pmcZNA = zdc.energyCommonZNA();
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
        double pmqZNC[4] = {0, 0, 0, 0,};
        double pmqZNA[4] = {0, 0, 0, 0,};
        //
        if(isZNChit) {
          for(int it=0; it<4; it++) {
            pmqZNC[it] = (zdc.energySectorZNC())[it];
            sumZNC += pmqZNC[it];
          }
          registry.get<TH1>(HIST("ZNCpmc"))->Fill(pmcZNC);
          registry.get<TH1>(HIST("ZNCpm1"))->Fill(pmqZNC[0]);
          registry.get<TH1>(HIST("ZNCpm2"))->Fill(pmqZNC[1]);
          registry.get<TH1>(HIST("ZNCpm3"))->Fill(pmqZNC[2]);
          registry.get<TH1>(HIST("ZNCpm4"))->Fill(pmqZNC[3]);
          registry.get<TH1>(HIST("ZNCsumq"))->Fill(sumZNC);
          //
          cumulate(0, pmcZNC, pmqZNC[0], pmqZNC[1], pmqZNC[2], pmqZNC[3], 1.);
        }
        if(isZNAhit) {
          for(int it=0; it<4; it++) {
            pmqZNA[it] = (zdc.energySectorZNA())[it];
            sumZNA += pmqZNA[it];
          }
          //
          registry.get<TH1>(HIST("ZNApmc"))->Fill(pmcZNA);
          registry.get<TH1>(HIST("ZNApm1"))->Fill(pmqZNA[0]);
          registry.get<TH1>(HIST("ZNApm2"))->Fill(pmqZNA[1]);
          registry.get<TH1>(HIST("ZNApm3"))->Fill(pmqZNA[2]);
          registry.get<TH1>(HIST("ZNApm4"))->Fill(pmqZNA[3]);
          registry.get<TH1>(HIST("ZNAsumq"))->Fill(sumZNA);
          //
          cumulate(1, pmcZNA, pmqZNA[0], pmqZNA[1], pmqZNA[2], pmqZNA[3], 1.);
        }
        if (isZNAhit || isZNChit) (pmcZNA, pmqZNA[0], pmqZNA[1], pmqZNA[2], pmqZNA[3], pmcZNC, pmqZNC[0], pmqZNC[1], pmqZNC[2], pmqZNC[3]);
      }
    }
    //
    /*for (int ih = 0; ih < ndet; ih++) {
      int ierr = minimization(ih);
      if (mSum[ih][5][5] >= 0.) {
        if (ierr) {
          printf(" -- FAILED processing data for detector %d\n\n", ih);
        } else {
          printf(" +++ Calculating intercalibration coeff. for det. %d", ih);
        }
      } else {
        printf(" -- FAILED processing data for det. %d, reason might be TOO FEW EVENTS! [nev = %1.0f]", ih, mSum[ih][5][5]);
      }
    }*/
  }
  

  void cumulate(int ih, double tc, double t1, double t2, double t3, double t4, double w = 1)
  {

    if ((ih<0 || ih>1)){
      return;
    }

    double val[npar] = {0, 0, 0, 0, 0, 1};
    val[0] = tc;
    val[1] = t1;
    val[2] = t2;
    val[3] = t3;
    val[4] = t4;
    //printf(" det.%d  mSum \n ", ih);
    for (int32_t i = 0; i < npar; i++) {
      for (int32_t j = i; j < npar; j++) {
        mSum[ih][i][j] += val[i] * val[j] * w;
        //printf(" %1.0f ", mSum[ih][i][j]);
      }
      //printf(" \n");
    }
    //printf(" \n");

  }

  static void fcn(int &npar, double *gin, double &f, double *par, int iflag)  
  {
      //calculate chisquare
      Double_t chi = 0;
      for (int32_t i = 0; i < npar; i++) {
        for (int32_t j = 0; j < npar; j++) {
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
    
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist , 1, ierflg);
    
    // We introduce the calibration of the common PM in separate workflows
    // Calibration cvoefficient is forced to and step is forced to zero
    gMinuit->mnparm(0, "c0", 1., 0., 1., 1., ierflg);
    //
    gMinuit->mnparm(1, "c1", start, step, lowbnd, uppbnd, ierflg);
    gMinuit->mnparm(2, "c2", start, step, lowbnd, uppbnd, ierflg);
    gMinuit->mnparm(3, "c3", start, step, lowbnd, uppbnd, ierflg);
    gMinuit->mnparm(4, "c4", start, step, lowbnd, uppbnd, ierflg);

    double coeff[npar], coefferr[npar];
    for ( int i = 0; i < npar; i++) {
      gMinuit->GetParameter(i, coeff[i], coefferr[i]);
      //if(ih==0)  registry.get<TH1>(HIST("ZNCcalibcoeff"))->Fill(coeff[i]);
      //else if (ih==1) registry.get<TH1>(HIST("ZNAcalibcoeff"))->Fill(coeff[i]);
    }
    //
    printf("\n\t _________________________ CALIBRATION COEFFICIENTS _________________________\n");
    if(ih==0) printf("\t - ZNC:  pmq1 %f  pmq2 %f  pmq3 %f  pmq4 %f \n", coeff[1], coeff[2], coeff[3], coeff[4]);
    else if(ih==1) printf("\t - ZNA:  pmq1 %f  pmq2 %f  pmq3 %f  pmq4 %f \n", coeff[1], coeff[2], coeff[3], coeff[4]);
    //double sumqcorr = c1 * x[1] + c2 * x[2] + c3 * x[3] + c4 * x[4];
    
    // OFFSET NON IMPLEMENTED FOR THE TIME BEING!
    // see $O2/Detectors/ZDC/calib/src/InterCalib.cxx to see how it is done on raw data

     return ierflg;
   }
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<zdcInterCalib>(cfgc) 
  };
}
