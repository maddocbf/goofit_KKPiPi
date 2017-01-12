//THIS PROGRAM GENERATES MONTECARLO DATA GIVEN AN AMPLITUDE MODEL

//ROOT
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>

// GooFit stuff
#include "goofit/Variable.h" 
#include "goofit/PDFs/PolynomialPdf.h" 
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/DP4Pdf.h"
#include "goofit/PDFs/TruthResolution_Aux.h" 
#include "goofit/PDFs/Tddp4Pdf.h"
#include <thrust/count.h>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <random>
#include <ctime>
#include <functional>
#include <mcbooster/functors/FlagAcceptReject.h>

using namespace std;

// Constants used in more than one PDF component. 
const fptype _mD0 = 1.8645; 
const fptype piPlusMass = 0.13957018;
const fptype piMinusMass = 0.13957018; 
const fptype kPlusMass = 0.493677; 
const fptype kMinusMass = .493677;

int main (int argc, char** argv) {
  DecayInfo_DP* DKKPP_DI = new DecayInfo_DP();
  DKKPP_DI->meson_radius =1.5;
  DKKPP_DI->particle_masses.push_back(_mD0);
  DKKPP_DI->particle_masses.push_back(piPlusMass);
  DKKPP_DI->particle_masses.push_back(piMinusMass);
  DKKPP_DI->particle_masses.push_back(kPlusMass);
  DKKPP_DI->particle_masses.push_back(kMinusMass);
 
  Variable* Kstar892M   =   new Variable("K892M", 0.89581);
  Variable* Kstar892W   =   new Variable("K892W", 0.0474); 

  Variable* K11270M  = new Variable("K1_1270M", 1.272);
  Variable* K11270W  = new Variable("K1_1270W", 0.09); 
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<SpinFactor*> SFK1P2Kstar;//K1(1270)+(Kstar0 pi+)K- 
  SFK1P2Kstar.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 3,0,1,2));

  //////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<Lineshape*> LSK1P2Kstar;
  LSK1P2Kstar.push_back( new Lineshape("K1(1270)+", K11270M, K11270W, 1, M_23_1, LS::BW, FF::BL2) );    
  LSK1P2Kstar.push_back( new Lineshape("Kstar(1430)", Kstar892M, Kstar892W, 1, M_23, LS::BW, FF::BL2) ); 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// the very last parameter means that we have two permutations. so the first half of the Lineshapes 
  // and the first half of the spinfactors are amplitude 1, rest is amplitude 2
  // This means that it is important for symmetrized amplitueds that the spinfactors and lineshapes are in the "right" order
  
  //RS Model
  Amplitude* AMP_K1P2Kstar       = new Amplitude( "K1P2Kstar",   new Variable("K1P2Kstar_real", 1.0),     new Variable("K1P2Kstar_imag", 0.0  ), LSK1P2Kstar, SFK1P2Kstar, 1);
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


  DKKPP_DI->amplitudes.push_back(AMP_K1P2Kstar);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////i////////////////////////////////////////////////////////////////////////////

  Variable* m12 = new Variable("m12", 0, 3);
  Variable* m34 = new Variable("m34", 0, 3); 
  Variable* cos12 = new Variable("cos12", -1, 1);
  Variable* cos34 = new Variable("cos34", -1, 1);
  Variable* phi = new Variable("phi", -3.5, 3.5);
  Variable* eventNumber = new Variable("eventNumber", 0, INT_MAX);

  Variable* constantOne = new Variable("constantOne", 1); 
  Variable* constantZero = new Variable("constantZero", 0);
 
  std::vector<Variable*> vars;
  vars.push_back(m12);
  vars.push_back(m34);
  vars.push_back(cos12);
  vars.push_back(cos34);
  vars.push_back(phi);
  vars.push_back(eventNumber); 

  UnbinnedDataSet currData(vars); 

  vector<Variable*> observables;
  vector<Variable*> coefficients; 
  vector<Variable*> offsets;

  observables.push_back(m12);
  observables.push_back(m34);
  observables.push_back(cos12);
  observables.push_back(cos34);
  observables.push_back(phi);
  observables.push_back(eventNumber);

  offsets.push_back(constantZero);
  offsets.push_back(constantZero);
  coefficients.push_back(constantOne); 

  fprintf(stderr, "I'm here zero"); 

  PolynomialPdf* eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);

  DPPdf* dp = new DPPdf("test_TD", observables, DKKPP_DI, eff, 1);
 
  fprintf(stderr,"I'm here one\n"); 

  TFile *file = new TFile( argv[4] , "RECREATE");
  TTree *tree = new TTree("events", "events");

  double tm12_2,tm34_2,tc12_2,tc34_2,tphi_2,tdtime_2;

  tree->Branch("m12",          &tm12_2,         "m12/D");
  tree->Branch("m34",          &tm34_2,         "m34/D");
  tree->Branch("c12",          &tc12_2,         "c12/D");
  tree->Branch("c34",          &tc34_2,         "c34/D");
  tree->Branch("phi",          &tphi_2,         "phi/D");

  fprintf(stderr, "I'm here two\n"); 

  //mcbooster::FlagAcceptReject(0,0);
  int generatedEvents = 0;
  int RunNum = 0;
  fprintf(stderr, "I'm here three\n");
  int BatchSize = strtoul(argv[1], NULL,0);
  fprintf(stderr, "I'm here three\n"); 
  unsigned int offi = strtoul(argv[3], NULL,0);
  unsigned int genEvts =strtoul(argv[2], NULL,0);

  double wmax = 0;
  //mcbooster::FlagAcceptReject FlagIt =1;// mcbooster::FlagAcceptReject(0.1,5);

  fprintf(stderr, "I'm here"); 
  while(generatedEvents < genEvts )
  {
    fprintf(stderr,"I'm here"); 
    unsigned int keptEvts = 0;
    dp->setGenerationOffset(offi);
    auto tuple = dp->GenerateSig(BatchSize);
    fprintf(stderr,"after gen\n");
    auto particles = std::get<0>(tuple);
    auto variables = std::get<1>(tuple);
    auto weights = std::get<2>(tuple);
    auto flags = std::get<3>(tuple);
    // int accepted = thrust::count_if(flags.begin(), flags.end(), thrust::identity<bool>());
    ++RunNum;
    // generatedEvents += accepted;
    fprintf(stderr,"after gen\n");

    for (int i = 0; i < weights.size(); ++i)
    {
      if (wmax<weights[i]) wmax = weights[i];
      if (generatedEvents < genEvts && flags[i]==1){
        ++generatedEvents;
        ++keptEvts;
        // printf("PF %i: %s %.5g %.5g %.5g %.5g %.5g %.5g\n",i, (bool)flags[i] ? "true" : "false", weights[i], (*(variables[0]))[i], (*(variables[1]))[i], (*(variables[2]))[i], (*(variables[3]))[i], (*(variables[4]))[i]);
        tm12_2 = (*(variables[0]))[i];
        tm34_2 = (*(variables[1]))[i];
        tc12_2 = (*(variables[2]))[i];
        tc34_2 = (*(variables[3]))[i];
        tphi_2 = (*(variables[4]))[i];
        tree->Fill();
        // printf("Buffer %i: %.5g %.5g %.5g %.5g %.5g %.5g \n",i, (*myweights)[i],(*Buffer_m12)[i], (*Buffer_m34)[i], (*Buffer_c12)[i], (*Buffer_c34)[i], (*Buffer_phi)[i], (*Buffer_dt)[i]);
      }
    }
    fprintf(stderr,"Run # %i:Using accept-reject method leaves you with %i out of %i events.  %.4g %% of Total offset: %u\n",RunNum, keptEvts, BatchSize, generatedEvents*100.0/genEvts, offi);
    offi += BatchSize;
    delete variables[0];
    delete variables[1];
    delete variables[2];
    delete variables[3];
    delete variables[4];

    delete particles[0];
    delete particles[1];
    delete particles[2];
    delete particles[3];
  }
  // printf("start\n");
  // int i = 0;
  // printf("Buffer %i: %.5g %.5g %.5g %.5g %.5g %.5g \n",i, (*myweights)[i],(*Buffer_m12)[i], (*Buffer_m34)[i], (*Buffer_c12)[i], (*Buffer_c34)[i], (*Buffer_phi)[i], (*Buffer_dt)[i]);

  // printf("start2\n");
  std::ofstream out;
  string outname ="Max_observed_weights.txt";
  out.open(outname.c_str(), std::ios::app);
  out.precision(10);

  out << wmax <<endl;

  tree->Write();
  file->Close();
  // printf("overall wmax %f, keept %u evts, reweight ratio %.5g\n",wmax, keptEvts, (double)keptEvts/genEvts );
  printf("%i\n",offi);
  //std::cout<<amps[2]<<std::endl;
  return 0; 

}
