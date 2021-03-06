//THIS PROGRAM GENERATES MONTECARLO DATA GIVEN AN AMPLITUDE MODEL


//ROOT
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>

// GooFit stuff
#include "goofit/Variable.h" 
#include "goofit/PDFs/basic/PolynomialPdf.h" 
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/physics/DP4Pdf.h"
#include "goofit/PDFs/physics/TruthResolution_Aux.h" 
#include "goofit/PDFs/physics/Tddp4Pdf.h"
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

  // cudaSetDevice(0);

  DecayInfo_DP* DKKPP_DI = new DecayInfo_DP();
  DKKPP_DI->meson_radius =5;
  DKKPP_DI->particle_masses.push_back(_mD0);
  DKKPP_DI->particle_masses.push_back(piPlusMass);
  DKKPP_DI->particle_masses.push_back(piMinusMass);
  DKKPP_DI->particle_masses.push_back(kPlusMass);
  DKKPP_DI->particle_masses.push_back(kMinusMass);
 
  Variable* RhoMass  =  new Variable("rho_mass", 0.77526);
  Variable* RhoWidth =  new Variable("rho_width", 0.1478); 
  Variable* Kstar892M   =   new Variable("K892M", 0.89581);
  Variable* Kstar892W   =   new Variable("K892W", 0.0474); 
  Variable* f600M  =    new Variable("f600M", 0.519);
  Variable* f600W  =    new Variable("f600W", 0.454); 
  Variable* a1M  =      new Variable("a1M", 1.237);
  Variable* a1W  =      new Variable("a1W", 0.526); 
  Variable* K11270M  = new Variable("K1_1270M", 1.28241);
  Variable* K11270W  = new Variable("K1_1270W", 0.06596); 
  Variable* K1430M  = new Variable("K0_1430M", 1.425);
  Variable* K1430W  = new Variable("K0_1430W", 0.27);

  Variable* Kstar1410M    = new Variable("K1410M", 1.414);
  Variable* Kstar1410W    = new Variable("K1410W", 0.232); 

  Variable* rho1450M  = new Variable("rho1450M", 1.465);
  Variable* rho1450W  = new Variable("rho1450W", 0.400); 

  Variable* K1460M    = new Variable("K1460M", 1.351);
  Variable* K1460W    = new Variable("K1460W", 0.281); 

  Variable* f0_1370M  = new Variable("f0_1370M", 1.350);
  Variable* f0_1370W  = new Variable("f0_1370W", 0.35); 

  Variable* K1_1400M  = new Variable("K1_1400M", 1.403);
  Variable* K1_1400W  = new Variable("K1_1400W", 0.174); 

  Variable* K2_1430M  = new Variable("K2_1430M", 1.4256);
  Variable* K2_1430W  = new Variable("K2_1430W", 0.0985); 
  
  Variable* phi1020M = new Variable("phi1020M", 0.1019); 
  Variable* phi1020W = new Variable("phi1020W", 0.0004); 

  std::vector<Variable*> LassVars;
  LassVars.push_back( new Variable("lass_a",2.07) );
  LassVars.push_back( new Variable("lass_r",3.32) );
  LassVars.push_back( new Variable("lass_pf",0.0) );
  LassVars.push_back( new Variable("lass_pr",0.0) );
  LassVars.push_back( new Variable("lass_F",1.0) );
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<SpinFactor*> SFK1P2Kstar;//K1(1270)+(Kstar0 pi+)K- 
  SFK1P2Kstar.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 3,0,1,2));

  std::vector<SpinFactor*> SFK1M2Kstar;//K1(1270)-(Kstar0bar pi-)K+ 
  SFK1M2Kstar.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 2,1,0,3));

  std::vector<SpinFactor*> SFK1P2Rho;//K1(1270)+(rho K+)K- 
  SFK1P2Rho.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 3,2,0,1)); 

  std::vector<SpinFactor*> SFK1M2Rho;//K1(1270)-(rho K-)K+ 
  SFK1M2Rho.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 2,3,1,0)); 

  std::vector<SpinFactor*> SFKstarP2Kstar;//Kstar(1410)+(Kstar0 pi+)K- 
  SFKstarP2Kstar.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 3,0,1,2));

  std::vector<SpinFactor*> SFKstarM2Kstar;//Kstar(1410)-(Kstar0 pi-)K+ 
  SFKstarM2Kstar.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4,2,1,0,3));  

  std::vector<SpinFactor*> SFKstarKstar;
  SFKstarKstar.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, 0,3,1,2));  

  std::vector<SpinFactor*> SFPhiRhoS; 
  SFPhiRhoS.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, 2,3,0,1)); 

  std::vector<SpinFactor*> SFPhiRhoD;
  SFPhiRhoD.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D, 2,3,0,1)); 

 // std::vector<SpinFactor*> SFPhipipi; 
 // SFPhipipi.push_back( new SpinFactor("SF", SF_4Body::DtoVP1P2_VtoP3P4, 0,1,2,3)); 

  std::vector<SpinFactor*> SFNonRes1;
  SFNonRes1.push_back( new SpinFactor("SF", SF_4Body::ONE, 0,1,2,3)); 

  std::vector<SpinFactor*> SFNonRes2; 
  SFNonRes2.push_back( new SpinFactor("SF", SF_4Body::ONE, 0,3,1,2));   

  //////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<Lineshape*> LSK1P2Kstar;
  LSK1P2Kstar.push_back( new Lineshape("K1(1270)+", K11270M, K11270W, 2, M_23_4, LS::BW, FF::BL2) );    
  LSK1P2Kstar.push_back( new Lineshape("Kstar(1410)", Kstar1410M, Kstar1410W, 2, M_23, LS::BW, FF::BL2) ); 

  std::vector<Lineshape*> LSK1M2Kstar; 
  LSK1M2Kstar.push_back( new Lineshape("K1(1270)-", K11270M, K11270W, 2, M_14_2, LS::BW, FF::BL2) ); 
  LSK1M2Kstar.push_back( new Lineshape("Kstar(1410)", Kstar1410M, Kstar1410W,1, M_14, LS::BW, FF::BL2) );  

  std::vector<Lineshape*> LSK1P2Rho;
  LSK1P2Rho.push_back( new Lineshape("K1(1270)+", K11270M, K11270W, 1, M_12_3, LS::BW, FF::BL2) );
  LSK1P2Rho.push_back( new Lineshape("Rho ", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) );   

  std::vector<Lineshape*> LSK1M2Rho; 
  LSK1M2Rho.push_back( new Lineshape("K1(1270)-", K11270M, K11270W, 1, M_12_4, LS::BW, FF::BL2) ); 
  LSK1M2Rho.push_back( new Lineshape("Rho ", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) ); 
 
  std::vector<Lineshape*> LSKstarP2Kstar;
  LSKstarP2Kstar.push_back( new Lineshape("KstarP", Kstar1410M, Kstar1410W, 1, M_23_1, LS::BW, FF::BL2) ); 
  LSKstarP2Kstar.push_back( new Lineshape("Kstar", Kstar892M, Kstar892W, 1, M_23, LS::BW, FF::BL2) );  

  std::vector<Lineshape*> LSKstarM2Kstar; 
  LSKstarM2Kstar.push_back( new Lineshape("KstarM", Kstar1410M, Kstar1410W, 1, M_14_2, LS::BW, FF::BL2) ); 
  LSKstarM2Kstar.push_back( new Lineshape("Kstar", Kstar892M, Kstar892W, 1, M_14, LS::BW, FF::BL2) ); 
 
  std::vector<Lineshape*> LSKstarKstarbar; 
  LSKstarKstarbar.push_back( new Lineshape("Kstar", Kstar892M, Kstar892W, 1, M_23, LS::BW, FF::BL2) ); 
  LSKstarKstarbar.push_back( new Lineshape("Kstarbar", Kstar892M, Kstar892W, 1, M_14, LS::BW, FF::BL2) ); 
  
  std::vector<Lineshape*> LSPhiRhoS; 
  LSPhiRhoS.push_back( new Lineshape("Phi ", phi1020M, phi1020W, 1, M_34, LS::BW, FF::BL2) ); 
  LSPhiRhoS.push_back( new Lineshape("Rho ", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) ); 


  std::vector<Lineshape*> LSPhiRhoD; 
  LSPhiRhoD.push_back( new Lineshape("Phi ", phi1020M, phi1020W, 1, M_34, LS::BW, FF::BL2) ); 
  LSPhiRhoD.push_back( new Lineshape("Rho ", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) ); 

 // std::vector<Lineshape*> LSPhipipi; 
 // LSPhipipi.push_back( new Lineshape("Phi ",phi1020M,phi1020W, 1, M_34, LS::BW, FF::BL2) ); 
 // LSPhipipi.push_back( new Lineshape("pipi ", new Variable("NR5", 0.0), new Variable("NR6",0.0), 1, M_12, LS::BW, FF::BL2) ); 

  std::vector<Lineshape*> LSNonRes1; 
  LSNonRes1.push_back( new Lineshape("NonRes1 ", new Variable("NR1", 0.0), new Variable("NR2", 0.0),1, M_34_2, LS::nonRes, FF::BL2) ); 

  std::vector<Lineshape*> LSNonRes2; 
  LSNonRes2.push_back( new Lineshape("NonRes2 ", new Variable("NR3", 0.0), new Variable("NR4", 0.0),1, M_34_2, LS::nonRes, FF::BL2) ); 
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




// the very last parameter means that we have two permutations. so the first half of the Lineshapes 
  // and the first half of the spinfactors are amplitude 1, rest is amplitude 2
  // This means that it is important for symmetrized amplitueds that the spinfactors and lineshapes are in the "right" order
  
  //RS Model
  Amplitude* AMP_K1P2Kstar       = new Amplitude( "K892_rho770_S",   new Variable("K892_rho770_S_real", 1.0),     new Variable("K892_rho770_S_imag", 0.0), LSK1P2Kstar, SFK1P2Kstar, 1);

  Amplitude* AMP_K1M2Kstar       = new Amplitude( "K892_rho770_P",   new Variable("K892_rho770_P_real",   1429.17),   new Variable("K892_rho770_P_imag", -446.544), LSK1M2Kstar, SFK1M2Kstar , 1);

  Amplitude* AMP_K1P2Rho       = new Amplitude( "K892_rho770_D",   new Variable("K892_rho770_D_real",    -55.6514), new Variable("K892_rho770_D_imag",24.9602), LSK1P2Rho, SFK1P2Rho, 1);

  Amplitude* AMP_K1M2Rho      = new Amplitude( "K1410_rho770",    new Variable("K1410_rho770_P_real",  -127.111),  new Variable("K1410_rho770_P_imag",18.9476), LSK1M2Rho, SFK1M2Rho, 1);
 
 Amplitude* AMP_KstarP2Kstar         = new Amplitude( "K892_f0600",      new Variable("K892_f0600_real",      11.7278),  new Variable("K892_f0600_imag",-20.2763  ), LSKstarP2Kstar, SFKstarP2Kstar, 1);
 
 Amplitude* AMP_KstarM2Kstar     = new Amplitude( "rho1450_K0_1430", new Variable("rho1450_K0_1430_real", 9.06177),   new Variable("rho1450_K0_1430_imag", 6.99855 ), LSKstarM2Kstar  , SFKstarM2Kstar , 1);
 
 Amplitude* AMP_KstarKstarbar          = new Amplitude( "K1460_K892",      new Variable("K1460_K892_real",      -0.507829),  new Variable("K1460_K892_imag",-0.421295), LSKstarKstarbar  , SFKstarKstar , 1);
 
 Amplitude* AMP_PhiRhoS       = new Amplitude( "K1460_f0_1370",   new Variable("K1460_f0_1370_real",   1336.66),  new Variable("K1460_f0_1370_imag", 273.186 ), LSPhiRhoS  , SFPhiRhoS , 1);


  Amplitude* AMP_PhiRhoD        = new Amplitude( "K1_1270_K892",    new Variable("K1_1270_K892_real",  3760.8),   new Variable("K1_1270_K892_imag",-246.697 ), LSPhiRhoD  , SFPhiRhoD , 1);
 
 Amplitude* AMP_NonRes1      = new Amplitude( "K1_1270_rho770",  new Variable("K1_1270_rho770_real", 0.0723525),  new Variable("K1_1270_rho770_imag", 0.0673346 ), LSNonRes1  , SFNonRes1 , 1);

  Amplitude* AMP_NonRes2     = new Amplitude( "K1_1270_K0_1430", new Variable("K1_1270_K0_1430_real",-0.0496293 ),  new Variable("K1_1270_K0_1430_imag", 0.0673346 ), LSNonRes2  , SFNonRes2 , 1);


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


  DKKPP_DI->amplitudes.push_back(AMP_K1P2Kstar);
  DKKPP_DI->amplitudes.push_back(AMP_K1M2Kstar);
  DKKPP_DI->amplitudes.push_back(AMP_K1P2Rho);
  DKKPP_DI->amplitudes.push_back(AMP_K1M2Rho);
  DKKPP_DI->amplitudes.push_back(AMP_KstarP2Kstar);
  DKKPP_DI->amplitudes.push_back(AMP_KstarM2Kstar);
  DKKPP_DI->amplitudes.push_back(AMP_KstarKstarbar);
  DKKPP_DI->amplitudes.push_back(AMP_PhiRhoS);
  DKKPP_DI->amplitudes.push_back(AMP_PhiRhoD);
  DKKPP_DI->amplitudes.push_back(AMP_NonRes1);
  DKKPP_DI->amplitudes.push_back(AMP_NonRes2);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  DKKPP_DI->_tau = new Variable("tau", 0.4101);
  DKKPP_DI->_xmixing = new Variable("xmixing", 0.0049);
  DKKPP_DI->_ymixing = new Variable("ymixing", 0.0061);
  // DK3P_DI->_xmixing = new Variable("xmixing", 0, 0.00001, -.15, .15);
  // DK3P_DI->_ymixing = new Variable("ymixing", 0, 0.00001, -.15, .15);
  DKKPP_DI->_SqWStoRSrate = new Variable("SqWStoRSrate", 1.0/sqrt(300.0));  


  Variable* m12 = new Variable("m12", 0, 3);
  Variable* m34 = new Variable("m34", 0, 3); 
  Variable* cos12 = new Variable("cos12", -1, 1);
  Variable* cos34 = new Variable("m12", -1, 1);
  Variable* phi = new Variable("phi", -3.5, 3.5);
  Variable* eventNumber = new Variable("eventNumber", 0, INT_MAX);
  //Variable* dtime = new Variable("dtime", 0, 10);
  //Variable* sigmat = new Variable("sigmat",-3,3);
  Variable* constantOne = new Variable("constantOne", 1); 
  Variable* constantZero = new Variable("constantZero", 0);
 

  std::vector<Variable*> vars;
  vars.push_back(m12);
  vars.push_back(m34);
  vars.push_back(cos12);
  vars.push_back(cos34);
  vars.push_back(phi);
  vars.push_back(eventNumber); 
  //vars.push_back(dtime); 
  //vars.push_back(sigmat); 
  UnbinnedDataSet currData(vars); 

 
  *//DKKPP_DI->_xmixing = strtof(argv[5], NULL);
  *//DKKPP_DI->_ymixing = strtof(argv[6], NULL);

  vector<Variable*> observables;
  vector<Variable*> coefficients; 
  vector<Variable*> offsets;

  observables.push_back(m12);
  observables.push_back(m34);
  observables.push_back(cos12);
  observables.push_back(cos34);
  observables.push_back(phi);
  observables.push_back(eventNumber);
  //observables.push_back(dtime);
 // observables.push_back(sigmat);
  offsets.push_back(constantZero);
  offsets.push_back(constantZero);
  coefficients.push_back(constantOne); 
  fprintf(stderr, "I'm here zero"); 
  TruthResolution* dat = new TruthResolution();
  PolynomialPdf* eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
  DPPdf* dp = new DPPdf("test_TD", observables, DKKPP_DI, eff, 1);
 
 // dp->setGenDecayTimeLimit(0,3.5); // this corresponds to rougly 97% of the exponential. So this should be ok. And speeds up Generation significantly compared to [0,5]  

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
    fprintf(stderr,"Run # %i: x=%.6g y=%.6g Using accept-reject method leaves you with %i out of %i events.  %.4g %% of Total offset: %u\n",RunNum, *DKKPP_DI->_xmixing, *DKKPP_DI->_ymixing, keptEvts, BatchSize, generatedEvents*100.0/genEvts, offi);
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
  return 0; 

}
