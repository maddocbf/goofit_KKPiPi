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
const fptype KmMass = .493677;
int main (int argc, char** argv) {

  // cudaSetDevice(0);

  DecayInfo_DP* DK3P_DI = new DecayInfo_DP();
  DK3P_DI->meson_radius =5;
  DK3P_DI->particle_masses.push_back(_mD0);
  DK3P_DI->particle_masses.push_back(piPlusMass);
  DK3P_DI->particle_masses.push_back(piPlusMass);
  DK3P_DI->particle_masses.push_back(KmMass);
  DK3P_DI->particle_masses.push_back(piPlusMass);
 
  Variable* RhoMass  =  new Variable("rho_mass", 0.77526);
  Variable* RhoWidth =  new Variable("rho_width", 0.1478); 
  Variable* K892M   =   new Variable("K892M", 0.89581);
  Variable* K892W   =   new Variable("K892W", 0.0474); 
  Variable* f600M  =    new Variable("f600M", 0.519);
  Variable* f600W  =    new Variable("f600W", 0.454); 
  Variable* a1M  =      new Variable("a1M", 1.237);
  Variable* a1W  =      new Variable("a1W", 0.526); 
  Variable* K1_1270M  = new Variable("K1_1270M", 1.28241);
  Variable* K1_1270W  = new Variable("K1_1270W", 0.06596); 
  Variable* K0_1430M  = new Variable("K0_1430M", 1.425);
  Variable* K0_1430W  = new Variable("K0_1430W", 0.27);

  Variable* K1410M    = new Variable("K1410M", 1.414);
  Variable* K1410W    = new Variable("K1410W", 0.232); 

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
  
  std::vector<Variable*> LassVars;
  LassVars.push_back( new Variable("lass_a",2.07) );
  LassVars.push_back( new Variable("lass_r",3.32) );
  LassVars.push_back( new Variable("lass_pf",0.0) );
  LassVars.push_back( new Variable("lass_pr",0.0) );
  LassVars.push_back( new Variable("lass_F",1.0) );

 //Spin factors: we have two due to the bose symmetrization of the two pi+
  std::vector<SpinFactor*> SF_K892_rho770_S;
  SF_K892_rho770_S.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, 0, 1, 2, 3) );
  SF_K892_rho770_S.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, 3, 1, 2, 0) );
  //Lineshapes, also for both pi+ configurations
  std::vector<Lineshape*> LS_K892_rho770_S;
  LS_K892_rho770_S.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_34, LS::BW, FF::BL2) );
  LS_K892_rho770_S.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) );
  LS_K892_rho770_S.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_13, LS::BW, FF::BL2) );
  LS_K892_rho770_S.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS::BW, FF::BL2) );

  std::vector<SpinFactor*> SF_K892_rho770_P;
  SF_K892_rho770_P.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P, 0, 1, 2, 3) );
  SF_K892_rho770_P.push_back( new SpinFactor("SF", SF_4Body::FF_12_34_L1, 0, 1, 2, 3) );
  SF_K892_rho770_P.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P, 3, 1, 2, 0) );
  SF_K892_rho770_P.push_back( new SpinFactor("SF", SF_4Body::FF_12_34_L1, 3, 1, 2, 0) );
  std::vector<Lineshape*> LS_K892_rho770_P;
  LS_K892_rho770_P.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) );
  LS_K892_rho770_P.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_34, LS::BW, FF::BL2) );
  LS_K892_rho770_P.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS::BW, FF::BL2) );
  LS_K892_rho770_P.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_13, LS::BW, FF::BL2) );

  std::vector<SpinFactor*> SF_K892_rho770_D;
  SF_K892_rho770_D.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D, 0, 1, 2, 3) );
  SF_K892_rho770_D.push_back( new SpinFactor("SF", SF_4Body::FF_12_34_L2, 0, 1, 2, 3) );
  SF_K892_rho770_D.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D, 3, 1, 2, 0) );
  SF_K892_rho770_D.push_back( new SpinFactor("SF", SF_4Body::FF_12_34_L2, 3, 1, 2, 0) );
  std::vector<Lineshape*> LS_K892_rho770_D;
  LS_K892_rho770_D.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) );
  LS_K892_rho770_D.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_34, LS::BW, FF::BL2) );
  LS_K892_rho770_D.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS::BW, FF::BL2) );
  LS_K892_rho770_D.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_13, LS::BW, FF::BL2) );

  std::vector<SpinFactor*> SF_K1410_rho770_S;
  SF_K1410_rho770_S.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, 0, 1, 2, 3) );
  SF_K1410_rho770_S.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, 3, 1, 2, 0) );
  std::vector<Lineshape*> LS_K1410_rho770_S;
  LS_K1410_rho770_S.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) );
  LS_K1410_rho770_S.push_back( new Lineshape("K*(1410)", K1410M, K1410W, 1, M_34, LS::BW, FF::BL2) );
  LS_K1410_rho770_S.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS::BW, FF::BL2) );
  LS_K1410_rho770_S.push_back( new Lineshape("K*(1410)", K1410M, K1410W, 1, M_13, LS::BW, FF::BL2) );


  std::vector<SpinFactor*> SF_K1410_rho770_P;
  SF_K1410_rho770_P.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P, 0, 1, 2, 3) );
  SF_K1410_rho770_P.push_back( new SpinFactor("SF", SF_4Body::FF_12_34_L1, 0, 1, 2, 3) );
  SF_K1410_rho770_P.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P, 3, 1, 2, 0) );
  SF_K1410_rho770_P.push_back( new SpinFactor("SF", SF_4Body::FF_12_34_L1, 3, 1, 2, 0) );
  std::vector<Lineshape*> LS_K1410_rho770_P;
  LS_K1410_rho770_P.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) );
  LS_K1410_rho770_P.push_back( new Lineshape("K*(1410)", K1410M, K1410W, 1, M_34, LS::BW, FF::BL2) );
  LS_K1410_rho770_P.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS::BW, FF::BL2) );
  LS_K1410_rho770_P.push_back( new Lineshape("K*(1410)", K1410M, K1410W, 1, M_13, LS::BW, FF::BL2) );

  std::vector<SpinFactor*> SF_K892_f0_600;
  SF_K892_f0_600.push_back( new SpinFactor("SF", SF_4Body::DtoVS_VtoP1P2_StoP3P4, 2, 3, 0, 1) );
  SF_K892_f0_600.push_back( new SpinFactor("SF", SF_4Body::FF_12_34_L1, 2, 3, 0, 1) );
  SF_K892_f0_600.push_back( new SpinFactor("SF", SF_4Body::DtoVS_VtoP1P2_StoP3P4, 2, 0, 3, 1) );
  SF_K892_f0_600.push_back( new SpinFactor("SF", SF_4Body::FF_12_34_L1, 2, 0, 3, 1) );
  std::vector<Lineshape*> LS_K892_f0_600;
  LS_K892_f0_600.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_34, LS::BW, FF::BL2) );
  LS_K892_f0_600.push_back( new Lineshape("f600", f600M, f600W, 0, M_12, LS::Bugg3, FF::BL2) );
  LS_K892_f0_600.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_13, LS::BW, FF::BL2) );
  LS_K892_f0_600.push_back( new Lineshape("f600", f600M, f600W, 0, M_24, LS::Bugg3, FF::BL2) );

  std::vector<SpinFactor*> SF_rho1450_K0_1430;
  SF_rho1450_K0_1430.push_back( new SpinFactor("SF", SF_4Body::DtoVS_VtoP1P2_StoP3P4, 0, 1, 2, 3) );
  SF_rho1450_K0_1430.push_back( new SpinFactor("SF", SF_4Body::FF_12_34_L1, 0, 1, 2, 3) );
  SF_rho1450_K0_1430.push_back( new SpinFactor("SF", SF_4Body::DtoVS_VtoP1P2_StoP3P4, 3, 1, 2, 0) );
  SF_rho1450_K0_1430.push_back( new SpinFactor("SF", SF_4Body::FF_12_34_L1, 3, 1, 2, 0) );
  std::vector<Lineshape*> LS_rho1450_K0_1430;
  LS_rho1450_K0_1430.push_back( new Lineshape("rho(1450)", rho1450M, rho1450W, 1, M_12, LS::BW, FF::BL2) );
  LS_rho1450_K0_1430.push_back( new Lineshape("K(0)*(1430)", K0_1430M, K0_1430W, 0, M_34, LS::Lass_M3, FF::BL2, 1.5, LassVars) );
  LS_rho1450_K0_1430.push_back( new Lineshape("rho(1450)", rho1450M, rho1450W, 1, M_24, LS::BW, FF::BL2) );
  LS_rho1450_K0_1430.push_back( new Lineshape("K(0)*(1430)", K0_1430M, K0_1430W, 0, M_13, LS::Lass_M3, FF::BL2, 1.5, LassVars) );

  std::vector<SpinFactor*> SF_K1460_K892;
  SF_K1460_K892.push_back( new SpinFactor("SF", SF_4Body::DtoPP1_PtoVP2_VtoP3P4, 0, 1, 2, 3) );
  SF_K1460_K892.push_back( new SpinFactor("SF", SF_4Body::DtoPP1_PtoVP2_VtoP3P4, 3, 1, 2, 0) );
  std::vector<Lineshape*> LS_K1460_K892;
  LS_K1460_K892.push_back( new Lineshape("K1460", K1460M, K1460W, 1, M_34_2, LS::BW, FF::BL2) );
  LS_K1460_K892.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_34, LS::BW, FF::BL2) );
  LS_K1460_K892.push_back( new Lineshape("K1460", K1460M, K1460W, 1, M_13_2, LS::BW, FF::BL2) );
  LS_K1460_K892.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_13, LS::BW, FF::BL2) );

  std::vector<SpinFactor*> SF_K1460_f0_1370;
  SF_K1460_f0_1370.push_back( new SpinFactor("SF", SF_4Body::DtoPP1_PtoSP2_StoP3P4, 0, 1, 2, 3) );
  SF_K1460_f0_1370.push_back( new SpinFactor("SF", SF_4Body::DtoPP1_PtoSP2_StoP3P4, 3, 1, 2, 0) );
  std::vector<Lineshape*> LS_K1460_f0_1370;
  LS_K1460_f0_1370.push_back( new Lineshape("K1460", K1460M, K1460W, 0, M_12_3, LS::BW, FF::BL2) );
  LS_K1460_f0_1370.push_back( new Lineshape("f0_1370", f0_1370M, f0_1370W, 0, M_12, LS::BW, FF::BL2) );
  LS_K1460_f0_1370.push_back( new Lineshape("K1460", K1460M, K1460W, 0, M_24_3, LS::BW, FF::BL2) );
  LS_K1460_f0_1370.push_back( new Lineshape("f0_1370", f0_1370M, f0_1370W, 0, M_24, LS::BW, FF::BL2) );

  std::vector<SpinFactor*> SF_K1_1270_K892;
  SF_K1_1270_K892.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 0, 1, 2, 3) );
  SF_K1_1270_K892.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1, 1, 2, 3, 0) );
  SF_K1_1270_K892.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 3, 1, 2, 0) );
  SF_K1_1270_K892.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1, 1, 2, 0, 3) );
  std::vector<Lineshape*> LS_K1_1270_K892;
  LS_K1_1270_K892.push_back( new Lineshape("K1_1270", K1_1270M, K1_1270W, 0, M_34_2, LS::BW, FF::BL2) );
  LS_K1_1270_K892.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_34, LS::BW, FF::BL2) );
  LS_K1_1270_K892.push_back( new Lineshape("K1_1270", K1_1270M, K1_1270W, 0, M_13_2, LS::BW, FF::BL2) );
  LS_K1_1270_K892.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_13, LS::BW, FF::BL2) );

  std::vector<SpinFactor*> SF_K1_1270_rho770;
  SF_K1_1270_rho770.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 3, 2, 0, 1) );
  SF_K1_1270_rho770.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1, 0, 1, 2, 3) );
  SF_K1_1270_rho770.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 0, 2, 3, 1) );
  SF_K1_1270_rho770.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1, 1, 2, 3, 0) );
  std::vector<Lineshape*> LS_K1_1270_rho770;
  LS_K1_1270_rho770.push_back( new Lineshape("K1_1270", K1_1270M, K1_1270W, 0, M_12_3, LS::BW, FF::BL2) );
  LS_K1_1270_rho770.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) );
  LS_K1_1270_rho770.push_back( new Lineshape("K1_1270", K1_1270M, K1_1270W, 0, M_24_3, LS::BW, FF::BL2) );
  LS_K1_1270_rho770.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS::BW, FF::BL2) );

  std::vector<SpinFactor*> SF_K1_1270_K0_1430;
  SF_K1_1270_K0_1430.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoSP2_StoP3P4, 0, 1, 2, 3) );
  SF_K1_1270_K0_1430.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1, 1, 2, 3, 0) );
  SF_K1_1270_K0_1430.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoSP2_StoP3P4, 3, 1, 2, 0) );
  SF_K1_1270_K0_1430.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1, 0, 1, 2, 3) );
  std::vector<Lineshape*> LS_K1_1270_K0_1430;
  LS_K1_1270_K0_1430.push_back( new Lineshape("K(1)(1270)bar", K1_1270M, K1_1270W, 1, M_34_2 , LS::BW, FF::BL2) );
  LS_K1_1270_K0_1430.push_back( new Lineshape("K(0)*(1430)", K0_1430M, K0_1430W, 0, M_34 , LS::Lass_M3, FF::BL2, 1.5, LassVars) );
  LS_K1_1270_K0_1430.push_back( new Lineshape("K(1)(1270)bar2", K1_1270M, K1_1270W, 1, M_13_2 , LS::BW, FF::BL2) );
  LS_K1_1270_K0_1430.push_back( new Lineshape("K(0)*1430)", K0_1430M, K0_1430W, 0, M_13 , LS::Lass_M3, FF::BL2, 1.5, LassVars) );

  std::vector<SpinFactor*> SF_K1_1400_K892;
  SF_K1_1400_K892.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 0, 1, 2, 3) );
  SF_K1_1400_K892.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1, 1, 2, 3, 0) );
  SF_K1_1400_K892.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 3, 1, 2, 0) );
  SF_K1_1400_K892.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1, 1, 2, 0, 3) );
  std::vector<Lineshape*> LS_K1_1400_K892;
  LS_K1_1400_K892.push_back( new Lineshape("K1_1400", K1_1400M, K1_1400W, 0, M_34_2, LS::BW, FF::BL2) );
  LS_K1_1400_K892.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_34, LS::BW, FF::BL2) );
  LS_K1_1400_K892.push_back( new Lineshape("K1_1400", K1_1400M, K1_1400W, 0, M_13_2, LS::BW, FF::BL2) );
  LS_K1_1400_K892.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_13, LS::BW, FF::BL2) );

  std::vector<SpinFactor*> SF_K2_1430_K892;
  SF_K2_1430_K892.push_back( new SpinFactor("SF", SF_4Body::DtoTP1_TtoVP2_VtoP3P4, 0, 1, 2, 3) );
  SF_K2_1430_K892.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L2, 1, 2, 3, 0) );
  SF_K2_1430_K892.push_back( new SpinFactor("SF", SF_4Body::DtoTP1_TtoVP2_VtoP3P4, 3, 1, 2, 0) );
  SF_K2_1430_K892.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L2, 1, 2, 0, 3) );
  std::vector<Lineshape*> LS_K2_1430_K892;
  LS_K2_1430_K892.push_back( new Lineshape("K2_1430", K2_1430M, K2_1430W, 2, M_34_2, LS::BW, FF::BL2) );
  LS_K2_1430_K892.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_34, LS::BW, FF::BL2) );
  LS_K2_1430_K892.push_back( new Lineshape("K2_1430", K2_1430M, K2_1430W, 2, M_13_2, LS::BW, FF::BL2) );
  LS_K2_1430_K892.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_13, LS::BW, FF::BL2) );

  std::vector<SpinFactor*> SF_K2_1430_rho770;
  SF_K2_1430_rho770.push_back( new SpinFactor("SF", SF_4Body::DtoTP1_TtoVP2_VtoP3P4, 3, 2, 0, 1) );
  SF_K2_1430_rho770.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L2, 0, 1, 2, 3) );
  SF_K2_1430_rho770.push_back( new SpinFactor("SF", SF_4Body::DtoTP1_TtoVP2_VtoP3P4, 0, 2, 3, 1) );
  SF_K2_1430_rho770.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L2, 3, 1, 2, 0) );
  std::vector<Lineshape*> LS_K2_1430_rho770;
  LS_K2_1430_rho770.push_back( new Lineshape("K2_1430", K2_1430M, K2_1430W, 2, M_12_3, LS::BW, FF::BL2) );
  LS_K2_1430_rho770.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) );
  LS_K2_1430_rho770.push_back( new Lineshape("K2_1430", K2_1430M, K2_1430W, 2, M_24_3, LS::BW, FF::BL2) );
  LS_K2_1430_rho770.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS::BW, FF::BL2) );

  std::vector<SpinFactor*> SF_a1_f0_600;
  SF_a1_f0_600.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoSP2_StoP3P4, 2, 3, 0, 1) );
  SF_a1_f0_600.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1, 0, 1, 3, 2) );
  SF_a1_f0_600.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoSP2_StoP3P4, 2, 0, 3, 1) );
  SF_a1_f0_600.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1, 0, 1, 3, 2) );
  std::vector<Lineshape*> LS_a1_f0_600;
  LS_a1_f0_600.push_back( new Lineshape("a(1)(1260)+", a1M, a1W, 1, M_12_4, LS::BW, FF::BL2, 5.71) );
  LS_a1_f0_600.push_back( new Lineshape("f600", f600M, f600W, 0, M_12, LS::Bugg3, FF::BL2) );
  LS_a1_f0_600.push_back( new Lineshape("a(1)(1260)+", a1M, a1W, 1, M_24_1, LS::BW, FF::BL2, 5.71) );
  LS_a1_f0_600.push_back( new Lineshape("f600", f600M, f600W, 0, M_24, LS::Bugg3, FF::BL2) );

  std::vector<SpinFactor*> SF_a1_rho770;
  SF_a1_rho770.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 2, 3, 0, 1) );
  SF_a1_rho770.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1, 0, 1, 3, 2) );
  SF_a1_rho770.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 2, 0, 3, 1) );
  SF_a1_rho770.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1, 0, 1, 3, 2) );
  std::vector<Lineshape*> LS_a1_rho770;
  LS_a1_rho770.push_back( new Lineshape("a(1)(1260)+", a1M, a1W, 0, M_12_4, LS::BW, FF::BL2, 5.71) );
  LS_a1_rho770.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) );
  LS_a1_rho770.push_back( new Lineshape("a(1)(1260)+", a1M, a1W, 0, M_24_1, LS::BW, FF::BL2, 5.71) );
  LS_a1_rho770.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS::BW, FF::BL2) );

  std::vector<SpinFactor*> SF_a1_rho770_D;
  SF_a1_rho770_D.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4, 2, 3, 0, 1) );
  SF_a1_rho770_D.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1, 0, 1, 3, 2) );
  SF_a1_rho770_D.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4, 2, 0, 3, 1) );
  SF_a1_rho770_D.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1, 0, 1, 3, 2) );
  std::vector<Lineshape*> LS_a1_rho770_D;
  LS_a1_rho770_D.push_back( new Lineshape("a(1)(1260)+", a1M, a1W, 2, M_12_4, LS::BW, FF::BL2, 5.71) );
  LS_a1_rho770_D.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) );
  LS_a1_rho770_D.push_back( new Lineshape("a(1)(1260)+", a1M, a1W, 2, M_24_1, LS::BW, FF::BL2, 5.71) );
  LS_a1_rho770_D.push_back( new Lineshape("rho(770)", RhoMass, RhoWidth, 1, M_24, LS::BW, FF::BL2) );

  std::vector<SpinFactor*> SF_nonRes;
  SF_nonRes.push_back( new SpinFactor("SF", SF_4Body::ONE, 2, 3, 0, 1) );
  SF_nonRes.push_back( new SpinFactor("SF", SF_4Body::ONE, 2, 0, 3, 1) );
  std::vector<Lineshape*> LS_nonRes;
  LS_nonRes.push_back( new Lineshape("nonRes", a1M,     a1W,          0, M_12, LS::ONE, FF::BL2) );
  LS_nonRes.push_back( new Lineshape("nonRes", RhoMass, RhoWidth,     0, M_34, LS::ONE, FF::BL2) );
  LS_nonRes.push_back( new Lineshape("nonRes", a1M,     a1W,          0, M_12, LS::ONE, FF::BL2) );
  LS_nonRes.push_back( new Lineshape("nonRes", RhoMass, RhoWidth,     0, M_34, LS::ONE, FF::BL2) );

  std::vector<SpinFactor*> SF_NonResA_K892;
  SF_NonResA_K892.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4, 0, 1, 2, 3) );
  SF_NonResA_K892.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1, 1, 2, 3, 0) );
  SF_NonResA_K892.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4, 3, 1, 2, 0) );
  SF_NonResA_K892.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1, 1, 2, 0, 3) );
  std::vector<Lineshape*> LS_NonResA_K892;
  LS_NonResA_K892.push_back( new Lineshape("K1_1400", new Variable("NR1",0.0), new Variable("NR2",0.0), 2, M_34_2, LS::nonRes, FF::BL2) );
  LS_NonResA_K892.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_34, LS::BW, FF::BL2) );
  LS_NonResA_K892.push_back( new Lineshape("K1_1400", new Variable("NR3",0.0), new Variable("NR4",0.0), 2, M_13_2, LS::nonRes, FF::BL2) );
  LS_NonResA_K892.push_back( new Lineshape("K*(892)bar", K892M, K892W, 1, M_13, LS::BW, FF::BL2) );

  // the very last parameter means that we have two permutations. so the first half of the Lineshapes 
  // and the first half of the spinfactors are amplitude 1, rest is amplitude 2
  // This means that it is important for symmetrized amplitueds that the spinfactors and lineshapes are in the "right" order
  
  //RS Model
  Amplitude* amp_K892_rho770_S       = new Amplitude( "K892_rho770_S",   new Variable("K892_rho770_S_real",   1.0),     new Variable("K892_rho770_S_imag", 0.0), LS_K892_rho770_S, SF_K892_rho770_S, 2);
  Amplitude* amp_K892_rho770_P       = new Amplitude( "K892_rho770_P",   new Variable("K892_rho770_P_real",   1.0),   new Variable("K892_rho770_P_imag", 0.0), LS_K892_rho770_P, SF_K892_rho770_P , 2);
  Amplitude* amp_K892_rho770_D       = new Amplitude( "K892_rho770_D",   new Variable("K892_rho770_D_real",    1.0), new Variable("K892_rho770_D_imag",0.0), LS_K892_rho770_D, SF_K892_rho770_D, 2);
  Amplitude* amp_K1410_rho770_P      = new Amplitude( "K1410_rho770",    new Variable("K1410_rho770_P_real",   4.001),  new Variable("K1410_rho770_P_imag",-2.620), LS_K1410_rho770_P, SF_K1410_rho770_P, 2);
  Amplitude* amp_K892_f0_600         = new Amplitude( "K892_f0600",      new Variable("K892_f0600_real",      -0.770),  new Variable("K892_f0600_imag",  -1.530), LS_K892_f0_600, SF_K892_f0_600, 2);
  Amplitude* amp_rho1450_K0_1430     = new Amplitude( "rho1450_K0_1430", new Variable("rho1450_K0_1430_real", -0.110),   new Variable("rho1450_K0_1430_imag",  1.814), LS_rho1450_K0_1430  , SF_rho1450_K0_1430 , 2);
  Amplitude* amp_K1460_K892          = new Amplitude( "K1460_K892",      new Variable("K1460_K892_real",      -0.696),  new Variable("K1460_K892_imag",  0.326), LS_K1460_K892  , SF_K1460_K892 , 2);
  Amplitude* amp_K1460_f0_1370       = new Amplitude( "K1460_f0_1370",   new Variable("K1460_f0_1370_real",   -0.849),  new Variable("K1460_f0_1370_imag",  0.972), LS_K1460_f0_1370  , SF_K1460_f0_1370 , 2);
  Amplitude* amp_K1_1270_K892        = new Amplitude( "K1_1270_K892",    new Variable("K1_1270_K892_real",    0.601),   new Variable("K1_1270_K892_imag",  -0.182), LS_K1_1270_K892  , SF_K1_1270_K892 , 2);
  Amplitude* amp_K1_1270_rho770      = new Amplitude( "K1_1270_rho770",  new Variable("K1_1270_rho770_real",  -1.523),  new Variable("K1_1270_rho770_imag",  1.244), LS_K1_1270_rho770  , SF_K1_1270_rho770 , 2);
  Amplitude* amp_K1_1270_K0_1430     = new Amplitude( "K1_1270_K0_1430", new Variable("K1_1270_K0_1430_real", 0.248),  new Variable("K1_1270_K0_1430_imag",  -0.088), LS_K1_1270_K0_1430  , SF_K1_1270_K0_1430 , 2);
  Amplitude* amp_K1_1400_K892        = new Amplitude( "K1_1400_K892",    new Variable("K1_1400_K892_real",    -0.808),  new Variable("K1_1400_K892_imag",  -0.358), LS_K1_1400_K892  , SF_K1_1400_K892 , 2);
  Amplitude* amp_NonResA_K892        = new Amplitude( "NonResA_K892",    new Variable("NonResA_K892_real",    -15.322), new Variable("NonResA_K892_imag",  -12.089), LS_NonResA_K892, SF_NonResA_K892, 2);
  Amplitude* amp_K2_1430_K892        = new Amplitude( "K2_1430_K892",    new Variable("K2_1430_K892_real",    17.008),  new Variable("K2_1430_K892_imag",  -5.014), LS_K2_1430_K892  , SF_K2_1430_K892 , 2);
  Amplitude* amp_K2_1430_rho770      = new Amplitude( "K2_1430_rho770",  new Variable("K2_1430_rho770_real",  13.039),  new Variable("K2_1430_rho770_imag",  -1.935), LS_K2_1430_rho770  , SF_K2_1430_rho770 , 2);
  Amplitude* amp_a1_rho770           = new Amplitude( "a1_rho770",       new Variable("a1_rho770_real",        -0.639), new Variable("a1_rho770_imag", -6.801), LS_a1_rho770, SF_a1_rho770, 2);
  Amplitude* amp_a1_f0_600           = new Amplitude( "a1_f0_600",       new Variable("a1_f0_600_real",       -0.062),  new Variable("a1_f0_600_imag",  2.872), LS_a1_f0_600  , SF_a1_f0_600 , 2);
  Amplitude* amp_a1_rho770_D         = new Amplitude( "a1_rho770_D",     new Variable("a1_rho770_D_real",     -9.465), new Variable("a1_rho770_D_imag",  15.390), LS_a1_rho770_D, SF_a1_rho770_D, 2);
  Amplitude* amp_nonRes              = new Amplitude( "nonRes",          new Variable("nonRes_real",     -0.265),       new Variable("nonRes_imag",  -0.003), LS_nonRes, SF_nonRes, 2);


  Amplitude* amp_WS_K892_rho770_S       = new Amplitude("WS_K892_rho770_S",   new Variable("WS_K892_rho770_S_real",      1.0),     new Variable("WS_K892_rho770_S_imag",   0.0), LS_K892_rho770_S, SF_K892_rho770_S, 2);
  Amplitude* amp_WS_K892_rho770_P       = new Amplitude("WS_K892_rho770_P",   new Variable("WS_K892_rho770_P_real",      -0.109),   new Variable("WS_K892_rho770_P_imag",   1.653), LS_K892_rho770_P, SF_K892_rho770_P , 2);
  Amplitude* amp_WS_K892_rho770_D       = new Amplitude("WS_K892_rho770_D",   new Variable("WS_K892_rho770_D_real",       25.463), new Variable("WS_K892_rho770_D_imag",     2.662), LS_K892_rho770_D, SF_K892_rho770_D, 2);
  Amplitude* amp_WS_rho1450_K0_1430     = new Amplitude("WS_rho1450_K0_1430", new Variable("WS_rho1450_K0_1430_real",  2.353),   new Variable("WS_rho1450_K0_1430_imag",     -0.234), LS_rho1450_K0_1430  , SF_rho1450_K0_1430 , 2);
  Amplitude* amp_WS_K1_1270_K892        = new Amplitude("WS_K1_1270_K892",    new Variable("WS_K1_1270_K892_real",        -0.035),   new Variable("WS_K1_1270_K892_imag",    -1.405), LS_K1_1270_K892  , SF_K1_1270_K892 , 2);
  Amplitude* amp_WS_K1_1270_rho770      = new Amplitude("WS_K1_1270_rho770",  new Variable("WS_K1_1270_rho770_real",    2.42),  new Variable("WS_K1_1270_rho770_imag",       2.471), LS_K1_1270_rho770  , SF_K1_1270_rho770 , 2);
  Amplitude* amp_WS_K1_1270_K0_1430     = new Amplitude("WS_K1_1270_K0_1430", new Variable("WS_K1_1270_K0_1430_real",  -1.990),  new Variable("WS_K1_1270_K0_1430_imag",     -2.105), LS_K1_1270_K0_1430  , SF_K1_1270_K0_1430 , 2);
  Amplitude* amp_WS_K1_1400_K892        = new Amplitude("WS_K1_1400_K892",    new Variable("WS_K1_1400_K892_real",        -3.347),  new Variable("WS_K1_1400_K892_imag",     -2.237), LS_K1_1400_K892  , SF_K1_1400_K892 , 2);
  Amplitude* amp_WS_nonRes              = new Amplitude("WS_nonRes",            new Variable("WS_nonRes_real",  -0.456), new Variable("WS_nonRes_imag",              -0.041), LS_nonRes, SF_nonRes, 2);

  //DK3P_DI->amplitudes_B.push_back(amp_K892_rho770_S);
  DK3P_DI->amplitudes_B.push_back(amp_K892_rho770_P);
  DK3P_DI->amplitudes_B.push_back(amp_K892_rho770_D);
  //DK3P_DI->amplitudes_B.push_back(amp_K1410_rho770_P);
  //DK3P_DI->amplitudes_B.push_back(amp_K892_f0_600);
  //DK3P_DI->amplitudes_B.push_back(amp_rho1450_K0_1430);
  //DK3P_DI->amplitudes_B.push_back(amp_K1460_K892);
  //DK3P_DI->amplitudes_B.push_back(amp_K1460_f0_1370);
  //DK3P_DI->amplitudes_B.push_back(amp_K1_1270_K892);
  //DK3P_DI->amplitudes_B.push_back(amp_K1_1270_rho770);
  //DK3P_DI->amplitudes_B.push_back(amp_K1_1270_K0_1430);
  //DK3P_DI->amplitudes_B.push_back(amp_K1_1400_K892);
  //DK3P_DI->amplitudes_B.push_back(amp_NonResA_K892);
  //DK3P_DI->amplitudes_B.push_back(amp_K2_1430_K892);
  //DK3P_DI->amplitudes_B.push_back(amp_K2_1430_rho770);
  //DK3P_DI->amplitudes_B.push_back(amp_a1_rho770);
  //DK3P_DI->amplitudes_B.push_back(amp_a1_f0_600);
  //DK3P_DI->amplitudes_B.push_back(amp_a1_rho770_D);
  //DK3P_DI->amplitudes_B.push_back(amp_nonRes);

  //DK3P_DI->amplitudes.push_back(amp_WS_K892_rho770_S);
  //DK3P_DI->amplitudes.push_back(amp_WS_K892_rho770_P);
  //DK3P_DI->amplitudes.push_back(amp_WS_K892_rho770_D);
  //DK3P_DI->amplitudes.push_back(amp_WS_rho1450_K0_1430);
  //DK3P_DI->amplitudes.push_back(amp_WS_K1_1270_K892);
  //DK3P_DI->amplitudes.push_back(amp_WS_K1_1270_rho770);
  //DK3P_DI->amplitudes.push_back(amp_WS_K1_1270_K0_1430);
  //DK3P_DI->amplitudes.push_back(amp_WS_K1_1400_K892);
  //DK3P_DI->amplitudes.push_back(amp_WS_nonRes);

  DK3P_DI->_tau = new Variable("tau", 0.4101);
  DK3P_DI->_xmixing = new Variable("xmixing", 0.0049);
  DK3P_DI->_ymixing = new Variable("ymixing", 0.0061);
  // DK3P_DI->_xmixing = new Variable("xmixing", 0, 0.00001, -.15, .15);
  // DK3P_DI->_ymixing = new Variable("ymixing", 0, 0.00001, -.15, .15);
  DK3P_DI->_SqWStoRSrate = new Variable("SqWStoRSrate", 1.0/sqrt(300.0));  


  Variable* m12 = new Variable("m12", 0, 3);
  Variable* m34 = new Variable("m34", 0, 3); 
  Variable* cos12 = new Variable("cos12", -1, 1);
  Variable* cos34 = new Variable("m12", -1, 1);
  Variable* phi = new Variable("phi", -3.5, 3.5);
  Variable* eventNumber = new Variable("eventNumber", 0, INT_MAX);
  Variable* dtime = new Variable("dtime", 0, 10);
  Variable* sigmat = new Variable("sigmat",-3,3);
  Variable* constantOne = new Variable("constantOne", 1); 
  Variable* constantZero = new Variable("constantZero", 0);
 

  std::vector<Variable*> vars;
  vars.push_back(m12);
  vars.push_back(m34);
  vars.push_back(cos12);
  vars.push_back(cos34);
  vars.push_back(phi);
  vars.push_back(eventNumber); 
  vars.push_back(dtime); 
  vars.push_back(sigmat); 
  UnbinnedDataSet currData(vars); 

 
  DK3P_DI->_xmixing->value = strtof(argv[5], NULL);
  DK3P_DI->_ymixing->value = strtof(argv[6], NULL);

  vector<Variable*> observables;
  vector<Variable*> coefficients; 
  vector<Variable*> offsets;

  observables.push_back(m12);
  observables.push_back(m34);
  observables.push_back(cos12);
  observables.push_back(cos34);
  observables.push_back(phi);
  observables.push_back(eventNumber);
  observables.push_back(dtime);
  observables.push_back(sigmat);
  offsets.push_back(constantZero);
  offsets.push_back(constantZero);
  coefficients.push_back(constantOne); 

  TruthResolution* dat = new TruthResolution();
  PolynomialPdf* eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
  TDDP4* dp = new TDDP4("test_TD", observables, DK3P_DI, dat, eff, 0, 1);
 
  //dp->setGenDecayTimeLimit(0,3.5); // this corresponds to rougly 97% of the exponential. So this should be ok. And speeds up Generation significantly compared to [0,5] 

  TFile *file = new TFile( argv[4] , "RECREATE");
  TTree *tree = new TTree("events", "events");

  double tm12_2,tm34_2,tc12_2,tc34_2,tphi_2,tdtime_2;

  tree->Branch("m12",          &tm12_2,         "m12/D");
  tree->Branch("m34",          &tm34_2,         "m34/D");
  tree->Branch("c12",          &tc12_2,         "c12/D");
  tree->Branch("c34",          &tc34_2,         "c34/D");
  tree->Branch("phi",          &tphi_2,         "phi/D");
  tree->Branch("dtime",        &tdtime_2,       "dtime/D");

  mcbooster::FlagAcceptReject(0,0);
  int generatedEvents = 0;
  int RunNum = 0;
  int BatchSize = strtoul(argv[1], NULL,0);
  unsigned int offi = strtoul(argv[3], NULL,0);
  unsigned int genEvts =strtoul(argv[2], NULL,0);

  double wmax = 0;
  mcbooster::FlagAcceptReject FlagIt = mcbooster::FlagAcceptReject(0.1,5);

  
  while(generatedEvents < genEvts )
  {
    unsigned int keptEvts = 0;
    dp->setGenerationOffset(offi);
    auto tuple = dp->GenerateSig(BatchSize);
    auto particles = std::get<0>(tuple);
    auto variables = std::get<1>(tuple);
    auto weights = std::get<2>(tuple);
    auto flags = std::get<3>(tuple);
    // int accepted = thrust::count_if(flags.begin(), flags.end(), thrust::identity<bool>());
    ++RunNum;
    // generatedEvents += accepted;
    for (int i = 0; i < weights.size(); ++i)
    {
      if (wmax<weights[i]) wmax = weights[i];
      if (generatedEvents < genEvts && FlagIt(i,weights[i])){
        ++generatedEvents;
        ++keptEvts;
        // printf("PF %i: %s %.5g %.5g %.5g %.5g %.5g %.5g\n",i, (bool)flags[i] ? "true" : "false", weights[i], (*(variables[0]))[i], (*(variables[1]))[i], (*(variables[2]))[i], (*(variables[3]))[i], (*(variables[4]))[i]);
        tm12_2 = (*(variables[0]))[i];
        tm34_2 = (*(variables[1]))[i];
        tc12_2 = (*(variables[2]))[i];
        tc34_2 = (*(variables[3]))[i];
        tphi_2 = (*(variables[4]))[i];
        tdtime_2 = (*(variables[5]))[i];
        tree->Fill();
        // printf("Buffer %i: %.5g %.5g %.5g %.5g %.5g %.5g \n",i, (*myweights)[i],(*Buffer_m12)[i], (*Buffer_m34)[i], (*Buffer_c12)[i], (*Buffer_c34)[i], (*Buffer_phi)[i], (*Buffer_dt)[i]);
      }
    }
    fprintf(stderr,"Run # %i: x=%.6g y=%.6g Using accept-reject method leaves you with %i out of %i events.  %.4g %% of Total offset: %u\n",RunNum, DK3P_DI->_xmixing->value, DK3P_DI->_ymixing->value, keptEvts, BatchSize, generatedEvents*100.0/genEvts, offi);
    offi += BatchSize;
    delete variables[0];
    delete variables[1];
    delete variables[2];
    delete variables[3];
    delete variables[4];
    delete variables[5];

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
