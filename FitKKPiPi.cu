//THIS PROGRAM FITS A SET OF DATA WITH THE RESONANCES YOU GIVE IT
//2/23/2017 
#include <fstream>

// GooFit stuff
#include "goofit/Application.h"
#include "goofit/Log.h"
#include "goofit/Variable.h" 
#include "goofit/fitting/FitManagerMinuit1.h"
#include "goofit/fitting/FitManagerMinuit2.h"
#include "goofit/PDFs/basic/PolynomialPdf.h" 
#include "goofit/PDFs/combine/AddPdf.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/physics/DP4Pdf.h"

#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/FunctionMinimum.h>

using namespace std;
using namespace GooFit;

const fptype _mD0 = 1.8645; 
const fptype piPlusMass = 0.13957018;
const fptype piMinusMass = 0.13957018; 
const fptype kPlusMass = 0.493677; 
const fptype kMinusMass = 0.493677;
// Constants used in more than one PDF component. 

int main (int argc, char** argv) {

  // Place this at the beginning of main
  GooFit::Application app{"Optional discription", argc, argv};
  
  // Command line options can be added here.
  bool minuit1;
  app.add_flag("--minuit1", minuit1, "Use Minuit 1 instead of 2");

  bool minuit2; 
  app.add_flag("--minuit2", minuit2, "Use explicit minuit2 instance");

  std::string filename;
  app.add_option("filename", filename, "The file to run")
    ->required()
    ->check(CLI::ExistingFile);
  
  try {
      app.run();
  } catch(const GooFit::ParseError &e) {
      return app.exit(e);
  }


  Variable* m12 = new Variable("m12", 0, 3);
  Variable* m34 = new Variable("m34", 0, 3); 
  Variable* cos12 = new Variable("cos12", -1, 1);
  Variable* cos34 = new Variable("cos34", -1, 1);
  Variable* phi = new Variable("phi", 0.0, 2*M_PI);
  Variable* eventNumber = new Variable("eventNumber", 0, INT_MAX);

  double Amplitudes[30]; 

  std::vector<Variable*> vars;
  vars.push_back(m12);
  vars.push_back(m34);
  vars.push_back(cos12);
  vars.push_back(cos34);
  vars.push_back(phi);
  vars.push_back(eventNumber); 
  UnbinnedDataSet currData(vars); 

  unsigned int MCevents = 0;
//Load in nTuple and give it to currData/addevent 
  fstream input(filename, std::ios_base::in);
  while(input >> *m12 >> *m34 >> *cos12 >> *cos34 >> *phi){
    //if(!*m12 || !*m34 || !*cos12 || !*cos34 || !*phi)
    //    continue;
    *eventNumber = MCevents++; 
    currData.addEvent();
  }

  printf("done reading in %i events\n", MCevents );

  DecayInfo_DP* DKKPP_DI = new DecayInfo_DP();
  DKKPP_DI->meson_radius =1.5;
  DKKPP_DI->particle_masses.push_back(_mD0);
  DKKPP_DI->particle_masses.push_back(piPlusMass);
  DKKPP_DI->particle_masses.push_back(piMinusMass);
  DKKPP_DI->particle_masses.push_back(kPlusMass);
  DKKPP_DI->particle_masses.push_back(kMinusMass);

  //Need to add K1(1270), phi, kstar0 and kstar0bar?? (this might be K1430 christoph has listed), k
  Variable* RhoMass  = new Variable("rho_mass" , 0.77526);
  Variable* RhoWidth = new Variable("rho_width", 0.1478 ); 
  Variable* K11270M = new Variable("K11270M", 1.272); 
  Variable* K11270W = new Variable("K11270W", 0.09); 
  Variable* phi1020M = new Variable("phi1020M", 1.019); 
  Variable* phi1020W = new Variable("phi1020W", 0.004); 
  Variable* K1430M   = new Variable("K1430M"   , 1.425  );//not used?
  Variable* K1430W   = new Variable("K1430W"   , 0.27   );
  Variable* Kstar1410M = new Variable("Kstar1410M" ,   1.414);
  Variable* Kstar1410W = new Variable("Kstar1410W",   0.232);
  //Not Used
  Variable* FZeroMass    = new Variable("f600M"    , 0.519  );
  Variable* FZeroWidth    = new Variable("f600W"    , 0.454  ); 
  Variable* a1M      = new Variable("a1M"      , 1.23   );
  Variable* a1W      = new Variable("a1W"      , 0.42   ); 
  Variable* K1M      = new Variable("K1M"      , 1.272  );
  Variable* K1W      = new Variable("K1W"      , 0.09   ); 
  Variable* Kstar892M   = new Variable("Kstar892M"   , 0.89581);
  Variable* Kstar892W   = new Variable("Kstar892W"   , 0.0474 ); 
  

  //Spin factors: we have two due to the bose symmetrization of the two pi+
  //K11270->Kstar0 (2?)----
  //K11270->rho (2?)----- 
  //Kstar1410-> Kstar0 (2?)-----
  //phi rho (D and S)
  //phi  (s) non resonant?
  // non resonant?

/////////////----------------------NOTE----------------
// p1-> pi+
// p2-> pi-
// p3-> K+
// p4-> K-

  std::vector<SpinFactor*> SFK1P2Kstar1430;//K1(1270)+(Kstar0 pi+)K- 
  SFK1P2Kstar1430.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 3,0,1,2));
  SFK1P2Kstar1430.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1,2,0,1,3));

  std::vector<SpinFactor*> SFK1M2Kstar1430;//K1(1270)-(Kstar0bar pi-)K+ 
  SFK1M2Kstar1430.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 2,1,0,3)); 
  SFK1M2Kstar1430.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1,3,1,0,2));

  std::vector<SpinFactor*> SFK1P2Kstar;//K1(1270)+(Kstar0 pi+)K- 
  SFK1P2Kstar.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 3,0,1,2));
  SFK1P2Kstar.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1,2,0,1,3));

  std::vector<SpinFactor*> SFK1M2Kstar;//K1(1270)-(Kstar0bar pi-)K+ 
  SFK1M2Kstar.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 2,1,0,3));
  SFK1M2Kstar.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1,3,1,0,2));

  std::vector<SpinFactor*> SFK1P2Rho;//K1(1270)+(rho K+)K- 
  SFK1P2Rho.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 3,2,0,1)); 
  SFK1P2Rho.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1,1,2,0,3));

  std::vector<SpinFactor*> SFK1M2Rho;//K1(1270)-(rho K-)K+ 
  SFK1M2Rho.push_back( new SpinFactor("SF", SF_4Body::DtoAP1_AtoVP2_VtoP3P4, 2,3,1,0)); 
  SFK1M2Rho.push_back( new SpinFactor("SF",SF_4Body::FF_123_4_L1,0,3,1,2));


  std::vector<SpinFactor*> SFKstarP2Kstar;//Kstar(1410)+(Kstar0 pi+)K- 
  SFKstarP2Kstar.push_back( new SpinFactor("SF", SF_4Body::DtoPP1_PtoVP2_VtoP3P4, 3,0,1,2));
  SFKstarP2Kstar.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1,2,0,1,3));

  std::vector<SpinFactor*> SFKstarM2Kstar;//Kstar(1410)-(Kstar0 pi-)K+ 
  SFKstarM2Kstar.push_back( new SpinFactor("SF", SF_4Body::DtoPP1_PtoVP2_VtoP3P4,2,1,0,3));  
  SFKstarM2Kstar.push_back( new SpinFactor("SF", SF_4Body::FF_123_4_L1,3,1,0,2));

  std::vector<SpinFactor*> SFKstarKstarS;
  SFKstarKstarS.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, 0,3,1,2));  

  std::vector<SpinFactor*> SFKstarKstarP;
  SFKstarKstarP.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P, 0,3,1,2));  
  SFKstarKstarP.push_back(new SpinFactor("SF",SF_4Body::FF_12_34_L1, 0,3,1,2));
  
  std::vector<SpinFactor*> SFKstarKstarD;
  SFKstarKstarD.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D, 0,3,1,2));  
  SFKstarKstarD.push_back(new SpinFactor("SF",SF_4Body::FF_12_34_L2, 0,3,1,2));
 

  std::vector<SpinFactor*> SFPhiRhoS; 
  SFPhiRhoS.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S, 2,3,0,1)); 
  //SFPhiRhoD.push_back(new SpinFactor("SF",SF_4Body::FF_12_34_L2, 2,3,0,1));
 

  std::vector<SpinFactor*> SFPhiRhoP; 
  SFPhiRhoP.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P, 2,3,0,1)); 
  SFPhiRhoP.push_back(new SpinFactor("SF",SF_4Body::FF_12_34_L1, 2,3,0,1));
 

  std::vector<SpinFactor*> SFPhiRhoD;
  SFPhiRhoD.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D, 2,3,0,1)); 
  SFPhiRhoD.push_back(new SpinFactor("SF",SF_4Body::FF_12_34_L2, 2,3,0,1));
 
  std::vector<SpinFactor*> SFPhiFZero;
  SFPhiFZero.push_back( new SpinFactor("SF", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D, 2,3,0,1)); 
  SFPhiFZero.push_back( new SpinFactor("SF",SF_4Body::FF_12_34_L2,2,3,0,1));

 // std::vector<SpinFactor*> SFPhipipi; 
 // SFPhipipi.push_back( new SpinFactor("SF", SF_4Body::DtoVP1P2_VtoP3P4, 0,1,2,3)); 

  std::vector<SpinFactor*> SFNonRes1;
  SFNonRes1.push_back( new SpinFactor("SF", SF_4Body::ONE, 0,1,2,3)); 

  std::vector<SpinFactor*> SFNonRes2; 
  SFNonRes2.push_back( new SpinFactor("SF", SF_4Body::ONE, 1,2,0,3));   

  //////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<Lineshape*> LSK1P2Kstar1430;
  LSK1P2Kstar1430.push_back( new Lineshape("K1(1270)+", K11270M, K11270W, 1, M_23_1, LS::BW, FF::BL2) );    
  LSK1P2Kstar1430.push_back( new Lineshape("Kstar(1430)", K1430M, K1430W, 0, M_23, LS::BW, FF::BL2) ); 

  std::vector<Lineshape*> LSK1M2Kstar1430; 
  LSK1M2Kstar1430.push_back( new Lineshape("K1(1270)-", K11270M, K11270W, 1, M_14_2, LS::BW, FF::BL2) ); 
  LSK1M2Kstar1430.push_back( new Lineshape("Kstar(1430)", K1430M, K1430W,0, M_14, LS::BW, FF::BL2) );  

  std::vector<Lineshape*> LSK1P2Kstar;
  LSK1P2Kstar.push_back( new Lineshape("K1(1270)+", K11270M, K11270W, 0, M_23_1, LS::BW, FF::BL2) );    
  LSK1P2Kstar.push_back( new Lineshape("Kstar(892)", Kstar892M, Kstar892W, 1, M_23, LS::BW, FF::BL2) ); 

  std::vector<Lineshape*> LSK1M2Kstar; 
  LSK1M2Kstar.push_back( new Lineshape("K1(1270)-", K11270M, K11270W, 0, M_14_2, LS::BW, FF::BL2) ); 
  LSK1M2Kstar.push_back( new Lineshape("Kstar(892)", Kstar892M, Kstar892W,1, M_14, LS::BW, FF::BL2) );  

  std::vector<Lineshape*> LSK1P2Rho;
  LSK1P2Rho.push_back( new Lineshape("K1(1270)+", K11270M, K11270W, 0, M_12_3, LS::BW, FF::BL2) );
  LSK1P2Rho.push_back( new Lineshape("Rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) );   

  std::vector<Lineshape*> LSK1M2Rho; 
  LSK1M2Rho.push_back( new Lineshape("K1(1270)-", K11270M, K11270W, 0, M_12_4, LS::BW, FF::BL2) ); 
  LSK1M2Rho.push_back( new Lineshape("Rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) ); 
 
  std::vector<Lineshape*> LSKstarP2Kstar;
  LSKstarP2Kstar.push_back( new Lineshape("Kstar1410", Kstar1410M, Kstar1410W, 1, M_23_1, LS::BW, FF::BL2) ); 
  LSKstarP2Kstar.push_back( new Lineshape("Kstar892", Kstar892M, Kstar892W, 1, M_23, LS::BW, FF::BL2) );  

  std::vector<Lineshape*> LSKstarM2Kstar; 
  LSKstarM2Kstar.push_back( new Lineshape("Kstar1410", Kstar1410M, Kstar1410W, 1, M_14_2, LS::BW, FF::BL2) ); 
  LSKstarM2Kstar.push_back( new Lineshape("Kstar(892)", Kstar892M, Kstar892W, 1, M_14, LS::BW, FF::BL2) ); 
   
  std::vector<Lineshape*> LSKstarKstarbarS; 
  LSKstarKstarbarS.push_back( new Lineshape("Kstar(892) ", Kstar892M, Kstar892W, 1, M_23, LS::BW, FF::BL2) ); 
  LSKstarKstarbarS.push_back( new Lineshape("Kstarbar(892)", Kstar892M, Kstar892W, 1, M_14, LS::BW, FF::BL2) ); 
  
  std::vector<Lineshape*> LSKstarKstarbarP; 
  LSKstarKstarbarP.push_back( new Lineshape("Kstar(892) ", Kstar892M, Kstar892W, 1, M_23, LS::BW, FF::BL2) ); 
  LSKstarKstarbarP.push_back( new Lineshape("Kstarbar(892)", Kstar892M, Kstar892W, 1, M_14, LS::BW, FF::BL2) ); 
  
  std::vector<Lineshape*> LSKstarKstarbarD; 
  LSKstarKstarbarD.push_back( new Lineshape("Kstar(892)", Kstar892M, Kstar892W, 1, M_23, LS::BW, FF::BL2) ); 
  LSKstarKstarbarD.push_back( new Lineshape("Kstarbar(892)", Kstar892M, Kstar892W, 1, M_14, LS::BW, FF::BL2) ); 
  
  std::vector<Lineshape*> LSPhiRhoS; 
  LSPhiRhoS.push_back( new Lineshape("Phi(1020) ", phi1020M, phi1020W, 1, M_34, LS::BW, FF::BL2) ); 
  LSPhiRhoS.push_back( new Lineshape("Rho(770)", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) ); 

  std::vector<Lineshape*> LSPhiRhoP; 
  LSPhiRhoP.push_back( new Lineshape("Phi(1020) ", phi1020M, phi1020W, 1, M_34, LS::BW, FF::BL2) ); 
  LSPhiRhoP.push_back( new Lineshape("Rho(770) ", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) ); 


  std::vector<Lineshape*> LSPhiRhoD; 
  LSPhiRhoD.push_back( new Lineshape("Phi ", phi1020M, phi1020W, 1, M_34, LS::BW, FF::BL2) ); 
  LSPhiRhoD.push_back( new Lineshape("Rho ", RhoMass, RhoWidth, 1, M_12, LS::BW, FF::BL2) ); 

  std::vector<Lineshape*> LSPhiFZero; 
  LSPhiFZero.push_back( new Lineshape("Phi ", phi1020M, phi1020W, 1, M_34, LS::BW, FF::BL2) ); 
  LSPhiFZero.push_back( new Lineshape("Rho ", FZeroMass, FZeroWidth, 0, M_12, LS::BW, FF::BL2) ); 



//  std::vector<Lineshape*> LSPhiRhoP;  // std::vector<Lineshape*> LSPhipipi; 
 // LSPhipipi.push_back( new Lineshape("Phi ",phi1020M,phi1020W, 1, M_34, LS::BW, FF::BL2) ); 
 // LSPhipipi.push_back( new Lineshape("pipi ", new Variable("NR5", 0.0), new Variable("NR6",0.0), 1, M_12, LS::BW, FF::BL2) ); 

  std::vector<Lineshape*> LSNonRes1; 
  LSNonRes1.push_back( new Lineshape("NonRes1 ", new Variable("NR1", 0.0), new Variable("NR2", 0.0),1, M_34_2, LS::nonRes, FF::BL2) ); 

  std::vector<Lineshape*> LSNonRes2; 
  LSNonRes2.push_back( new Lineshape("NonRes2 ", new Variable("NR3", 0.0), new Variable("NR4", 0.0),1, M_34_2, LS::nonRes, FF::BL2) ); 
 /////////////////////////////////////////////////////////////////////////////////////////////////

 // Amplitude* AMP_K1P2Kstar = new Amplitude("K1(1270)+(Kstar0 pi+)K- ", new Variable("amp_real1", -0.1, 0.001,0,0), new Variable("amp_imag1", 0.1, 0.001, 0,0), LSK1P2Kstar, SFK1P2Kstar, 1);

  Amplitude* AMP_K1P2Kstar1430 = new Amplitude("K1(1270)+(Kstar1430) pi+)K- ",new Variable("AmPK1P2Kstar1430_R", -0.1, 0.001, 0,0),new Variable("AmpK1P2Kstar1430_I", 0.1, 0.001, 0, 0), LSK1P2Kstar1430, SFK1P2Kstar1430, 1);  

  Amplitude* AMP_K1M2Kstar1430 = new Amplitude("K1(1270)-(Kstar1430 pi-)K+ ", new Variable("AmpK1M2Kstar1430_R", -0.1, 0.001,0,0), new Variable("AmpK1M2Kstar1430_I", 0.1, 0.001, 0,0), LSK1M2Kstar1430, SFK1M2Kstar1430, 1); 

  Amplitude* AMP_K1P2Kstar = new Amplitude("K1(1270)+(Kstar0) pi+)K- ",new Variable("AmPK1P2Kstar_R", -0.1, 0.001, 0,0),new Variable("AmpK1P2Kstar_I", 0.1, 0.001, 0, 0), LSK1P2Kstar, SFK1P2Kstar, 1);  

  Amplitude* AMP_K1M2Kstar = new Amplitude("K1(1270)-(Kstar0 pi-)K+ ", new Variable("AmpK1M2Kstar_R", -0.1, 0.001,0,0), new Variable("AmpK1M2Kstar_I", 0.1, 0.001, 0,0), LSK1M2Kstar, SFK1M2Kstar, 1); 

  Amplitude* AMP_K1P2Rho = new Amplitude("K1(1270)+(Rho K+) K- ", new Variable("AmpK1P2Rho_R", -0.1, 0.001,0,0), new Variable("AmpK1P2Rho_I", 0.1, 0.001, 0,0), LSK1P2Rho, SFK1P2Rho, 1);

  Amplitude* AMP_K1M2Rho = new Amplitude("K1(1270)-(Rho k-)K+ ", new Variable("AmpK1M2Rho_R", -0.1, 0.001,0,0), new Variable("AmpK1M2Rho_I", 0.1, 0.001, 0,0), LSK1M2Rho, SFK1M2Rho, 1);

  Amplitude* AMP_KstarP2Kstar = new Amplitude("Kstar(1410)+(Kstar pi+)K- ", new Variable("AmpKstarP2Kstar_R", -0.1, 0.001,0,0), new Variable("AmpKstarP2Kstar_I", 0.1, 0.001, 0,0), LSKstarP2Kstar, SFKstarP2Kstar, 1);

  Amplitude* AMP_KstarM2Kstar = new Amplitude("Kstar(1410)-(Kstarbar pi-) K+ ", new Variable("AMpKstarM2Kstar_R", -0.1, 0.001,0,0), new Variable("AmpKstarM2Kstar_I", 0.1, 0.001, 0,0), LSKstarM2Kstar, SFKstarM2Kstar, 1);

  Amplitude* AMP_KstarKstarbarS = new Amplitude("KstarKstarS ", new Variable("AmpKstarKstarbarS_R", -0.1, 0.001,0,0), new Variable("AmpKstarKstarbarS_I",  -0.1,0.001, 0,0), LSKstarKstarbarS, SFKstarKstarS, 1);


  Amplitude* AMP_KstarKstarbarP = new Amplitude("KstarKstarP ", new Variable("AmpKstarKstarbarP_R", -0.1, 0.001,0,0), new Variable("AmpKstarKstarbarP_I",  -0.1,0.001, 0,0), LSKstarKstarbarP, SFKstarKstarP, 1);


  Amplitude* AMP_KstarKstarbarD = new Amplitude("KstarKstarD ", new Variable("AmpKstarKstarbarD_R", -0.1, 0.001,0,0), new Variable("AmpKstarKstarbarD_I",  -0.1,0.001, 0,0), LSKstarKstarbarD, SFKstarKstarD, 1);

  Amplitude* AMP_PhiRhoSFix = new Amplitude("PhiRhoS", new Variable("AmpPhiRhoS_R", 1), new Variable("AmpPhiRhoS_I", 0), LSPhiRhoS, SFPhiRhoS, 1);

  Amplitude* AMP_PhiRhoP = new Amplitude("PhiRhoP", new Variable("AmpPhiRhoP_R", -0.1, 0.001,0,0), new Variable("AmpPhiRhoP_I", 0.1, 0.001, 0,0), LSPhiRhoP, SFPhiRhoP, 1);


  Amplitude* AMP_PhiRhoD = new Amplitude("PhiRhoD", new Variable("AmpPhiRhoD_R", -0.1, 0.001,0,0), new Variable("AmpPhiRhoD_I", 0.1, 0.001, 0,0), LSPhiRhoD, SFPhiRhoD, 1);

 
  Amplitude* AMP_PhiFZero = new Amplitude("PhiFZero", new Variable("AmpPhiFZero_R", -0.1, 0.001,0,0), new Variable("AmpPhiFZero_I", 0.1, 0.001, 0,0), LSPhiFZero, SFPhiFZero, 1);

// Amplitude* AMP_PhiPiPi = new Amplitude("PhiPiPi", new Variable("amp_real10", -0.1, 0.001,0,0), new Variable("amp_imag10", 0.1, 0.001, 0,0), LSPhipipi, SFPhipipi, 1);

  Amplitude* AMP_NonRes1 = new Amplitude("NonRes1", new Variable("amp_real11", -0.1, 0.001,0,0), new Variable("amp_imag11", 0.1, 0.001, 0,0), LSNonRes1, SFNonRes1, 1);

  Amplitude* AMP_NonRes2 = new Amplitude("NonRes2", new Variable("amp_real12", -0.1, 0.001,0,0), new Variable("amp_imag12", 0.1, 0.001, 0,0), LSNonRes1, SFNonRes2, 1);


  DKKPP_DI->amplitudes.push_back(AMP_K1P2Kstar1430);
  DKKPP_DI->amplitudes.push_back(AMP_K1M2Kstar1430);
  DKKPP_DI->amplitudes.push_back(AMP_K1P2Kstar);
  DKKPP_DI->amplitudes.push_back(AMP_K1M2Kstar);
  DKKPP_DI->amplitudes.push_back(AMP_K1P2Rho);
  DKKPP_DI->amplitudes.push_back(AMP_K1M2Rho);
  DKKPP_DI->amplitudes.push_back(AMP_KstarP2Kstar);
  DKKPP_DI->amplitudes.push_back(AMP_KstarM2Kstar);
  DKKPP_DI->amplitudes.push_back(AMP_KstarKstarbarS); 
  DKKPP_DI->amplitudes.push_back(AMP_KstarKstarbarP);
  DKKPP_DI->amplitudes.push_back(AMP_KstarKstarbarD);
  DKKPP_DI->amplitudes.push_back(AMP_PhiRhoSFix);
  DKKPP_DI->amplitudes.push_back(AMP_PhiRhoP);
  DKKPP_DI->amplitudes.push_back(AMP_PhiRhoD);
  DKKPP_DI->amplitudes.push_back(AMP_PhiFZero); 
  //DKKPP_DI->amplitudes.push_back(AMP_PhiPiPi);
  //DKKPP_DI->amplitudes.push_back(AMP_NonRes1);
  //DKKPP_DI->amplitudes.push_back(AMP_NonRes2);




  for (auto res = LSK1P2Kstar1430.begin(); res != LSK1P2Kstar1430.end(); ++res) {
    (*res)->setParameterConstantness(true); 
  }

  for (auto res = LSK1M2Kstar1430.begin(); res != LSK1M2Kstar1430.end(); ++res) {
    (*res)->setParameterConstantness(true); 
  }

  for (auto res = LSK1P2Kstar.begin(); res != LSK1P2Kstar.end(); ++res) {
    (*res)->setParameterConstantness(true); 
  }

  for (auto res = LSK1M2Kstar.begin(); res != LSK1M2Kstar.end(); ++res) {
    (*res)->setParameterConstantness(true); 
  }

  for (auto res = LSK1P2Rho.begin(); res != LSK1P2Rho.end(); ++res) {
    (*res)->setParameterConstantness(true); 
  }

  for (auto res = LSK1M2Rho.begin(); res != LSK1M2Rho.end(); ++res) {
    (*res)->setParameterConstantness(true); 
  }
  for (auto res = LSKstarP2Kstar.begin(); res != LSKstarP2Kstar.end(); ++res) {
    (*res)->setParameterConstantness(true); 
  }

  for (auto res = LSKstarM2Kstar.begin(); res != LSKstarM2Kstar.end(); ++res) {
    (*res)->setParameterConstantness(true); 
  }

  for (auto res = LSKstarKstarbarS.begin(); res != LSKstarKstarbarS.end(); ++res) {
    (*res)->setParameterConstantness(true); 
  }


  for (auto res = LSKstarKstarbarP.begin(); res != LSKstarKstarbarP.end(); ++res) {
    (*res)->setParameterConstantness(true); 
  }


  for (auto res = LSKstarKstarbarD.begin(); res != LSKstarKstarbarD.end(); ++res) {
    (*res)->setParameterConstantness(true); 
  }

  for (auto res = LSPhiRhoS.begin(); res != LSPhiRhoS.end(); ++res) {
    (*res)->setParameterConstantness(true); 
  }


  for (auto res = LSPhiRhoP.begin(); res != LSPhiRhoP.end(); ++res) {
    (*res)->setParameterConstantness(true); 
  }


  for (auto res = LSPhiRhoD.begin(); res != LSPhiRhoD.end(); ++res) {
    (*res)->setParameterConstantness(true); 
  }

  for (auto res = LSPhiFZero.begin(); res != LSPhiFZero.end(); ++res) {
    (*res)->setParameterConstantness(true); 
  }

 // for (auto res = LSPhipipi.begin(); res != LSPhipipi.end(); ++res) {
 //   (*res)->setParameterConstantness(true); 
 // }


 // for (auto res = LSNonRes1.begin(); res != LSNonRes1.end(); ++res) {
//    (*res)->setParameterConstantness(true); 
//  }



 // for (auto res = LSNonRes2.begin(); res != LSNonRes2.end(); ++res) {
  //  (*res)->setParameterConstantness(true); 
 // }



  Variable* constantOne = new Variable("constantOne", 1); 
  Variable* constantZero = new Variable("constantZero", 0);

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

  PolynomialPdf* eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
  DPPdf* dp = new DPPdf("test", observables, DKKPP_DI, eff,1e6);

  Variable* constant = new Variable("constant", 0.1); 
  Variable* constant2 = new Variable("constant2", 1.0); 
  vars.clear();
  vars.push_back(constant);
  PolynomialPdf backgr("backgr", m12, vars); 
  AddPdf* signal = new AddPdf("signal",constant2,dp, &backgr);

  signal->setData(&currData);
  dp->setDataSize(currData.getNumEvents(), 6); 

  if(minuit1) {
    GooFit::FitManagerMinuit1 datapdf(signal);
    datapdf.useHesseBefore(false);
    datapdf.fit();
    return 0; 
  } else {

    if(minuit2) {
      GooFit::Params upar{*signal};
      GooFit::FCN fcn{upar};
      Minuit2::MnPrint::SetLevel(3);
      Minuit2::MnMigrad migrad{fcn, upar};
      Minuit2::FunctionMinimum min = migrad(10000);
      cout << min << endl;
      return 0;
    } else {
      GooFit::FitManagerMinuit2 datapdf(signal);
      datapdf.setMaxCalls(10000);
      datapdf.fit();
      return datapdf; 
    }
  }
}
