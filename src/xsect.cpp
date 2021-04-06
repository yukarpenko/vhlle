#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <TGraph.h>
#include "xsect.h"

using namespace std;

CrossSections::CrossSections(void)
{
 const double mPi = 0.1396;
 const double mN = 0.938;
 // E_kin [GeV]
 float Ekin [] =
 {0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.5,
  2.0,  2.5,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,  10.0};
  // sigma_E [fm^2]
 float sigmaE [] =
 {0.0,  0.0,  0.011, 0.04,  0.105, 0.16,  0.22,  0.25,  0.28,  0.4,
  0.49,  0.57,  0.64,  0.73,  0.81,  0.865, 0.91,  0.94,  0.97,  0.99};
  // sigma_T [fm^2]
 float sigmaT [] =
 {1.67, 1.48, 1.32, 1.23, 1.26, 1.37, 1.25, 1.10, 1.01, 0.86,
  0.75,  0.68,  0.62,  0.55,  0.49,  0.45,  0.41,  0.38,  0.35,  0.33};
 gSigmaE = new TGraph(sizeof(Ekin)/sizeof(float), Ekin, sigmaE);
 gSigmaT = new TGraph(sizeof(Ekin)/sizeof(float), Ekin, sigmaT);
 // reading tabular data for pi-p total cross section
 ifstream fin("tables/rpp2014-pimp_total.dat") ;
 if(!fin){ cout << "cannot read tables/rpp2014-pimp_total.dat\n"; exit(1) ; }
 string line ;
 istringstream instream ;
 int n; // line number
 double plab, plab_min, plab_max;
 const int dimTb = 1000;
 double s [dimTb], sigpip [dimTb];
 int i=0;
 do{
  getline(fin, line) ;
  instream.str(line) ;
  instream.seekg(0) ;
  instream.clear() ; // does not work with gcc 4.1 otherwise
  instream >> n >> plab >> plab_min >> plab_max >> sigpip[i];
  sigpip[i] = sigpip[i] * 0.1; // convers [mb] -> [fm^2]
  s[i] = mPi*mPi + mN*mN + 2.*sqrt(mPi*mPi + plab*plab)*mN;
  i++;
  if(i>dimTb-1){ cout << "CrossSections: please increase dimTb\n"; exit(1); }
 }while(!instream.fail());
 cout << "pi-p cross section: " << n << " lines read.\n";
 gSigmaPiN = new TGraph(i, s, sigpip);
}

// moments of cross section [fm^2] as a function of kinetic energy Ekin
void CrossSections::NN(double Ekin, double& sigmaT,
 double& sigmaE, double& sigmaP)
{
 const double mN = 0.94; // nucleon mass [GeV]
 if(Ekin<0.){ // catching numerical errors
  cout << "p-t friction: Ekin<0 \n";
  exit(1);
 }
 if(Ekin<=0.2){
  sigmaT = 1.239 + 0.0448/Ekin + 0.00831/(Ekin*Ekin);
  sigmaE = 0.;
 }
 else if(Ekin>0.2 && Ekin<10){
  sigmaT = gSigmaT->Eval(Ekin);
  sigmaE = gSigmaE->Eval(Ekin);
 }
 else if(Ekin>=10. && Ekin<100.){
  sigmaT = 1.365 * pow(Ekin, -0.623);
  sigmaE = 0.464 + (0.257 - 0.0125*log(Ekin)) * log(Ekin);
 }
 else if(Ekin>=100.){
  sigmaT = 0.865 * pow(Ekin, -0.525);
  sigmaE = 1.403 + (-0.15 + 0.0317 * log(Ekin)) * log(Ekin);
 }
 // ... and same formula for sigmaP for all Ekin
 sigmaP = sigmaT + (1. + 2.*mN/Ekin)*sigmaE;
}

double CrossSections::piN(double s)
{
 return gSigmaPiN->Eval(s);
}
