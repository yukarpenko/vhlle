/******************************************************************************
*                                                                             *
*            vHLLE : a 3D viscous hydrodynamic code                           *
*            by Iurii Karpenko                                                *
*  contact:  yu.karpenko@gmail.com                                            *
*  For the detailed description please refer to:                              *
*  Comput. Phys. Commun. 185 (2014), 3016   arXiv:1312.4160                   *
*                                                                             *
*  This code can be freely used and redistributed, provided that this         *
*  copyright appear in all the copies. If you decide to make modifications    *
*  to the code, please contact the authors, especially if you plan to publish *
* the results obtained with such modified code. Any publication of results    *
* obtained using this code must include the reference to                      *
* arXiv:1312.4160 [nucl-th] or the published version of it.                   *
*                                                                             *
*******************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <sstream>
#include "fld.h"
#include "hdo.h"
#include "ic.h"
#include "ic3F.h"
#include "ickw.h"
#include "icPartUrqmd.h"
#include "icPartSMASH.h"
#include "icGlauber.h"
#include "icGubser.h"
#include "icGlissando.h"
#include "icTrento.h"
#include "eos.h"
#include "eo3.h"
#include "eo1.h"
#include "eoChiral.h"
#include "eoAZH.h"
#include "eoHadron.h"
#include "eoSmash.h"
#include "trancoeff.h"
#include "multiHydro.h"

using namespace std;

// program parameters, to be read from file
int nx, ny, nz, nevents, eosType;
int eosTypeHadron = 0;
double xmin, xmax, ymin, ymax, etamin, etamax, tau0, tauMax, tauResize, dtau;
string collSystem, outputDir, isInputFile;
double etaS, zetaS, eCrit;
int icModel,
    glauberVariable =
        1;  // icModel=1 for pure Glauber, 2 for table input (Glissando etc)
double epsilon0, Rgt, Rgz, impactPar, s0ScaleFactor;

double snn, b_min, b_max;
int projA, targA, projZ, targZ;

void setDefaultParameters() {
 tauResize = 4.0;
}

void readParameters(char *parFile) {
 char parName[255], parValue[255];
 ifstream fin(parFile);
 if (!fin.is_open()) {
  cout << "cannot open parameters file " << parFile << endl;
  exit(1);
 }
 cout << "vhlle: reading parameters from " << parFile << endl;
 while (fin.good()) {
  string line;
  getline(fin, line);
  istringstream sline(line);
  sline >> parName >> parValue;
  if (strcmp(parName, "eosType") == 0)
   eosType = atoi(parValue);
  else if (strcmp(parName, "eosTypeHadron") == 0)
   eosTypeHadron = atoi(parValue);
  else if (strcmp(parName, "nx") == 0)
   nx = atoi(parValue);
  else if (strcmp(parName, "ny") == 0)
   ny = atoi(parValue);
  else if (strcmp(parName, "nz") == 0)
   nz = atoi(parValue);
  else if (strcmp(parName, "icModel") == 0)
   icModel = atoi(parValue);
  else if (strcmp(parName, "glauberVar") == 0)
   glauberVariable = atoi(parValue);
  else if (strcmp(parName, "xmin") == 0)
   xmin = atof(parValue);
  else if (strcmp(parName, "xmax") == 0)
   xmax = atof(parValue);
  else if (strcmp(parName, "ymin") == 0)
   ymin = atof(parValue);
  else if (strcmp(parName, "ymax") == 0)
   ymax = atof(parValue);
  else if (strcmp(parName, "etamin") == 0)
   etamin = atof(parValue);
  else if (strcmp(parName, "etamax") == 0)
   etamax = atof(parValue);
  else if (strcmp(parName, "tau0") == 0)
   tau0 = atof(parValue);
  else if (strcmp(parName, "tauMax") == 0)
   tauMax = atof(parValue);
  else if (strcmp(parName, "tauGridResize") == 0)
   tauResize = atof(parValue);
  else if (strcmp(parName, "dtau") == 0)
   dtau = atof(parValue);
  else if (strcmp(parName, "e_crit") == 0)
   eCrit = atof(parValue);
  else if (strcmp(parName, "etaS") == 0)
   etaS = atof(parValue);
  else if (strcmp(parName, "zetaS") == 0)
   zetaS = atof(parValue);
  else if (strcmp(parName, "epsilon0") == 0)
   epsilon0 = atof(parValue);
  else if (strcmp(parName, "Rg") == 0)
   Rgt = atof(parValue);
  else if (strcmp(parName, "Rgz") == 0)
   Rgz = atof(parValue);
  else if (strcmp(parName, "impactPar") == 0)
   impactPar = atof(parValue);
  else if (strcmp(parName, "s0ScaleFactor") == 0)
   s0ScaleFactor = atof(parValue);
  else if (strcmp(parName, "nevents") == 0)
   nevents = atoi(parValue);
  else if (strcmp(parName, "snn") == 0)
   snn = atof(parValue);
  else if (strcmp(parName, "b_min") == 0)
   b_min = atof(parValue);
  else if (strcmp(parName, "b_max") == 0)
   b_max = atof(parValue);
  else if (strcmp(parName, "projA") == 0)
   projA = atoi(parValue);
  else if (strcmp(parName, "targA") ==0)
   targA = atoi(parValue);
  else if (strcmp(parName, "projZ") == 0)
   projZ = atoi(parValue);
  else if (strcmp(parName, "targZ") ==0)
   targZ = atoi(parValue);
  else if (parName[0] == '!')
   cout << "CCC " << sline.str() << endl;
  else
   cout << "UUU " << sline.str() << endl;
 }
}

void printParameters() {
 cout << "====== parameters ======\n";
 cout << "outputDir = " << outputDir << endl;
 cout << "eosType = " << eosType << endl;
 cout << "eosTypeHadron = " << eosTypeHadron << endl;
 cout << "nevents = " << nevents << endl;
 cout << "snn = " << snn << endl;
 cout << "projA = " << projA << endl;
 cout << "targA = " << targA << endl;
 cout << "projZ = " << projZ << endl;
 cout << "targZ = " << targZ << endl;
 cout << "nx = " << nx << endl;
 cout << "ny = " << ny << endl;
 cout << "nz = " << nz << endl;
 cout << "icModel = " << icModel << endl;
 cout << "glauberVar = " << glauberVariable << "   ! 0=epsilon,1=entropy"
      << endl;
 cout << "xmin = " << xmin << endl;
 cout << "xmax = " << xmax << endl;
 cout << "ymin = " << ymin << endl;
 cout << "ymax = " << ymax << endl;
 cout << "etamin = " << etamin << endl;
 cout << "etamax = " << etamax << endl;
 cout << "tau0 = " << tau0 << endl;
 cout << "tauMax = " << tauMax << endl;
 cout << "tauGridResize = " << tauResize << endl;
 cout << "dtau = " << dtau << endl;
 cout << "e_crit = " << eCrit << endl;
 cout << "eta/s = " << etaS << endl;
 cout << "zeta/s = " << zetaS << endl;
 cout << "epsilon0 = " << epsilon0 << endl;
 cout << "Rgt = " << Rgt << "  Rgz = " << Rgz << endl;
 cout << "impactPar = " << impactPar << endl;
 cout << "s0ScaleFactor = " << s0ScaleFactor << endl;
 cout << "======= end parameters =======\n";
}

void readCommandLine(int argc, char** argv)
{
  if(argc==1){
  cout << "no CL params - exiting.\n"; exit(1) ;
 }
 else{
  for(int iarg=1; iarg<argc-1; iarg++){
   if(strcmp(argv[iarg],"-system")==0) collSystem = argv[iarg+1];
   if(strcmp(argv[iarg],"-params")==0) readParameters(argv[iarg+1]);
   if(strcmp(argv[iarg],"-ISinput")==0) isInputFile = argv[iarg+1];
   if(strcmp(argv[iarg],"-outputDir")==0) outputDir = argv[iarg+1];
  }
  cout << "vhlle: command line parameters are:\n";
  cout << "collision system:  " << collSystem << endl;
  cout << "ini.state input:  " << isInputFile << endl;
  cout << "output directory:  " << outputDir << endl;
 }
}


Fluid* expandGrid2x(Hydro* h, EoS* eos, EoS* eosH, TransportCoeff *trcoeff) {
 Fluid* f = h->getFluid();
 if(f->getX(0) + f->getX(f->getNX()-1)>0.001
  || f->getY(0) + f->getY(f->getNY()-1)>0.001) {
   cout << "this grid expansion works only with symmetric min/max ranges\n";
   return f;
 }
 // creating a new fluid, twice the transverse size of the current one
 Fluid* fnew = new Fluid(eos, eosH, trcoeff, f->getNX(), f->getNY(), f->getNZ(),
   2.0*f->getX(0), 2.0*f->getX(f->getNX()-1), 2.0*f->getY(0), 2.0*f->getY(f->getNY()-1),
   f->getZ(0), f->getZ(f->getNZ()-1), 2.0*h->getDtau(), f->geteCrit());
 // filling the new fluid
 for(int ix=0; ix<f->getNX(); ix++)
  for(int iy=0; iy<f->getNY(); iy++)
   for(int iz=0; iz<f->getNZ(); iz++) {
    fnew->getCell(ix, iy, iz)->importVars(f->getCell(2*(ix - f->getNX()/2) + f->getNX()/2,
      2*(iy - f->getNY()/2) + f->getNY()/2, iz));
   }
 h->setFluid(fnew);  // now Hydro object operates on the new fluid
 h->setDtau(2.0*h->getDtau());
 delete f;
 return fnew;
}

// program parameters, to be read from file
// int nx, ny, nz, eosType ;
// double xmin, xmax, ymin, ymax, zmin, zmax, tau0, tauMax, dtau ;
// char outputDir[255], eosFile[255], chiBfile[255], chiSfile[255] ;
// char icInputFile [255] ;
// double T_ch, mu_b, mu_q, mu_s, gammaS, gammaFactor, exclVolume ;
// int icModel, NPART, glauberVariable=1 ;
// double epsilon0, alpha, impactPar, s0ScaleFactor ;

int main(int argc, char **argv) {
 // pointers to all the main objects
 EoS *eos;
 EoS *eosH;
 TransportCoeff *trcoeff;
 Fluid *f_p, *f_t, *f_f;
 Hydro *h_p, *h_t, *h_f;
 MultiHydro *mh;
 time_t start = 0, end;

 time(&start);

 // read parameters from file
 readCommandLine(argc, argv);
 setDefaultParameters();
 printParameters();

 // EoS for hydro evolution
 if (eosType == 0)
  eos = new EoSs("eos/Laine_nf3.dat", 3);
 else if (eosType == 1)
  eos = new EoSChiral();
 else if (eosType == 2)
  eos = new EoSAZH();
 else {
  cout << "eosType != 0,1,2\n";
  return 0;
 }

 // hadronic EoS for hypersurface creation
 if (eosTypeHadron == 0) {
   eosH = new EoSHadron((char*)"eos/eosHadronLog.dat"); //PDG hadronic EoS
 } else if (eosTypeHadron == 1) {
   eosH = new EoSSmash((char*)"eos/hadgas_eos_SMASH.dat", 101, 101, 101); //SMASH hadronic EoS
 } else {
   cout << "Unknown haronic EoS type for hypersurface creation.\n" <<
           "eosTypeHadron should be either \"0\" (PDG hadronic EoS) or " <<
           "\"1\" (SMASH hadronic EoS).\n";
   return 0;
 }

 // CFL criterion
 double dx = (xmax - xmin) / (nx - 1);
 double dy = (ymax - ymin) / (ny - 1);
 double deta = (etamax - etamin) / (nz - 1);
 cout << "Checking x-coordinate for CFL criterion:";
 if (dx > dtau && 0.1*dx < dtau) {
  cout << "OK" << endl;
 } else {
  nx = (int)0.5*(xmax - xmin) / dtau;
  if (nx % 2 == 0) nx++;
  cout << "Not OK, resizing nx to " << nx << endl;
 }

 cout << "Checking y-coordinate for CFL criterion:";
 if (dy > dtau && 0.1*dy < dtau) {
  cout << "OK" << endl;
 } else {
  ny = (int)0.5*(ymax - ymin) / dtau;
  if (ny % 2 == 0) ny++;
  cout << "Not OK, resizing ny to " << ny << endl;
 }

 cout << "Checking eta-coordinate for CFL criterion:";
 if (deta*tau0 > dtau && 0.1*deta*tau0 < dtau) {
  cout << "OK" << endl;
 } else {
  nz = (int)0.5*(etamax - etamin)*tau0 / dtau;
  if (nz % 2 == 0) nz++;
  cout << "Not OK, resizing nz to " << nz << endl;
 }




 // transport coefficients
 trcoeff = new TransportCoeff(etaS, zetaS, eos);

 f_p = new Fluid(eos, eosH, trcoeff, nx, ny, nz, xmin, xmax, ymin, ymax, etamin,
               etamax, dtau, eCrit);
 f_t = new Fluid(eos, eosH, trcoeff, nx, ny, nz, xmin, xmax, ymin, ymax, etamin,
               etamax, dtau, eCrit);
 f_f = new Fluid(eos, eosH, trcoeff, nx, ny, nz, xmin, xmax, ymin, ymax, etamin,
               etamax, dtau, eCrit);
 cout << "fluid allocation done\n";

/*
 // initilal conditions
 if (icModel == 1) {  // optical Glauber
  ICGlauber *ic = new ICGlauber(epsilon0, impactPar, tau0);
  ic->setIC(f, eos);
  delete ic;
 } else if (icModel == 2) {  // Glauber_table + parametrized rapidity dependence
  IC *ic = new IC(isInputFile.c_str(), s0ScaleFactor);
  ic->setIC(f, eos, tau0);
  delete ic;
 } else if (icModel == 3) {  // UrQMD IC
  IcPartUrqmd *ic = new IcPartUrqmd(f, isInputFile.c_str(), Rgt, Rgz, tau0);
  ic->setIC(f, eos);
  delete ic;
 } else if (icModel == 4) {  // analytical Gubser solution
  ICGubser *ic = new ICGubser();
  ic->setIC(f, eos, tau0);
  delete ic;
  }else if(icModel==5){ // IC from GLISSANDO + rapidity dependence
   IcGlissando *ic = new IcGlissando(f, isInputFile.c_str(), tau0, collSystem.c_str());
   ic->setIC(f, eos);
   delete ic;
 } else if (icModel == 6){ // SMASH IC
   IcPartSMASH *ic = new IcPartSMASH(f, isInputFile.c_str(), Rgt, Rgz, tau0);
   ic->setIC(f, eos);
   delete ic;
 } else if(icModel==7){ // IC from Trento
   IcTrento *ic = new IcTrento(f, isInputFile.c_str(), tau0, collSystem.c_str());
   ic->setIC(f, eos);
   delete ic;
 } else {
  cout << "icModel = " << icModel << " not implemented\n";
 }
 cout << "IC done\n";
*/

 IC3F *ic = new IC3F(f_p, f_t, tau0, nevents, snn, b_min, b_max, projA, targA, projZ, targZ, Rgt, tau0);
 ic->setIC(f_p, f_t, eos);
 delete ic;
 cout << "IC done\n";

 // For calculating initial anisotropy without running full hydro, uncomment following line
 //f->InitialAnisotropies(tau0) ;

 time_t tinit = 0;
 time(&tinit);
 float diff = difftime(tinit, start);
 cout << "Init time = " << diff << " [sec]" << endl;

 // hydro init
 h_p = new Hydro(f_p, eos, trcoeff, tau0, dtau);
 h_t = new Hydro(f_t, eos, trcoeff, tau0, dtau);
 h_f = new Hydro(f_f, eos, trcoeff, tau0, dtau);
 int maxstep = ceil((tauMax - tau0) / dtau);
 start = 0;
 time(&start);
 // h->setNSvalues() ; // initialize viscous terms

 mh = new MultiHydro(f_p, f_t, f_f, h_p, h_t, h_f, eos, trcoeff, dtau, eCrit);

 f_p->initOutput(outputDir.c_str(), tau0, "proj");
 f_t->initOutput(outputDir.c_str(), tau0, "targ");
 f_f->initOutput(outputDir.c_str(), tau0, "fire");
 mh->initOutput(outputDir.c_str());
 f_p->outputCorona(tau0);
 f_t->outputCorona(tau0);
 f_f->outputCorona(tau0);
 mh->getEnergyDensity();

 do {
  mh->performStep();
  f_p->outputGnuplot(h_p->getTau());
  f_t->outputGnuplot(h_t->getTau());
  f_f->outputGnuplot(h_t->getTau());
  f_p->outputSurface(h_p->getTau());
  f_t->outputSurface(h_t->getTau());
  f_f->outputSurface(h_f->getTau());
  mh->findFreezeout();
  cout << "step done, tau=" << h_p->getTau() << endl;
  if (0.1*f_p->getDz()*h_p->getTau() > h_p->getDtau()) {
   cout << "grid resize" << endl;
   f_p = expandGrid2x(h_p, eos, eosH, trcoeff);
   f_t = expandGrid2x(h_t, eos, eosH, trcoeff);
   f_f = expandGrid2x(h_f, eos, eosH, trcoeff);
   mh->setFluids(f_p, f_t, f_f, h_p, h_t, h_f);
  }
 } while(h_p->getTau() < tauMax+0.0001);

 /*bool resized = false; // flag if the grid has been resized
 do {
  // small tau: decrease timestep by makins substeps, in order
  // to avoid instabilities in eta direction (signal velocity ~1/tau)
  int nSubSteps = 1;
  while (dtau / nSubSteps >
         1.0 * h_p->getTau() * (etamax - etamin) / (nz - 1)) {
   nSubSteps *= 2;  // 0.02 in "old" coordinates
  }
  if(nSubSteps>1) {
   h_p->setDtau(h_p->getDtau() / nSubSteps);
   for (int j = 0; j < nSubSteps; j++)
    mh->performStep();
   h->setDtau(h->getDtau() * nSubSteps);
   cout << "timestep reduced by " << nSubSteps << endl;
  } else
   mh->performStep();
  f_p->outputGnuplot(h_p->getTau());
  f_t->outputGnuplot(h_t->getTau());
  f_f->outputGnuplot(h_f->getTau());
  f_p->outputSurface(h_p->getTau());
  f_t->outputSurface(h_t->getTau());
  f_f->outputSurface(h_f->getTau());
  if(h->getTau()>=tauResize and resized==false) {
   cout << "grid resize\n";
   f_p = expandGrid2x(h_p, eos, eosH, trcoeff);
   f_t = expandGrid2x(h_t, eos, eosH, trcoeff);
   f_f = expandGrid2x(h_f, eos, eosH, trcoeff);
   resized = true;
  }
 } while(h_p->getTau()<tauMax+0.0001);*/

 end = 0;
 time(&end);
 float diff2 = difftime(end, start);
 cout << "Execution time = " << diff2 << " [sec]" << endl;

 delete mh;
 delete f_p;
 delete f_t;
 delete f_f;
 delete h_p;
 delete h_t;
 delete h_f;
 delete eos;
 delete eosH;
}
