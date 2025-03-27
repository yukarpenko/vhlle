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

#include <cstring>
#include <ctime>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include "fld.h"
#include "hdo.h"
#include "ic.h"
#include "ickw.h"
#include "icPartUrqmd.h"
#include "icPartSMASH.h"
#include "icGlauber.h"
#include "icGubser.h"
#include "icGlissando.h"
#include "icSuperMC.h"
#include "icTrento.h"
#include "icTrento3d.h"
#include "eos.h"
#include "eo3.h"
#include "eo1.h"
#include "eoChiral.h"
#include "eoCMF.h"
#include "eoCMFe.h"
#include "eoAZH.h"
#include "eoHadron.h"
#include "eoSmash.h"
#include "trancoeff.h"
#include "vtk.h"
#include "colour.h"


using namespace std;

void checkGridDimension(int n, char _x) {
  if (n < 5) {
   std::cerr << red << "FATAL: The number of cells in the hydrodynamic grid is too small! \n";
   std::cerr << "n" << _x << " < 5 \n" << reset;
   exit(1);
  }
}

void checkGridBorders(double min, double max, std::string _x) {
  if (min >= max) {
   std::cerr << red << "FATAL: The hydrodynamic grid is not set up correctly! \n";
   std::cerr << _x << "min >= " << _x << "max \n" << reset;
   exit(1);
  }
}

int nx {100}, ny {100}, nz {100}, eosType {1}, etaSparam {0}, zetaSparam {0}, eosTypeHadron {0};
// only FO hypersurface output: {0,1};  freezeout output extended by e,nb: {0,1}
bool vtk_cartesian {false}, vtk {false}, freezeoutOnly {false}, freezeoutExtend {false}, vorticityOn {false}; 
double xmin {-5.0}, xmax {5.0}, ymin {-5.0}, ymax {5.0}, etamin {-5.0}, 
  etamax {5.0}, tau0 {1.0}, tauMax {20.0}, tauResize {4.0}, dtau {0.05},
  etaS {0.08}, zetaS {0.0}, eCrit {0.5}, etaSEpsilonMin {5.}, al {0.}, ah {0.}, aRho {0.}, T0 {0.15}, 
  etaSMin {0.08}, etaSShiftMuB {0.}, etaSScaleMuB {0.}, zetaSPeakEpsilon {5.}, 
  zetaSScaleBeta {0.103}, zetaSSigmaMinus {0.1}, zetaSSigmaPlus {0.1}, epsilon0, Rgt {1.0},
  Rgz {1.0}, impactPar, s0ScaleFactor;
string collSystem, outputDir {"data"}, isInputFile, vtk_values {""};
int icModel {1},glauberVariable  {1};  // icModel=1 for pure Glauber, 2 for table input (Glissando etc) 
int smoothingType {0}; // 0 for kernel contracted in eta, 1 for invariant kernel 

void readParameters(char *parFile) {
    char parName[255], parValue[255];
    ifstream fin(parFile);
    if (!fin.is_open()) {
        cout << "cannot open parameters file " << parFile << endl;
        exit(1);
    }
    cout << "vhlle: reading parameters from " << parFile << endl;

    map<string, function<void(const string&)>> handlers = {
        {"eosType", [](const string& value) { eosType = atoi(value.c_str()); }},
        {"eosTypeHadron", [](const string& value) { eosTypeHadron = atoi(value.c_str()); }},
        {"nx", [](const string& value) { nx = atoi(value.c_str()); }},
        {"ny", [](const string& value) { ny = atoi(value.c_str()); }},
        {"nz", [](const string& value) { nz = atoi(value.c_str()); }},
        {"icModel", [](const string& value) { icModel = atoi(value.c_str()); }},
        {"glauberVar", [](const string& value) { glauberVariable = atoi(value.c_str()); }},
        {"xmin", [](const string& value) { xmin = atof(value.c_str()); }},
        {"xmax", [](const string& value) { xmax = atof(value.c_str()); }},
        {"ymin", [](const string& value) { ymin = atof(value.c_str()); }},
        {"ymax", [](const string& value) { ymax = atof(value.c_str()); }},
        {"etamin", [](const string& value) { etamin = atof(value.c_str()); }},
        {"etamax", [](const string& value) { etamax = atof(value.c_str()); }},
        {"tau0", [](const string& value) { tau0 = atof(value.c_str()); }},
        {"tauMax", [](const string& value) { tauMax = atof(value.c_str()); }},
        {"tauGridResize", [](const string& value) { tauResize = atof(value.c_str()); }},
        {"dtau", [](const string& value) { dtau = atof(value.c_str()); }},
        {"e_crit", [](const string& value) { eCrit = atof(value.c_str()); }},
        {"etaS", [](const string& value) { etaS = atof(value.c_str()); }},
        {"zetaS", [](const string& value) { zetaS = atof(value.c_str()); }},
        {"zetaSparam", [](const string& value) { zetaSparam = atoi(value.c_str()); }},
        {"zetaSScaleBeta", [](const string& value) { zetaSScaleBeta = atof(value.c_str()); }},
        {"zetaSPeakEpsilon", [](const string& value) { zetaSPeakEpsilon = atof(value.c_str()); }},
        {"zetaSSigmaMinus", [](const string& value) { zetaSSigmaMinus = atof(value.c_str()); }},
        {"zetaSSigmaPlus", [](const string& value) { zetaSSigmaPlus = atof(value.c_str()); }},
        {"epsilon0", [](const string& value) { epsilon0 = atof(value.c_str()); }},
        {"Rg", [](const string& value) { Rgt = atof(value.c_str()); }},
        {"Rgz", [](const string& value) { Rgz = atof(value.c_str()); }},
        {"impactPar", [](const string& value) { impactPar = atof(value.c_str()); }},
        {"s0ScaleFactor", [](const string& value) { s0ScaleFactor = atof(value.c_str()); }},
        {"VTK_output", [](const string& value) { vtk = atoi(value.c_str()); }},
        {"VTK_output_values", [](const string& value) { vtk_values = value; }},
        {"VTK_cartesian", [](const string& value) { vtk_cartesian= atoi(value.c_str()); }},
        {"etaSparam", [](const string& value) { etaSparam = atoi(value.c_str()); }},
        {"aRho", [](const string& value) { aRho = atof(value.c_str()); }},
        {"ah", [](const string& value) { ah = atof(value.c_str()); }},
        {"al", [](const string& value) { al = atof(value.c_str()); }},
        {"T0", [](const string& value) { T0 = atof(value.c_str()); }},
        {"etaSEpsilonMin", [](const string& value) { etaSEpsilonMin = atof(value.c_str()); }},
        {"etaSMin", [](const string& value) { etaSMin = atof(value.c_str()); }},
        {"etaSShiftMuB", [](const string& value) { etaSShiftMuB = atof(value.c_str()); }},
        {"etaSScaleMuB", [](const string& value) { etaSScaleMuB = atof(value.c_str()); }},
        {"freezeoutOnly", [](const string& value) { freezeoutOnly = atoi(value.c_str()); }},
        {"freezeoutExtend", [](const string& value) { freezeoutExtend = atoi(value.c_str()); }},
        {"vorticity", [](const string& value) { vorticityOn = atoi(value.c_str()); }},
        {"smoothingType", [](const string& value) { smoothingType = atoi(value.c_str()); }},
    };

    while (fin.good()) {
        string line;
        getline(fin, line);
        istringstream sline(line);
        sline >> parName >> parValue;
        auto handler = handlers.find(parName);
        if (handler != handlers.end()) {
            handler->second(parValue);
        } else if (parName[0] == '!') {
            cout << "CCC " << sline.str() << endl;
        } else {
            cout << "UUU " << sline.str() << endl;
        }
    }
    checkGridBorders(xmin, xmax, "x");
    checkGridBorders(ymin, ymax, "y");
    checkGridBorders(etamin, etamax, "eta");
}

void printParameters() {
 cout << "====== parameters ======\n";
 cout << "outputDir = " << outputDir << endl;
 cout << "freezeoutOnly = " << freezeoutOnly << endl;
 cout << "freezeoutExtend = " << freezeoutExtend << endl;
 cout << "vorticity = " << vorticityOn << endl;
 cout << "eosType = " << eosType << endl;
 cout << "eosTypeHadron = " << eosTypeHadron << endl;
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
 cout << "zeta/s param : " << zetaSparam << endl;
 cout << "etaSparam = " << etaSparam << endl;
 if (etaSparam == 0){
    cout << "eta/s = " << etaS << endl;
 }
 else if (etaSparam == 1){
    cout << "al = " << al << endl;
    cout << "ah = " << ah << endl;
    cout << "etaSMin = " << etaSMin << endl;
    cout << "T0 = " << T0 << endl;
 }
 else if (etaSparam == 2){
    cout << "al = " << al << endl;
    cout << "ah = " << ah << endl;
    cout << "aRho = " << aRho << endl;
    cout << "etaSMin = " << etaSMin << endl;
    cout << "etaSEpsilonMin = " << etaSEpsilonMin << endl;
 }
 else if (etaSparam == 3){
    cout << "al = " << al << endl;
    cout << "ah = " << ah << endl;
    cout << "T0 = " << T0 << endl;
    cout << "etaSMin = " << etaSMin << endl;
    cout << "etaSShiftMuB = " << etaSShiftMuB << endl;
    cout << "etaSScaleMuB = " << etaSScaleMuB << endl;
 }
 cout << "zeta/s = " << zetaS << endl;
 if (zetaSparam == 4){
    cout << "zetaSPeakEpsilon = " << zetaSPeakEpsilon << endl;
    cout << "zetaSScaleBeta = " << zetaSScaleBeta << endl;
    cout << "zetaSSigmaMinus = " << zetaSSigmaMinus << endl;
    cout << "zetaSSigmaPlus = " << zetaSSigmaPlus << endl;
 }
 cout << "Rgt = " << Rgt << endl;
 cout << "Rgz = " << Rgz << endl;
 cout << "epsilon0 = " << epsilon0 << endl;
 cout << "smoothingType = " << smoothingType << endl;
 cout << "impactPar = " << impactPar << endl;
 cout << "s0ScaleFactor = " << s0ScaleFactor << endl;
 cout << "VTK output = " << vtk << endl;
 cout << "VTK output values = " << vtk_values << endl;
 cout << "VTK cartesian = " << vtk_cartesian << endl;
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
 if (vorticityOn) fnew->enableVorticity();
 // filling the new fluid
 for(int ix=0; ix<f->getNX(); ix++)
  for(int iy=0; iy<f->getNY(); iy++)
   for(int iz=0; iz<f->getNZ(); iz++) {
    fnew->getCell(ix, iy, iz)->importVars(f->getCell(2*(ix - f->getNX()/2) + f->getNX()/2,
      2*(iy - f->getNY()/2) + f->getNY()/2, iz));
    // enable vorticity in all cells if it was enabled in the original fluid



    // Can be deleted
    if (vorticityOn) fnew->getCell(ix, iy, iz)->enableVorticity();




   }
 h->setFluid(fnew);  // now Hydro object operates on the new fluid
 h->setDtau(2.0*h->getDtau());
 delete f;
 return fnew;
}

int main(int argc, char **argv) {
 // pointers to all the main objects
 EoS *eos;
 EoS *eosH;
 TransportCoeff *trcoeff;
 Fluid *f;
 Hydro *h;
 time_t start = 0, end;

 time(&start);

 // read parameters from file
 readCommandLine(argc, argv);
 printParameters();

 // EoS for hydro evolution
 if (eosType == 0)
  eos = new EoSs("eos/Laine_nf3.dat", 3);
 else if (eosType == 1)
  eos = new EoSChiral();
 else if (eosType == 2)
  eos = new EoSAZH();
 else if (eosType == 3)
     eos = new EoSCMF();
 else if (eosType == 4)
     eos = new EoSCMFe();
 else {
  cout << "eosType != 0,1,2,3,4\n";
  return 0;
 }

 // hadronic EoS for hypersurface creation
 if (eosTypeHadron == 0) {
   eosH = new EoSHadron((char*)"eos/eosHadronLog.dat"); //PDG hadronic EoS
 } else if (eosTypeHadron == 1) {
   eosH = new EoSSmash((char*)"eos/hadgas_eos_SMASH.dat", 101, 51, 51); //SMASH hadronic EoS
 } else {
   cout << "Unknown hadronic EoS type for hypersurface creation.\n" <<
           "eosTypeHadron should be either \"0\" (PDG hadronic EoS) or " <<
           "\"1\" (SMASH hadronic EoS).\n";
   return 0;
 }


 // transport coefficients
 trcoeff = new TransportCoeff(etaS, zetaS, ah, al, aRho, T0, etaSMin, etaSEpsilonMin, etaSShiftMuB,
  etaSScaleMuB, zetaSPeakEpsilon, zetaSScaleBeta, zetaSSigmaMinus, zetaSSigmaPlus,  eos, etaSparam, zetaSparam);

 f = new Fluid(eos, eosH, trcoeff, nx, ny, nz, xmin, xmax, ymin, ymax, etamin,
               etamax, dtau, eCrit);
 cout << "fluid allocation done\n";
 // initial conditions
 if (icModel == 1) { // optical Glauber
  ICGlauber *ic = new ICGlauber(epsilon0, impactPar, tau0);
  ic->setIC(f, eos);
  delete ic;
 } else if (icModel == 2) { // Glauber_table + parametrized rapidity dependence
  IC *ic = new IC(isInputFile.c_str(), s0ScaleFactor);
  ic->setIC(f, eos, tau0);
  delete ic;
 } else if (icModel == 3) { // UrQMD IC
  IcPartUrqmd *ic = new IcPartUrqmd(f, isInputFile.c_str(), Rgt, Rgz, tau0);
  ic->setIC(f, eos);
  delete ic;
 } else if (icModel == 4) { // analytical Gubser solution
  ICGubser *ic = new ICGubser();
  ic->setIC(f, eos, tau0);
  delete ic;
 } else if(icModel==5) { // IC from GLISSANDO + rapidity dependence
   IcGlissando *ic = new IcGlissando(f, isInputFile.c_str(), tau0, collSystem.c_str());
   ic->setIC(f, eos);
   delete ic;
 } else if (icModel == 6) { // SMASH IC
   IcPartSMASH *ic;
   ic = new IcPartSMASH(f, isInputFile.c_str(), Rgt, Rgz, smoothingType);
   tau0 = ic->getTau0();
   ic->setIC(f, eos);
   delete ic;
 } else if(icModel == 7) { // IC from Trento
   IcTrento *ic = new IcTrento(f, isInputFile.c_str(), tau0, collSystem.c_str());
   ic->setIC(f, eos);
   delete ic;  
 } else if(icModel == 8) { // IC from Trento
   IcTrento3d *ic = new IcTrento3d(f, isInputFile.c_str(), tau0, collSystem.c_str());
   ic->setIC(f, eos);
   delete ic;
 } else if(icModel == 9) { // IC from SuperMC
   //For this setup collSystem MUST be a path to a file containing the parameters for superMC
   ICSuperMC *ic = new ICSuperMC(isInputFile, tau0, collSystem);
   ic->setIC(f, eos);
   delete ic;
 } else {
   cout << "icModel = " << icModel << " not implemented\n";
 }
 cout << "IC done\n";

 // For calculating initial anisotropy without running full hydro, uncomment following line
 //f->InitialAnisotropies(tau0) ;

 time_t tinit = 0;
 time(&tinit);
 float diff = difftime(tinit, start);
 cout << "Init time = " << diff << " [sec]" << endl;

 // hydro init
 h = new Hydro(f, eos, trcoeff, tau0, dtau);

 // Enable vorticity if key is set in the config file
 if (vorticityOn) {
  h -> enableVorticity();
 }

 start = 0;
 time(&start);
 // h->setNSvalues() ; // initialize viscous terms

 f->initOutput(outputDir.c_str(), tau0, freezeoutOnly);
 f->outputCorona(tau0, freezeoutExtend);
 if (vorticityOn) f->printDbetaHeader();

 bool resized = false; // flag if the grid has been resized
 
 std::string dir=outputDir.c_str();
 VtkOutput vtk_out=VtkOutput(dir,eos,xmin,ymin,etamin, vtk_cartesian);

 int nelements = 0;
 do {
  // small tau: decrease timestep by making substeps, in order
  // to avoid instabilities in eta direction (signal velocity ~1/tau)
  if(vtk>0) {
    vtk_out.write(*h,vtk_values);
  }
  int nSubSteps = 1;
  while (dtau / nSubSteps >
         1.0 * h->getTau() * (etamax - etamin) / (nz - 1)) {
   nSubSteps *= 2;  // 0.02 in "old" coordinates
  }
  if(nSubSteps>1) {
   h->setDtau(h->getDtau() / nSubSteps);
   for (int j = 0; j < nSubSteps; j++){
    h->performStep();
   }
   h->setDtau(h->getDtau() * nSubSteps);
   cout << "timestep reduced by " << nSubSteps << endl;
  } else {
   h->performStep();
  }
  nelements = f->outputSurface(h->getTau(), freezeoutExtend);
  if (!freezeoutOnly)
   f->outputGnuplot(h->getTau());
  if(h->getTau()>=tauResize and resized==false) {
   cout << "grid resize\n";
   f = expandGrid2x(h, eos, eosH, trcoeff);
   resized = true;
  }
 } while((h->getTau()<tauMax+0.0001) and (nelements>0));

 end = 0;
 time(&end);
 float diff2 = difftime(end, start);
 cout << "Execution time = " << diff2 << " [sec]" << endl;

 f->renameOutput(outputDir.c_str());

 delete f;
 delete h;
 delete eos;
 delete eosH;
}
