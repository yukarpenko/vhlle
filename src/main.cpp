/******************************************************************************
*                                                                             *
*            vHLLE : a 3D viscous hydrodynamic code                           *
*            version 1.1,            October 2014                            *
*            by Iurii Karpenko                                                *
*  contact:  yu.karpenko@gmail.com                                            *
*  For the detailed description please refer to:                              *
*  http://arxiv.org/abs/1312.4160                                             *
*                                                                             *
*  This code can be freely used and redistributed, provided that this         *
*  copyright appear in all the copies. If you decide to make modifications    *
*  to the code, please contact the authors, especially if you plan to publish *
* the results obtained with such modified code. Any publication of results    *
* obtained using this code must include the reference to                      *
* arXiv:1312.4160 [nucl-th] or the published version of it, when available.   *
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
#include "ickw.h"
#include "icPartUrqmd.h"
#include "icGlauber.h"
#include "icGubser.h"
#include "icGlissando.h"
#include "icEpos.h"
#include "eos.h"
#include "eo3.h"
#include "eo1.h"
#include "eoChiral.h"
#include "eoAZH.h"
#include "eoHadron.h"
#include "trancoeff.h"
//--- parton cascade modules
#include "const.h"
#include "params.h"
#include "cascade.h"

//using namespace std;  // already declared in JT/const.h

// program parameters, to be read from file
int nx, ny, nz, eosType;
double xmin, xmax, ymin, ymax, etamin, etamax, tau0, tauMax, dtau;
char outputDir[255];
char icInputFile[255];
double etaS, zetaS, eCrit;
int icModel, nSmear = 1,
    glauberVariable =
        1;  // icModel=1 for pure Glauber, 2 for table input (Glissando etc)
double epsilon0, Rgt, Rgz, impactPar, s0ScaleFactor;
// ##### jet-related parameters
int eventNo; // enumerates the initial/input configuration
int jetOversampleFactor;
double jetMinPt;

void readParameters(char *parFile) {
 char parName[255], parValue[255];
 ifstream fin(parFile);
 if (!fin.is_open()) {
  cout << "cannot open parameters file " << parFile << endl;
  exit(1);
 }
 while (fin.good()) {
  string line;
  getline(fin, line);
  istringstream sline(line);
  sline >> parName >> parValue;
  if (strcmp(parName, "outputDir") == 0)
   strcpy(outputDir, parValue);
  else if (strcmp(parName, "eosType") == 0)
   eosType = atoi(parValue);
  else if (strcmp(parName, "icInputFile") == 0)
   strcpy(icInputFile, parValue);
  else if (strcmp(parName, "nx") == 0)
   nx = atoi(parValue);
  else if (strcmp(parName, "ny") == 0)
   ny = atoi(parValue);
  else if (strcmp(parName, "nz") == 0)
   nz = atoi(parValue);
  else if (strcmp(parName, "icModel") == 0)
   icModel = atoi(parValue);
  else if (strcmp(parName, "jetELsmear") == 0)
   nSmear = atoi(parValue);
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
  // ### JET-related parameters
  else if (strcmp(parName, "jetOversampleFactor") == 0)
   jetOversampleFactor = atoi(parValue);
  else if (strcmp(parName, "jetMinPt") == 0)
   jetMinPt = atof(parValue);
  else if (parName[0] == '!')
   cout << "CCC " << sline.str() << endl;
  else
   cout << "UUU " << sline.str() << endl;
 }
 fin.close();
}

void printParameters() {
 cout << "====== parameters ======\n";
 cout << "event #" << eventNo << endl;
 cout << "outputDir = " << outputDir << endl;
 cout << "eosType = " << eosType << endl;
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
 cout << "dtau = " << dtau << endl;
 cout << "e_crit = " << eCrit << endl;
 cout << "eta/s = " << etaS << endl;
 cout << "zeta/s = " << zetaS << endl;
 cout << "epsilon0 = " << epsilon0 << endl;
 cout << "Rgt = " << Rgt << "  Rgz = " << Rgz << endl;
 cout << "impactPar = " << impactPar << endl;
 cout << "s0ScaleFactor = " << s0ScaleFactor << endl;
 cout << "jetELsmear = " << nSmear << endl;
 cout << "jetOversampleFactor = " << jetOversampleFactor << endl;
 cout << "jetMinPt = " << jetMinPt << endl;
 cout << "======= end parameters =======\n";
}

// program parameters, to be read from file
// int nx, ny, nz, eosType ;
// double xmin, xmax, ymin, ymax, zmin, zmax, tau0, tauMax, dtau ;
// char outputDir[255], eosFile[255], chiBfile[255], chiSfile[255] ;
// char icInputFile [255] ;
// double T_ch, mu_b, mu_q, mu_s, gammaS, gammaFactor, exclVolume ;
// int icModel, NPART, glauberVariable=1 ;
// double epsilon0, alpha, impactPar, s0ScaleFactor ;

void readIniPartons(const char* file, vector<Jet*> &jets)
{
 int count;
 for(int ios=0; ios<jetOversampleFactor; ios++){
 ifstream fin(file);
 if (!fin.is_open()) {
  cout << "cannot open initial partons file " << file << endl;
  exit(1);
 }
 count=0;
 while (fin.good()) {  // reading the initial partons, line by line
  string line;
  getline(fin, line);
  istringstream sline(line);
  int type;
  double px, py, pz, E, Q2, x, y, z, t;
  sline >> type >> px >> py >> pz >> E >> Q2 >> x >> y >> z >> t; // new input format!
  Q2 = sqrt(px*px+py*py);  // input Q2 from EPOS ignored
  //Q2 = sqrt(Q2);
  if(Q2!=Q2) Q2 = 0.;
  double eta = 0.5*log((t+z)/(t-z));
  if(eta!=eta) {
   //cout << "etanan" << setw(14) << t << setw(14) << z << endl;
   eta = 0.;
  }
  // ONLY QUARK JETS!
  if(Q2>0.6 and px*px+py*py>jetMinPt*jetMinPt and abs(type)!=9) {
   // skip low-Q partons, qsuch() won't work
   // default/min: 0.6; test 22.0;
   // Event numbering: last 3 digits enumerate the oversampled jet evolution,
   // the upper digits enumerate initial state / hydro configuration.
   jets.push_back(new Jet(1000*eventNo+ios, count, type, px, py, pz, E, Q2, x, y, eta, tau0));
   count++;
  }
 } // file read loop
 fin.close();
 } // oversample loop
 cout << count << " jet partons read in.\n";
}

int main(int argc, char **argv) {
 // pointers to all the main objects
 EoS *eos;
 TransportCoeff *trcoeff;
 Fluid *f;
 Hydro *h;
 time_t start = 0, end;

 time(&start);
 // setting default values for some parameters
 jetOversampleFactor = 1;
 jetMinPt = 10.0;

 // read parameters from file
 char *parFile, *jetParFile, *iniPartonsFile;
 if (argc == 1) {
  cout << "NO PARAMETERS, exiting\n";
  cout << "usage: ./hlle_visc <input file> <jet_input> <eventId>\n";
  exit(1);
 } else if (argc==4){
  parFile = argv[1];
  //jetParFile = argv[2];
  iniPartonsFile = argv[2];
  eventNo = atoi(argv[3]);
 } else {
  cout << "wrong parameter list, exiting\n" ;
  exit(1);
 }
 readParameters(parFile);
 printParameters();

 // EoS
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
 EoS *eosH = new EoSHadron("eos/eosHadronLog.dat");

 // transport coefficients
 trcoeff = new TransportCoeff(etaS, zetaS, eos);

 f = new Fluid(eos, eosH, trcoeff, nx, ny, nz, xmin, xmax, ymin, ymax, etamin,
               etamax, dtau, eCrit, nSmear);
 cout << "fluid allocation done\n";

 // initilal conditions
 if (icModel == 1) {  // optical Glauber
  ICGlauber *ic = new ICGlauber(epsilon0, impactPar, tau0);
  ic->setIC(f, eos);
  delete ic;
 } else if (icModel == 2) {  // Glauber_table + parametrized rapidity dependence
  IC *ic = new IC(icInputFile, s0ScaleFactor);
  ic->setIC(f, eos, tau0);
  delete ic;
 } else if (icModel == 3) {  // UrQMD IC
  IcPartUrqmd *ic = new IcPartUrqmd(f, icInputFile, Rgt, Rgz, tau0);
  ic->setIC(f, eos);
  delete ic;
 } else if (icModel == 4) {  // analytical Gubser solution
  ICGubser *ic = new ICGubser();
  ic->setIC(f, eos, tau0);
  delete ic;
  }else if(icModel==5){ // IC from GLISSANDO + rapidity dependence
   IcGlissando *ic = new IcGlissando(f, icInputFile, tau0, argv[2]);
   ic->setIC(f, eos);
   delete ic;
  }else if(icModel==6){ // EPOS IS
   icEpos::loadIC("ic/epos/hydro_grid_params", "ic/epos/initial_state_table");
   icEpos::setIC(f, eos);
   icEpos::deleteIC();
 } else {
  cout << "icModel = " << icModel << " not implemented\n";
 }
 cout << "IC done\n";

 time_t tinit = 0;
 time(&tinit);
 float diff = difftime(tinit, start);
 cout << "Init time = " << diff << " [sec]" << endl;

 if(icModel==6) { // take tau0 and dt explicitly from EPOS IS output
  tau0 = icEpos::tau0;
  //dtau = icEpos::dt;
 }
 // hydro init
 h = new Hydro(f, eos, trcoeff, tau0, dtau);
 int maxstep = ceil((tauMax - tau0) / dtau);
 start = 0;
 time(&start);
 // h->setNSvalues() ; // initialize viscous terms
 f->initOutput(outputDir, maxstep, tau0, 2);
 //f->outputCorona(tau0);
 // Jet init
 vector<Jet*> jets;
 JetParameters* jparams = new JetParameters(parFile);
 srand(438468301);
 init_tables(438468301);
 jets.clear(); // clear the whole parton vector
 jets.reserve(100);
 readIniPartons(iniPartonsFile, jets);
 string sOutputDir(outputDir);
 ofstream fjetiniout ((sOutputDir+"/jets_initial").c_str());
 for(uint i=0; i<jets.size(); i++) {
  jets[i]->output(fjetiniout);
 }
 fjetiniout.close();

 // disable jet timestep-wise output: consumes too much space
 //ofstream fjetTimeStep((sOutputDir+"/jetTimeSteps").c_str());
 ofstream fjetInterm((sOutputDir+"/jetInermediate").c_str());

 for (int istep = 0; istep < maxstep; istep++) {
  // decrease timestep automatically, but use fixed dtau for output
  int nSubSteps = 1;
  while (dtau / nSubSteps >
         1.0 * (tau0 + dtau * istep) * (etamax - etamin) / (nz - 1))
   nSubSteps *= 2;  // 0.02 in "old" coordinates
  h->setDtau(dtau / nSubSteps);
  // cout<<"dtau = "<<dtau / nSubSteps<<endl;
  for (int j = 0; j < nSubSteps; j++) {
   h->performStep();
  }
  double t = tau0 + istep*dtau; // beginning of timestep (for jet evolution)
  //for(Jet* jet: jets) {
  uint i=0;
  while(i<jets.size()) {
   if(jets[i]->nprts()==0){
    vector<Jet*>::iterator it = jets.begin() + i ;
    jets.erase(it);
    continue;
   }
   jets[i]->makeStep(f, jparams, t, t+dtau, fjetInterm);
   // jets[i]->outputTimestep(t+dtau, i, fjetTimeStep);
   i++;
  };
  f->outputGnuplot(h->getTau());
  //f->outputSurface(h->getTau());
  cout << "tau=" << t << "  done\n";
 }
 //fjetTimeStep.close();
 cout << "finalizing jet branchings...\n";
 for(double tauF=tauMax; tauF<100.; tauF+=1.0) {
  uint i=0;
  while(i<jets.size()) {
   if(jets[i]->nprts()==0){
    vector<Jet*>::iterator it = jets.begin() + i ;
    jets.erase(it);
    continue;
   }
   jets[i]->makeStep(f, jparams, tauF, tauF+1.0, fjetInterm);
   i++;
  };
 }
 cout << "done.\n";
 fjetInterm.close();

 // printing final jets
 ofstream fjetout ((sOutputDir+"/jets_final").c_str());
 ofstream fjetoutGF ((sOutputDir+"/jets_final_GF").c_str());
 ofstream fjettot ((sOutputDir+"/jet_totals").c_str());
 for(uint i=0; i<jets.size(); i++) {
  jets[i]->output(fjetout);
  jets[i]->outputGlobalFrame(fjetoutGF);
  jets[i]->outputTotal(fjettot);
 }
 fjetout.close();
 fjetoutGF.close();
 fjettot.close();

 end = 0;
 time(&end);
 float diff2 = difftime(end, start);
 cout << "Execution time = " << diff2 << " [sec]" << endl;

 delete f;
 delete h;
 delete eos;
 delete eosH;
}
