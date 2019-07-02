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
#include <TFile.h>
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
#include "../JT/const.h"
#include "../JT/params.h"
#include "../JT/cascade.h"
//--- hadron sampling
#include "../HS/gen.h"
#include "../HS/cascade.h"
#include "../HS/particle.h"
#include "../HS/params.h"
#include "../HS/tree.h"
//--- PYTHIA
#include "Pythia8/Pythia.h"
#include "../JT/interfacePythia.h"
//--- ROOT
#include <TRandom3.h>

using namespace std;

// program parameters, to be read from file
int nx, ny, nz, eosType;
double xmin, xmax, ymin, ymax, etamin, etamax, tau0, tauMax, tauResize, dtau;
char outputDir[255];
char icInputFile[255], icGridParams[255], iniHardPartons[255];

double etaS, zetaS, eCrit;
int icModel, nSmear = 1,
    glauberVariable =
        1;  // icModel=1 for pure Glauber, 2 for table input (Glissando etc)
double epsilon0, Rgt, Rgz, impactPar, s0ScaleFactor;
// ##### jet-related parameters
int eventNo; // enumerates the initial/input configuration
int jetOversampleFactor, nJetSubSteps;
double jetMinPt, ptTrigger;

extern TRandom3 *rnd;

void setDefaultParameters() {
 tauResize = 4.0;
 ptTrigger = 10.0;
}

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
  else if (strcmp(parName, "icGridParams") == 0)
   strcpy(icGridParams, parValue);
  else if (strcmp(parName, "iniHardPartons") == 0)
   strcpy(iniHardPartons, parValue);
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
  // ### JET-related parameters
  else if (strcmp(parName, "jetOversampleFactor") == 0)
   jetOversampleFactor = atoi(parValue);
  else if (strcmp(parName, "nJetSubSteps") == 0)
   nJetSubSteps = atoi(parValue);
  else if (strcmp(parName, "jetMinPt") == 0)
   jetMinPt = atof(parValue);
  else if (strcmp(parName, "ptTrigger") == 0)
   ptTrigger = atof(parValue);
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
 cout << "icInputFile = " << icInputFile << endl;
 cout << "icGridParams = " << icGridParams << endl;
 cout << "iniHardPartons = " << iniHardPartons << endl;
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
 cout << "jetELsmear = " << nSmear << endl;
 cout << "jetOversampleFactor = " << jetOversampleFactor << endl;
 cout << "nJetSubSteps = " << nJetSubSteps << endl;
 cout << "jetMinPt = " << jetMinPt << endl;
 cout << "ptTrigger = " << ptTrigger << endl;
 cout << "======= end parameters =======\n";
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
   f->getZ(0), f->getZ(f->getNZ()-1), 2.0*h->getDtau(), f->geteCrit(), nSmear);
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

namespace IniPartons {

vector<double> xo, yo, etao, rapo; // original list of parton coordinates
vector<int> typeo; // original list of parton types
vector<std::pair<double, int>> ptOrder; // pt-index map
bool discardOrigPartons; // flag to discard the original partons
                         // if they don't match the pT trigger criterion

void readIniPartons(const char* file, vector<Jet*> &jets, vector<vector<Jet*> > &jetEvents)
{
 jetEvents.resize(jetOversampleFactor);
 discardOrigPartons = false;
 double ptMax = 0.0;
 ifstream fin(file);
 if (!fin.is_open()) {
  cout << "cannot open initial partons file " << file << endl;
  exit(1);
 }
 int count=0;
 xo.clear();
 yo.clear();
 etao.clear();
 rapo.clear();
 typeo.clear();
 ptOrder.clear();
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
   // --- adding the original jet partons
   const int ios = 0; // later, in the oversamplig loop, ios starts from 1.
   Jet* _jet = new Jet(1000*eventNo+ios, count, type, px, py, pz, E, Q2, x, y, eta, tau0);
   jets.push_back(_jet);
   jetEvents[ios].push_back(_jet);
   ptMax = std::max(ptMax, Q2);
   // --- filling the coordinate arrays for the oversampling procedure
   xo.push_back(x);
   yo.push_back(y);
   etao.push_back(eta);
   double rap = 0.5*log((E+pz)/(E-pz));
   if(rap!=rap) rap = 20.0;
   rapo.push_back(rap);
   typeo.push_back(type);
   std::pair<double, int> _pair (Q2, count);
   ptOrder.push_back(_pair);
   cout << "orig_pt" << setw(14) << Q2 << setw(14) << x << setw(14) << y <<endl;
   count++;
  }
 } // file read loop
 fin.close();
 cout << count << " jet partons read in.\n";
 // if the original hard partons from EPOS don't meet the pT trigger condition,
 // set the flag to discard them and replace with another oversampled distribution.
 if(ptMax<ptTrigger) discardOrigPartons = true;
 // pT sorting
 std::sort(ptOrder.begin(), ptOrder.end());
 for(int i=0; i<ptOrder.size(); i++) {
  cout << "sorted " << setw(14) << ptOrder[i].first << setw(14) << xo[ptOrder[i].second]
     << setw(14) << yo[ptOrder[i].second] << endl;
 }
}

void oversampleIniPartons(vector<Jet*> &jets, vector<vector<Jet*> > &jetEvents)
{
 int startOS;
 if(discardOrigPartons) {
  // if the flag is set, clear the jet vectors to discard the original partons
  jets.clear();
  jetEvents.clear();
  jetEvents.resize(jetOversampleFactor);
  startOS = 0;
 } else {
  startOS = 1;
 }
 for(int ios=startOS; ios<jetOversampleFactor; ios++) { // oversampling loop, STARTING FROM 1
  cout << "oversampling event " << ios << endl;
  vector<double> ptPool;
  double ptMax = 0.0;
  do {
  ptPool.clear();
   for(int i=0; i<ptOrder.size(); i++) {
    double _pt = jetMinPt/pow(1.0 - rnd->Rndm(), 1.0/(5.3-1.0)); // 5.3 is the power law
    ptPool.push_back(_pt);
    ptMax = std::max(ptMax, _pt);
   }
  } while(ptMax<ptTrigger);
  std::sort(ptPool.begin(), ptPool.end());
  cout << "oversampled_pt_max:" << setw(14) << ptPool[ptPool.size()-1] << endl;
  for(int i=0; i<ptPool.size(); i++) {
   const double _phi = 2.0*M_PI*rnd->Rndm();
   const double _px = ptPool[i]*cos(_phi);
   const double _py = ptPool[i]*sin(_phi);
   const double _pz = ptPool[i]*sinh(rapo[ptOrder[i].second]);
   const double _E = ptPool[i]*cosh(rapo[ptOrder[i].second]);
   if(i==ptPool.size()-1) {
    cout << setw(14) << xo[ptOrder[i].second] << setw(14) << yo[ptOrder[i].second] << endl;
   }
   // Event numbering: last 3 digits enumerate the oversampled jet evolution,
   // the upper digits enumerate initial state / hydro configuration.
   Jet* _jet = new Jet(1000*eventNo+ios, i, typeo[ptOrder[i].second],
    _px, _py, _pz, _E, ptPool[i], xo[ptOrder[i].second], yo[ptOrder[i].second],
    etao[ptOrder[i].second], tau0);
   jets.push_back(_jet);
   jetEvents[ios].push_back(_jet);
  }
 } // oversampling loop
}

} // end namespace IniPartons

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
 nJetSubSteps = 2;

 // read parameters from file
 char *parFile;
 if (argc==3) {
  parFile = argv[1];
  //jetParFile = argv[1];
  eventNo = atoi(argv[2]);
 } else {
  cout << "usage: ./hlle_visc <input file> <eventId>\n";
  exit(1);
 }
 setDefaultParameters();
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
   icEpos::loadIC(icGridParams, icInputFile);
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
 }
 // hydro init
 h = new Hydro(f, eos, trcoeff, tau0, dtau);
 start = 0;
 time(&start);
 // h->setNSvalues() ; // initialize viscous terms
 f->initOutput(outputDir, tau0);
 //f->outputCorona(tau0);
 // Jet init
 vector<Jet*> jets; // all jets in a plain vector
 vector<vector<Jet*> > jetEvents; // jets separated into jet events
 JetParameters* jparams = new JetParameters(parFile);
 srand(438468301);
 init_tables(438468301);
 jets.clear(); // clear the whole parton vector
 jets.reserve(100);
 IniPartons::readIniPartons(iniHardPartons, jets, jetEvents);
 IniPartons::oversampleIniPartons(jets, jetEvents);
 string sOutputDir(outputDir);
 ofstream fjetiniout ((sOutputDir+"/jets_initial").c_str());
 for(uint i=0; i<jets.size(); i++) {
  jets[i]->output(fjetiniout);
 }
 fjetiniout.close();
 // hadron sampler (particlization) init
 HSparams::readParams(parFile);
 HSparams::NEVENTS = jetOversampleFactor;
 gen::init();
 // interface to PYTHIA
 InterfacePythia pyInt;

 // disable jet timestep-wise output: consumes too much space
 //ofstream fjetTimeStep((sOutputDir+"/jetTimeSteps").c_str());
 ofstream fjetInterm((sOutputDir+"/jets_intermediate_GF").c_str());

 bool resized = false; // flag if the hydro grid has been resized
 do {
  // small tau: decrease timestep by makins substeps, in order
  // to avoid instabilities in eta direction (signal velocity ~1/tau)
  int nSubSteps = 1;
  while (dtau / nSubSteps >
         1.0 * h->getTau() * (etamax - etamin) / (nz - 1)) {
   nSubSteps *= 2;  // 0.02 in "old" coordinates
  }
  if(nSubSteps>1) {
   h->setDtau(h->getDtau() / nSubSteps);
   for (int j = 0; j < nSubSteps; j++)
    h->performStep();
   h->setDtau(h->getDtau() * nSubSteps);
   cout << "timestep reduced by " << nSubSteps << endl;
  } else
   h->performStep();
  uint i=0;
  while(i<jets.size()) {
   if(jets[i]->nprts()==0){
    vector<Jet*>::iterator it = jets.begin() + i ;
    jets.erase(it);
    continue;
   }
   for(int ijst=0; ijst<nJetSubSteps; ijst++) {
    const double jt0 = h->getTau()-h->getDtau();
    const double jdt = h->getDtau() / nJetSubSteps;
    jets[i]->makeStep(f, jparams, jt0+ijst*jdt, jt0+(ijst+1)*jdt, fjetInterm);
   }
   // jets[i]->outputTimestep(t+dtau, i, fjetTimeStep);
   i++;
  };
  f->outputGnuplot(h->getTau());
  f->outputSurface(h->getTau());
  if(h->getTau()>=tauResize and resized==false) {
   cout << "grid resize\n";
   f = expandGrid2x(h, eos, eosH, trcoeff);
   resized = true;
  }
 } while(h->getTau()<tauMax+0.0001);
 cout << "finalizing jet branchings...\n";
 jparams->nullifyTransportCoeff();
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

 // sampling final state medium hadrons
 gen::generate();
 // decaying resonances and filling ROOT trees
 TFile *outputFile = new TFile("output.root", "RECREATE"); 
 outputFile->cd();
 MyTree *tree = new MyTree("tree") ;
 // Cooper-Frye oversampling loop
 for(int iev=0; iev<HSparams::NEVENTS; iev++){
  tree->clear();
  tree->addMediumHadrons(iev);
  gen::urqmd(iev);   // it only does resonance decays
  tree->addMediumHadrons(iev);
  vector<gen::Particle> jetPartons, jetHadrons;
  jetPartons.clear();
  jetHadrons.clear();
  for(uint i=0; i<jetEvents[iev].size(); i++) {
   pyInt.do1JetHadronization(*(jetEvents[iev][i]), jetPartons, jetHadrons);
   tree->addJetParticles(jetPartons);   // assuming jetOversampling==1
   tree->addJetParticles(jetHadrons);
  }
  tree->fillTree();
 } // end events loop
 outputFile->Write() ;
 outputFile->Close() ;

 // printing final jets
 ofstream fjetout ((sOutputDir+"/jets_final").c_str());
 ofstream fjetoutGF ((sOutputDir+"/jets_final_GF").c_str());
 ofstream fjettot ((sOutputDir+"/jet_totals").c_str());
 ofstream fjetBothFr ((sOutputDir+"/jets_both_frames").c_str());
 for(uint i=0; i<jets.size(); i++) {
  jets[i]->output(fjetout);
  jets[i]->outputGlobalFrameColour(fjetoutGF);
  jets[i]->outputTotal(fjettot);
  jets[i]->outputBothFrames(fjetBothFr);
 }
 fjetout.close();
 fjetoutGF.close();
 fjettot.close();
 fjetBothFr.close();

 end = 0;
 time(&end);
 float diff2 = difftime(end, start);
 cout << "Execution time = " << diff2 << " [sec]" << endl;

 delete f;
 delete h;
 delete eos;
 delete eosH;
}
