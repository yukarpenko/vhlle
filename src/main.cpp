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
#include "icGlauberMC.h"
#include "icGubser.h"
#include "icBjorken.h"
#include "eos.h"
#include "eo3.h"
#include "eo1.h"
#include "eoChiral.h"
#include "eoAZH.h"
#include "eoHadron.h"
#include "eoSimpleSpline.h"
#include "eosLinearCombination.h"
#include "eosLinearCombinationTable.h"
#include "eosuQGP.h"
#include "trancoeff.h"
#include "photons.h"
#include "dileptons.h"

using namespace std;

// program parameters, to be read from file
int nx, ny, nz, eosType;
double xmin, xmax, ymin, ymax, etamin, etamax, tau0, tauMax, dtau, taus;
char outputDir[255];
char icInputFile[255];
char eosFile[255];
char eosDir[255];
double etaS, zetaS, eCrit;
double TCrit = -1.;
int icModel,
    glauberVariable =
        1;  // icModel=1 for pure Glauber, 2 for table input (Glissando etc)
double epsilon0, Rgt, Rgz, impactPar, s0ScaleFactor;
double Tcutphot, Tcutdilept;

// Photon spectrum parameters
char photonInputFile[255];
vector<double> pTphot, Yphot;

// Dilepton spectrum parameters
char dileptonInputFile[255];
vector<double> Mdilept, Ydilept;

void readPhotonParams(char *parFile) {
	ifstream fin(parFile);
    if (!fin.is_open()) {
      cout << "cannot open photon parameters file " << parFile << endl;
	  return;
      //exit(1);
    }
	fin >> Tcutphot;
	int entries = 0;
	fin >> entries;
	pTphot.resize(entries);
	Yphot.resize(entries);
	for(int i=0;i<entries;++i) {
		fin >> pTphot[i] >> Yphot[i];
	}
	fin.close();
}

void readDileptonParams(char *parFile) {
	ifstream fin(parFile);
    if (!fin.is_open()) {
      cout << "cannot open dileptons parameters file " << parFile << endl;
	  cout << "no dileptons calculation " << endl;
	  return;
      //exit(1);
    }
	fin >> Tcutdilept;
	int entries = 0;
	fin >> entries;
	Mdilept.resize(entries);
	Ydilept.resize(entries);
	for(int i=0;i<entries;++i) {
		fin >> Mdilept[i] >> Ydilept[i];
	}
	fin.close();
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
	else if (strcmp(parName, "eosFile") == 0)
       strcpy(eosFile, parValue);
	else if (strcmp(parName, "eosDir") == 0)
       strcpy(eosDir, parValue);
    else if (strcmp(parName, "icInputFile") == 0)
      strcpy(icInputFile, parValue);
	else if (strcmp(parName, "photonInputFile") == 0)
      strcpy(photonInputFile, parValue);
	else if (strcmp(parName, "dileptonInputFile") == 0)
      strcpy(dileptonInputFile, parValue);
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
    else if (strcmp(parName, "dtau") == 0)
      dtau = atof(parValue);
	else if (strcmp(parName, "taus") == 0)
      taus = atof(parValue);
    else if (strcmp(parName, "e_crit") == 0)
      eCrit = atof(parValue);
	else if (strcmp(parName, "T_crit") == 0)
      TCrit = atof(parValue);
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
  cout << "taus = " << taus << endl;
  cout << "dtau = " << dtau << endl;
  if (TCrit<0.) cout << "e_crit = " << eCrit << endl;
  else cout << "t_crit = " << TCrit << endl;
  cout << "eta/s = " << etaS << endl;
  cout << "zeta/s = " << zetaS << endl;
  cout << "epsilon0 = " << epsilon0 << endl;
  cout << "Rgt = " << Rgt << "  Rgz = " << Rgz << endl;
  cout << "impactPar = " << impactPar << endl;
  cout << "s0ScaleFactor = " << s0ScaleFactor << endl;
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

int main(int argc, char **argv) {
  // pointers to all the main objects
  EoS *eos;
  TransportCoeff *trcoeff;
  Fluid *f;
  Hydro *h;
  time_t start = 0, end;

  taus = 0.;

  time(&start);

  // read parameters from file
  char *parFile;
  if (argc == 1) {
	parFile = "song2DCPC.b0.ideal.test";
    //cout << "NO PARAMETERS, exiting\n";
    //exit(1);
  } else {
    parFile = argv[1];
  }
  readParameters(parFile);
  printParameters();

  readPhotonParams(photonInputFile);
  readDileptonParams(dileptonInputFile);

  // EoS
  // char * eosfile = "eos/Laine_nf3.dat" ;
  // int ncols = 3, nrows = 286 ;
  // eos = new EoSs(eosfile,ncols) ;
  /*if (eosType == 1)
    eos = new EoSChiral();
  else if (eosType == 2)
    eos = new EoSAZH();
  else {
    cout << "eosType != 1,2\n";
    return 0;
  }*/
  //char eosfile[] = "eos/Laine_nf3.dat";
  char eosfile[] = "eos/Lattice_BW_QGP.dat";
  //char eosfile[] = "eos/Lattice_BW_YM.dat";
  int ncols = 3;
  EoSimpleSpline *eos1 = new EoSimpleSpline((string(eosDir) + "/Lattice_BW_QGP.dat").c_str());
  EoSimpleSpline *eos2 = new EoSimpleSpline((string(eosDir) + "/Lattice_BW_YM.dat").c_str());
  //EoSimpleSpline *eos1 = new EoSuQGP(0., tau0);
  //EoSimpleSpline *eos2 = new EoSuQGP(1.e15, tau0);
  //eos = new EoSs(eosfile, ncols);
  if (eosType==1) eos = new EoSuQGP(taus, tau0);
  else if (eosType==2) {
	//eos = new EoSimpleSpline(eosFile, taus, tau0);
	//eos = new EoSLinearCombination(eos1, eos2, taus, tau0);
	eos = new EoSLinearCombinationTable(eos1, eos2, taus, tau0);
  }
  else {
	eos = new EoSimpleSpline(eosFile);
  }
  //double taus = 0.0;
  //eos = new EoSuQGP(taus, tau0);
  EoS *eosH = new EoSs(eosfile, ncols);
  //EoS *eosH = new EoSHadron("eos/eosHadronLog.dat");
  //EoS *eosH = new EoSHadron("eos/eosHadronLog.dat");

  // transport coefficients
  trcoeff = new TransportCoeff(etaS, zetaS, eos);

  if (TCrit<0.) f = new Fluid(eos, eosH, trcoeff, nx, ny, nz, xmin, xmax, ymin, ymax, etamin,
                etamax, dtau, eCrit);
  else f = new Fluid(eos, eosH, trcoeff, nx, ny, nz, xmin, xmax, ymin, ymax, etamin,
                etamax, dtau, TCrit, true);
  cout << "fluid allocation done\n";

  // initilal conditions
  if(icModel==1){ // optical Glauber
   ICGlauber *ic = new ICGlauber(epsilon0, impactPar, tau0, 208, 6.6, 0.545, 7.0);
   //ic->setIC(f, eos);
   if (glauberVariable==2) // Glauber Monte Carlo
   {
	   ic->init();
	   double rho0 = ic->getrho0();
	   delete ic;
	   ICGlauberMC *ic2 = new ICGlauberMC(epsilon0 / rho0, impactPar, tau0, 208, 6.6, 0.545, 7.0);
	   ic2->setIC(f, eos);
	   delete ic2;
   }
   else {
	   ic->setIC(f, eos);
	   delete ic;
   }
  }else if(icModel==2){ // Glauber_table + parametrized rapidity dependence
   IC *ic = new IC(icInputFile, s0ScaleFactor);
   ic->setIC(f, eos, tau0);
   delete ic;
  }else if(icModel==3){ // UrQMD IC
   IcPartUrqmd *ic = new IcPartUrqmd(f, icInputFile, Rgt, Rgz, tau0);
   ic->setIC(f, eos);
   delete ic;
  }else if(icModel==4){ // analytical Gubser solution
   ICGubser *ic = new ICGubser();
   ic->setIC(f, eos, tau0);
   delete ic;
  }else if(icModel==5){ // analytical Bjorken solution
   ICBjorken *ic = new ICBjorken(epsilon0);
   ic->setIC(f, eos, tau0);
   delete ic;
  }else {
   cout << "icModel = " << icModel << " not implemented\n";
  }
  cout << "IC done\n";

  time_t tinit = 0;
  time(&tinit);
  float diff = difftime(tinit, start);
  cout << "Init time = " << diff << " [sec]" << endl;

  // hydro init
  h = new Hydro(f, eos, trcoeff, tau0, dtau);
  int maxstep = ceil((tauMax - tau0) / dtau);
  start = 0;
  time(&start);
  h->setNSvalues() ; // initialize viscous terms
  h->setQfull();  // set Qfull in each cell, in order to output IC correctly

  //Photons *phot = new PhotonsQGP();
  Photons *phot = new PhotonsAMY(tau0, taus, Tcutphot);
  phot->setMode(0);
  for(int i=0;i<pTphot.size();++i)
	  phot->addPtY(pTphot[i], Yphot[i]);
  //phot->addPtY(0.1, 0.0);
  /*phot->addPtY(0.5, 0.0);
  phot->addPtY(1.0, 0.0);
  phot->addPtY(1.5, 0.0);
  phot->addPtY(2.0, 0.0);
  phot->addPtY(2.5, 0.0);
  phot->addPtY(3.0, 0.0);
  phot->addPtY(3.5, 0.0);
  phot->addPtY(4.0, 0.0);*/
  
  Photons *phot2 = new PhotonsAMY(tau0, taus, Tcutphot);
  phot2->setMode(1);
  for(int i=0;i<pTphot.size();++i)
	  phot2->addPtY(pTphot[i], Yphot[i]);
  
  Dileptons *dilept = new DileptonsQGP(tau0, taus, Tcutdilept);
  for(int i=0;i<Mdilept.size();++i)
	  dilept->addPtY(Mdilept[i], Ydilept[i]);

  
  int photmode = 0;
  if (icModel==5 || (icModel==1 && glauberVariable==1 && fabs(impactPar)<1e-6)) photmode = 1;
  //photmode = 0;
  //f->processPhotons(h->getTau(), dtau, phot, photmode);

  f->initOutput(outputDir, maxstep, tau0, 2);
  f->calcTotals(h->getTau());
  //f->outputCorona(tau0);
#ifndef SWAP_EOS
  //f->outputCorona(tau0);
#endif

  f->outputXTau(h->getTau());
  f->outputYTau(h->getTau());
  f->outputTau(h->getTau());

  clock_t tbeg  = clock();
  clock_t tpure = clock();
  double totpure = 0.;

  for (int istep = 0; istep < maxstep; istep++) {
    // decrease timestep automatically, but use fixed dtau for output
    int nSubSteps = 1;
    while (dtau / nSubSteps > 1.0 * (tau0 + dtau * istep) * (etamax - etamin) / (nz - 1))
      nSubSteps *= 2;  // 0.02 in "old" coordinates
    h->setDtau(dtau / nSubSteps);
    //cout<<"dtau = "<<dtau / nSubSteps<<endl;
    for (int j = 0; j < nSubSteps; j++) {
      h->performStep();
	  f->processPhotons(h->getTau(), dtau, phot, photmode);
	  f->processPhotons(h->getTau(), dtau, phot2, photmode);
	  f->processDileptons(h->getTau(), dtau, dilept, photmode);
    }
	totpure += (clock() - tpure);
	cout << "step= " << istep << "  tau= " << h->getTau() << "  dtau= " << dtau / nSubSteps << "\n"
         << endl; 
    //f->outputGnuplot(h->getTau());
	f->calcTotals(h->getTau());
	f->outputXTau(h->getTau());
	f->outputYTau(h->getTau());
	f->outputTau(h->getTau());
	//f->processPhotons(h->getTau(), dtau, phot, photmode);
#ifndef SWAP_EOS
    //f->outputSurface(h->getTau());
#endif
	{
		  PhotonSpectrum phspec = phot->GetSpectrum();
		  cout << "Photon spectrum 1" << endl;
		  cout << setw(14) << "pT(GeV)";
		  cout << setw(14) << "y";
		  cout << setw(14) << "yield";
		  cout << setw(14) << "v1";
		  cout << setw(14) << "v2";
		  cout << setw(14) << "v3";
		  cout << endl;
		  for(int i=0;i<phspec.Entries.size();++i) {
			  cout << setw(14) << phspec.Entries[i].pt;
			  cout << setw(14) << phspec.Entries[i].y;
			  cout << setw(14) << phspec.Entries[i].yield;
			  cout << setw(14) << phspec.Entries[i].v1;
			  cout << setw(14) << phspec.Entries[i].v2;
			  cout << setw(14) << phspec.Entries[i].v3;
			  cout << endl;
		  }
	}
	{
		  PhotonSpectrum phspec = phot2->GetSpectrum();
		  cout << "Photon spectrum 2" << endl;
		  cout << setw(14) << "pT(GeV)";
		  cout << setw(14) << "y";
		  cout << setw(14) << "yield";
		  cout << setw(14) << "v1";
		  cout << setw(14) << "v2";
		  cout << setw(14) << "v3";
		  cout << endl;
		  for(int i=0;i<phspec.Entries.size();++i) {
			  cout << setw(14) << phspec.Entries[i].pt;
			  cout << setw(14) << phspec.Entries[i].y;
			  cout << setw(14) << phspec.Entries[i].yield;
			  cout << setw(14) << phspec.Entries[i].v1;
			  cout << setw(14) << phspec.Entries[i].v2;
			  cout << setw(14) << phspec.Entries[i].v3;
			  cout << endl;
		  }
	}
	{
		  DileptonSpectrum dilspec = dilept->GetSpectrum();
		  cout << "Dilepton spectrum" << endl;
		  cout << setw(14) << "M(GeV)";
		  cout << setw(14) << "y";
		  cout << setw(14) << "yield";
		  cout << setw(14) << "v2";
		  cout << endl;
		  for(int i=0;i<dilspec.Entries.size();++i) {
			  cout << setw(14) << dilspec.Entries[i].M;
			  cout << setw(14) << dilspec.Entries[i].y;
			  cout << setw(14) << dilspec.Entries[i].yield;
			  cout << setw(14) << dilspec.Entries[i].v2;
			  cout << endl;
		  }
	}
	tpure = clock();
  }

  end = 0;
  time(&end);
  float diff2 = difftime(end, start);
  cout << "Execution time  = " << diff2 << " [sec]" << endl;
  cout << "       Run time = " << (clock() - tbeg) / (double)(CLOCKS_PER_SEC) << " [sec]" << endl;
  cout << "Simulation time = " << totpure / (double)(CLOCKS_PER_SEC) << " [sec]" << endl;
  cout << "  Time per step = " << totpure / (double)(CLOCKS_PER_SEC) / maxstep << " [sec]" << endl;

  PhotonSpectrum phspec = phot->GetSpectrum();

  cout << "Photon spectrum 1" << endl;
  cout << setw(14) << "pT(GeV)";
  cout << setw(14) << "y";
  cout << setw(14) << "yield";
  cout << setw(14) << "v1";
  cout << setw(14) << "v2";
  cout << setw(14) << "v3";
  cout << endl;
  for(int i=0;i<phspec.Entries.size();++i) {
	  cout << setw(14) << phspec.Entries[i].pt;
	  cout << setw(14) << phspec.Entries[i].y;
	  cout << setw(14) << phspec.Entries[i].yield;
	  cout << setw(14) << phspec.Entries[i].v1;
	  cout << setw(14) << phspec.Entries[i].v2;
	  cout << setw(14) << phspec.Entries[i].v3;
	  cout << endl;
  }
  
  PhotonSpectrum phspec2 = phot2->GetSpectrum();
  
  cout << "Photon spectrum 2" << endl;
  cout << setw(14) << "pT(GeV)";
  cout << setw(14) << "y";
  cout << setw(14) << "yield";
  cout << setw(14) << "v1";
  cout << setw(14) << "v2";
  cout << setw(14) << "v3";
  cout << endl;
  for(int i=0;i<phspec2.Entries.size();++i) {
	  cout << setw(14) << phspec2.Entries[i].pt;
	  cout << setw(14) << phspec2.Entries[i].y;
	  cout << setw(14) << phspec2.Entries[i].yield;
	  cout << setw(14) << phspec2.Entries[i].v1;
	  cout << setw(14) << phspec2.Entries[i].v2;
	  cout << setw(14) << phspec2.Entries[i].v3;
	  cout << endl;
  }
  
  DileptonSpectrum dilspec = dilept->GetSpectrum();
  cout << "Dilepton spectrum" << endl;
  cout << setw(14) << "M(GeV)";
  cout << setw(14) << "y";
  cout << setw(14) << "yield";
  cout << setw(14) << "v2";
  cout << endl;
  for(int i=0;i<dilspec.Entries.size();++i) {
	  cout << setw(14) << dilspec.Entries[i].M;
	  cout << setw(14) << dilspec.Entries[i].y;
	  cout << setw(14) << dilspec.Entries[i].yield;
	  cout << setw(14) << dilspec.Entries[i].v2;
	  cout << endl;
  }

  double mn = 1.;
  double Rn = 6.5;
  if (nx==5 && ny==5) // No transverse dynamics
  {
	  mn = C_PI * Rn * Rn / f->getDx() / f->getDy();
  }

  {
	  fstream fout((string(outputDir) + "-photons-1.dat").c_str(), fstream::out);

	  fout << setw(14) << "pT(GeV)";
	  fout << setw(14) << "y";
	  fout << setw(14) << "yield";
	  fout << setw(14) << "v1";
	  fout << setw(14) << "v2";
	  fout << setw(14) << "v3";
	  fout << endl;
	  for(int i=0;i<phspec.Entries.size();++i) {
		  fout << setw(14) << phspec.Entries[i].pt;
		  fout << setw(14) << phspec.Entries[i].y;
		  fout << setw(14) << phspec.Entries[i].yield * mn;
		  fout << setw(14) << phspec.Entries[i].v1;
		  fout << setw(14) << phspec.Entries[i].v2;
		  fout << setw(14) << phspec.Entries[i].v3;
		  fout << endl;
	  }
	  fout.close();

	  delete phot;
  }
  
  {
	  fstream fout((string(outputDir) + "-photons-2.dat").c_str(), fstream::out);

	  fout << setw(14) << "pT(GeV)";
	  fout << setw(14) << "y";
	  fout << setw(14) << "yield";
	  fout << setw(14) << "v1";
	  fout << setw(14) << "v2";
	  fout << setw(14) << "v3";
	  fout << endl;
	  for(int i=0;i<phspec2.Entries.size();++i) {
		  fout << setw(14) << phspec2.Entries[i].pt;
		  fout << setw(14) << phspec2.Entries[i].y;
		  fout << setw(14) << phspec2.Entries[i].yield * mn;
		  fout << setw(14) << phspec2.Entries[i].v1;
		  fout << setw(14) << phspec2.Entries[i].v2;
		  fout << setw(14) << phspec2.Entries[i].v3;
		  fout << endl;
	  }
	  fout.close();

	  delete phot2;
  }
  
  {
	  fstream fout((string(outputDir) + "-dileptons.dat").c_str(), fstream::out);

	  fout << setw(14) << "M(GeV)";
	  fout << setw(14) << "y";
	  fout << setw(14) << "yield";
	  fout << setw(14) << "v2";
	  fout << endl;
	  for(int i=0;i<dilspec.Entries.size();++i) {
		  fout << setw(14) << dilspec.Entries[i].M;
		  fout << setw(14) << dilspec.Entries[i].y;
		  fout << setw(14) << dilspec.Entries[i].yield * mn;
		  fout << setw(14) << dilspec.Entries[i].v2;
		  fout << endl;
	  }
	  fout.close();
  }

  delete dilept;

  delete f;
  delete h;
  delete eos;
  delete eosH;
  delete eos1;
  delete eos2;
}
