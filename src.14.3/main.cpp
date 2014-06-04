#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <sstream>
//#include <TMath.h>
#include "fld.h"
#include "hdo.h"
#include "ic.h"
#include "ickw.h"
#include "eos.h"
#include "eo3.h"
#include "eo1.h"
#include "eoChiral.h"
#include "eoHadron.h"
#include "trancoeff.h"

using namespace std ;

// program parameters, to be read from file
int nx, ny, nz, eosType ;
double xmin, xmax, ymin, ymax, etamin, etamax, tau0, tauMax, dtau ;
char outputDir[255], eosFile[255], chiBfile[255], chiSfile[255] ;
char icInputFile [255] ;
double T_ch, mu_b, mu_q, mu_s, gammaS, gammaFactor, exclVolume, etaS, zetaS, eCrit ;
int icModel, glauberVariable=1 ; // icModel=1 for pure Glauber, 2 for table input (Glissando etc)
double epsilon0, alpha, impactPar, s0ScaleFactor ;

void readParameters(char *parFile)
{
	char parName [255], parValue [255] ;
	ifstream fin(parFile) ;
	if(!fin.is_open()) { cout << "cannot open parameters file " << parFile << endl ; exit(1) ; }
	while(fin.good()){
	 string line ;
	 getline(fin, line) ;
	 istringstream sline (line) ;
	 sline >> parName >> parValue ;
	 if     (strcmp(parName,"outputDir")==0) strcpy(outputDir, parValue) ;
	 else if(strcmp(parName,"eosType")==0) eosType = atoi(parValue) ;
	 else if(strcmp(parName,"eosFile")==0) strcpy(eosFile, parValue) ;
	 else if(strcmp(parName,"chiBfile")==0) strcpy(chiBfile, parValue) ;
	 else if(strcmp(parName,"chiSfile")==0) strcpy(chiSfile, parValue) ;
	 else if(strcmp(parName,"icInputFile")==0) strcpy(icInputFile, parValue) ;
	 else if(strcmp(parName,"nx")==0) nx = atoi(parValue) ;
	 else if(strcmp(parName,"ny")==0) ny = atoi(parValue) ;
	 else if(strcmp(parName,"nz")==0) nz = atoi(parValue) ;
	 else if(strcmp(parName,"icModel")==0) icModel = atoi(parValue) ;
	 else if(strcmp(parName,"glauberVar")==0) glauberVariable = atoi(parValue) ;
	 else if(strcmp(parName,"xmin")==0) xmin = atof(parValue) ;
	 else if(strcmp(parName,"xmax")==0) xmax = atof(parValue) ;
	 else if(strcmp(parName,"ymin")==0) ymin = atof(parValue) ;
	 else if(strcmp(parName,"ymax")==0) ymax = atof(parValue) ;
	 else if(strcmp(parName,"etamin")==0) etamin = atof(parValue) ;
	 else if(strcmp(parName,"etamax")==0) etamax = atof(parValue) ;
	 else if(strcmp(parName,"tau0")==0) tau0 = atof(parValue) ;
	 else if(strcmp(parName,"tauMax")==0) tauMax = atof(parValue) ;
	 else if(strcmp(parName,"dtau")==0) dtau = atof(parValue) ;
	 else if(strcmp(parName,"T_ch")==0) T_ch = atof(parValue) ;
	 else if(strcmp(parName,"mu_b")==0) mu_b = atof(parValue) ;
	 else if(strcmp(parName,"mu_q")==0) mu_q = atof(parValue) ;
	 else if(strcmp(parName,"mu_s")==0) mu_s = atof(parValue) ;
   else if(strcmp(parName,"e_crit")==0) eCrit = atof(parValue) ;
   else if(strcmp(parName,"etaS")==0) etaS = atof(parValue) ;
   else if(strcmp(parName,"zetaS")==0) zetaS = atof(parValue) ;
	 else if(strcmp(parName,"gammaS")==0) gammaS = atof(parValue) ;
	 else if(strcmp(parName,"gammaFactor")==0) gammaFactor = atof(parValue) ;
	 else if(strcmp(parName,"exclVolume")==0) exclVolume = atof(parValue) ;
	 else if(strcmp(parName,"epsilon0")==0) epsilon0 = atof(parValue) ;
	 else if(strcmp(parName,"alpha")==0) alpha = atof(parValue) ;
	 else if(strcmp(parName,"impactPar")==0) impactPar = atof(parValue) ;
	 else if(strcmp(parName,"s0ScaleFactor")==0) s0ScaleFactor = atof(parValue) ;
	 else if(parName[0]=='!') cout << "CCC " << sline.str() << endl ;
	 else cout << "UUU " << sline.str() << endl ;
	}
}


void printParameters()
{
  cout << "====== parameters ======\n" ;
  cout << "outputDir = " << outputDir << endl ;
  cout << "eosType = " << eosType << endl ;
  cout << "eosFile = " << eosFile << endl ;
  cout << "chiBfile = " << chiBfile << endl ;
  cout << "chiSfile = " << chiSfile << endl ;
  cout << "nx = " <<  nx << endl ;
  cout << "ny = " << ny << endl ;
  cout << "nz = " << nz << endl ;
  cout << "icModel = " << icModel << endl ;
  cout << "glauberVar = " << glauberVariable << "   ! 0=epsilon,1=entropy" << endl ;
  cout << "xmin = " << xmin << endl ;
  cout << "xmax = " << xmax << endl ;
  cout << "ymin = " << ymin << endl ;
  cout << "ymax = " << ymax << endl ;
  cout << "etamin = " << etamin << endl ;
  cout << "etamax = " << etamax << endl ;
  cout << "tau0 = " << tau0 << endl ;
  cout << "tauMax = " << tauMax << endl ;
  cout << "dtau = " << dtau << endl ;
  cout << "T_ch = " << T_ch << endl ;
  cout << "mu_b = " << mu_b << endl ;
  cout << "mu_q = " << mu_q << endl ;
  cout << "mu_s = " << mu_s << endl ;
  cout << "e_crit = " << eCrit << endl ;
  cout << "eta/s = " << etaS << endl ;
  cout << "zeta/s = " << zetaS << endl ;
  cout << "gammaS = " << gammaS << endl ;
  cout << "gammaFactor = " << gammaFactor << endl ;
  cout << "exclVolume = " << exclVolume << endl ;
  cout << "epsilon0 = " << epsilon0 << endl ;
  cout << "alpha = " << alpha << endl ;
  cout << "impactPar = " << impactPar << endl ;
  cout << "s0ScaleFactor = " << s0ScaleFactor << endl ;
  cout << "======= end parameters =======\n" ;
}


// program parameters, to be read from file
//int nx, ny, nz, eosType ;
//double xmin, xmax, ymin, ymax, zmin, zmax, tau0, tauMax, dtau ;
//char outputDir[255], eosFile[255], chiBfile[255], chiSfile[255] ;
//char icInputFile [255] ;
//double T_ch, mu_b, mu_q, mu_s, gammaS, gammaFactor, exclVolume ;
//int icModel, NPART, glauberVariable=1 ;
//double epsilon0, alpha, impactPar, s0ScaleFactor ;


int main(int argc, char **argv)
{
  // pointers to all the main objects
  EoS *eos ;
  TransportCoeff *trcoeff ;
  Fluid *f ;
  Hydro* h ;
  time_t start=0, end ;
  
  time(&start);
  int memory = 0 ;
  //checkmemory_(&memory) ;

  // read parameters from file
	char* parFile ;
	if(argc==1){
		cout << "NO PARAMETERS, exiting\n" ;
		exit(1) ;
	}else{
		parFile = argv[1] ;
	}
	readParameters(parFile) ;
	printParameters() ;
  
  // EoS
  char * eosfile = "eos/Laine_nf3.dat" ;
  int ncols = 3, nrows = 286 ;
  eos = new EoSs(eosfile,ncols) ;
  //eos = new EoSChiral() ;
  EoS* eosH = new EoSHadron("eos/eosHadron3D.dat") ;
  
  
  // transport coefficients
  trcoeff = new TransportCoeff(etaS, zetaS, eos) ;
  
  
  f = new Fluid(eos, eosH, trcoeff, nx, ny, nz, xmin, xmax, ymin, ymax, etamin, etamax, dtau, eCrit) ;
  //checkmemory_(&memory) ;
  cout << "fluid allocation done\n" ;

  // initilal conditions
  IC *ic = new IC(epsilon0, impactPar, alpha) ; // pure Glauber
  //IC_KW *ic = new IC_KW("ic/zzz6.ico") ;
  ic->setIC(f,eos, tau0) ;
  delete ic ;
	cout<<"IC done\n" ;
	//checkmemory_(&memory) ;
	

  time_t tinit=0 ;
  time(&tinit); float diff = difftime(tinit, start);
  cout<<"Init time = "<<diff<< " [sec]" << endl;
  
  // hydro init
  h = new Hydro(f, eos, trcoeff, tau0, dtau) ;
  int maxstep = ceil((tauMax-tau0)/dtau) ;
  start = 0;
  time(&start);
  h->setNSvalues() ; // initialize viscous terms
  h->setQfull() ; // set Qfull in each cell, in order to output IC correctly

  // hllev321v1 = with pre-advection
  f->initOutput(outputDir, maxstep, tau0, 2) ;
  f->calcTotals(h->getTau()) ;

	for(int istep=0; istep<maxstep; istep++)
	{
    // decrease timestep automatically, but use fixed dtau for output
    int nSubSteps = 1 ;
    while(dtau/nSubSteps>0.5*(tau0+dtau*istep)) nSubSteps *= 2 ; // 0.02 in "old" coordinates
    h->setDtau(dtau/nSubSteps) ;
    for(int j=0; j<nSubSteps; j++){ h->performStep() ;
    }
    cout << "step= " << istep << "  dtau= " << dtau/nSubSteps << "\n" << endl ; // "\r" << flush
    f->outputPDirections(h->getTau());
    f->calcTotals(h->getTau()) ;
	}

  end=0 ;
  time(&end); float diff2 = difftime(end, start);
  cout<<"Execution time = "<<diff2<< " [sec]" << endl;
    
	//checkmemory_(&memory) ;
	
  delete f ;
  delete h ;
  delete eos ;
  delete eosH ;
  //checkmemory_(&memory) ;

	//checkmemory_(&memory) ;
}
