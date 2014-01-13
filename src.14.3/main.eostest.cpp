#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <sstream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>

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
#include "rmn.h"

using namespace std ;

// program parameters, to be read from file
int nx, ny, nz, eosType ;
double xmin, xmax, ymin, ymax, etamin, etamax, tau0, tauMax, dtau ;
char outputDir[255], eosFile[255], chiBfile[255], chiSfile[255] ;
char icInputFile [255] ;
double T_ch, mu_b, mu_q, mu_s, gammaS, gammaFactor, exclVolume, etaS, zetaS, eCrit ;
int icModel, glauberVariable=1 ; // icModel=1 for pure Glauber, 2 for table input (Glissando etc)
double epsilon0, alpha, impactPar, s0ScaleFactor ;


int main(int argc, char **argv)
{
  TApplication theApp("App", &argc, argv);
  // pointers to all the main objects
  EoS *eos ;

  char * eosfile = "eos/Laine_nf3.dat" ;
  int ncols = 3, nrows = 286 ;
  //eos = new EoSs(eosfile,ncols) ;
  eos = new EoSChiral() ;
  EoS* eosH = new EoSHadron("eos/eosHadron3D.dat") ;
  
  // ==== test, Jan 2, 2014
  double e, p, nb, nq, ns, vx, vy, vz ;
  double Q [7] = {5.84382/0.82, 0.2551/0.82, -0.09191/0.82, 3.10939/0.82, 0.555881/0.82, 0.06406486/0.82, 0} ;
  transformPV(eosH, Q, e, p, nb, nq, ns, vx, vy, vz) ;
  cout << "test passed\n" ;
  return 0 ;
  
  
  // ==== previous tests ====
  //double e=0.5, p, nb=0.2, nq=0., ns=0., vx=0.8, vy=0., vz=0. ;
  //double Q [7] ;
  transformCV(e, eos->p(e,nb,nq,ns), nb, nq, ns, vx, vy, vz, Q) ;
  transformPV(eosH, Q, e, p, nb, nq, ns, vx, vy, vz) ;
  cout<<"HelloWorld;\n" ;
  transformCV(e, eosH->p(e,nb,nq,ns), nb, nq, ns, vx, vy, vz, Q) ;
  cout<<"HelloWorld;\n" ;

  const int N = 500 ;
  double x[N], y[N] ;
  for(int i=0; i<N; i++){
   double Tt, mubt, muqt, must, pt ;
   double e = i*1.0/N ;
   double nb = i*0.2/N ;
   double nq = i*0.1/N ;
   eosH->eos(0.5, nb, 0., 0., Tt, mubt, muqt, must, pt) ;
   x[i] = nb;
   y[i] = mubt ;
  }
  TGraph *g = new TGraph(N,x,y) ;
  g->SetMarkerStyle(22) ;
  g->SetMarkerSize(0.8) ;
  g->Draw("AP") ;
  theApp.Run() ;

  delete eos ;
  delete eosH ;
}
