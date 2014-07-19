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
#include "eoAZH.h"
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
double epsilon0, Rgauss, impactPar, s0ScaleFactor ;


int main(int argc, char **argv)
{
  TApplication theApp("App", &argc, argv);
  // pointers to all the main objects
  EoS *eos ;

  eos = new EoSAZH() ;
  
  const int N = 5000 ;
  double x[N], y[N] ;
  int npoints=0 ;
  for(double e=0.; e<3.0; e+=0.01){
   double Tt, mubt, muqt, must, pt ;
   eos->eos(e, 0., 0., 0., Tt, mubt, muqt, must, pt) ;
   x[npoints] = e;
   y[npoints] = Tt ;
   npoints++ ;
  }
  TGraph *g = new TGraph(npoints,x,y) ;
  g->SetMarkerStyle(22) ;
  g->SetMarkerSize(0.8) ;
  g->Draw("AP") ;
  theApp.Run() ;

  delete eos ;
}
