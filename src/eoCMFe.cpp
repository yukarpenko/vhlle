#include <math.h>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>


#include "eos.h"
#include "eoCMFe.h"


/*
Equation of state based on 
J. Steinheimer, S. Schramm, and H. Stocker, Phys. Rev.
C 84, 045208 (2011), arXiv:1108.2596 [hep-ph].
A. Motornenko, J. Steinheimer, V. Vovchenko,
S. Schramm, and H. Stoecker, Phys. Rev. C 101,
034904 (2020), arXiv:1905.00866 [hep-ph].
E. R. Most, A. Motornenko, J. Steinheimer, V. Dex-
heimer, M. Hanauske, L. Rezzolla, and H. Stoecker,
Phys. Rev. D 107, 043034 (2023), arXiv:2201.13150
[nucl-th]
*/

using namespace std;


class EoSCMFeaux {
  const double EPS_0 = 146.51751415742*0.87326281;
  const double N_0 = 0.15891*0.87272727;
  const double EPS_CUTOFF = 2e-2;
 double emax, nmax, emin, nmin;
 int ne, nn;
 double **ptab, **Ttab, **mubtab, **mustab, **stab;

public:
 EoSCMFeaux(char* filename, int Ne, int Nn);
 ~EoSCMFeaux();
 double read_file(ifstream& stream);
 void get(double e, double nb, double& p, double& T, double& mub, double& mus);
 double p(double e, double nb);
};

double EoSCMFeaux::read_file(ifstream& stream){
  string text;
  stream >> text;
  double d;
  if(text == "NAN") {
    d = numeric_limits<double>::quiet_NaN();
  }else{
    d = atof(text.c_str());
  }
 
  return d;

}
EoSCMFeaux::EoSCMFeaux(char* filename, int Ne, int Nn) {
 ne = Ne;
 nn = Nn;
 ptab = new double*[ne];
 Ttab = new double*[ne];
 mubtab = new double*[ne];
 mustab = new double*[ne];
 stab = new double*[ne];
 for (int i = 0; i < ne; i++) {
  ptab[i] = new double[nn];
  Ttab[i] = new double[nn];
  mubtab[i] = new double[nn];
  mustab[i] = new double[nn];
  stab[i] = new double[nn];
 }
 double* etab = new double[ne];
 double* ntab = new double[nn];
 ifstream fin(filename);
 if (!fin.good()) {
  cout << "I/O error with " << filename << endl;
  exit(1);
 }
 fin.ignore(1000, '\n');
 for (int ie = 0; ie < ne; ie++)
  for (int in = 0; in < nn; in++) {
   etab[ie]=read_file(fin)*EPS_0/1000.0;    // --> e[GeV/fm3]
   ntab[in]=read_file(fin)*N_0;
   Ttab[ie][in]=read_file(fin)/1000.0; // -> T[GeV]
   ptab[ie][in]=read_file(fin)* EPS_0/1000.0;     // --> p[GeV/fm3]
   stab[ie][in]=read_file(fin)* N_0;      // --> s[1/fm3]
   mubtab[ie][in]=read_file(fin)/ 1000.0;  // --> mub[GeV]
   mustab[ie][in]=read_file(fin)/ 1000.0;  // --> mus[GeV]

  }
 emin = etab[0];
 emax = etab[ne - 1];
 nmin = ntab[0];
 nmax = ntab[nn - 1];
 cout << "EoSCMFeaux: table " << filename
      << " read, [emin,emax,nmin,nmax] = " << emin << "  " << emax << "  "
      << nmin << "  " << nmax << endl;
 delete[] etab;
 delete[] ntab;
}

EoSCMFeaux::~EoSCMFeaux() {
 for (int i = 0; i < ne; i++) {
  delete[] ptab[i];
  delete[] Ttab[i];
  delete[] mubtab[i];
  delete[] mustab[i];
  delete[] stab[i];
 }
 delete ptab;
 delete Ttab;
 delete mubtab;
 delete mustab;
 delete stab;
}

void EoSCMFeaux::get(double e, double nb, double& p, double& T, double& mub,
                 double& mus) {
 if (e < EPS_CUTOFF) {
  T = mub = mus = p = 0.;
  return;
 }
 const double de = (emax - emin) / (ne - 1);
 const double dn = (nmax - nmin) / (nn - 1);
 int ie = (int)((e - emin) / de);
 int in = (int)((nb - nmin) / dn);
 if (ie < 0) ie = 0;
 if (in < 0) in = 0;
 if (ie > ne - 2) ie = ne - 2;
 if (in > nn - 2) in = nn - 2;
 const double em = e - emin - ie * de;
 const double nm = nb - nmin - in * dn;

 double we[2] = {1. - em / de, em / de};
 double wn[2] = {1. - nm / dn, nm / dn};

 T = mub = mus = p = 0.0;
 for (int je = 0; je < 2; je++)
  for (int jn = 0; jn < 2; jn++) {
   p += we[je] * wn[jn] * ptab[ie + je][in + jn];
   T += we[je] * wn[jn] * Ttab[ie + je][in + jn];
   mub += we[je] * wn[jn] * mubtab[ie + je][in + jn];
   mus += we[je] * wn[jn] * mustab[ie + je][in + jn];
  }
 if(isnan(p)) {
   if(e> 0.8*EPS_0  || nb> 40.*N_0){
    p = 0.2964 * e;
    T = 0.15120476935 * pow(e, 0.25);
    mub = mus = 0.0;
   }else{
    T = mub = mus = p = 0.;
   }
   
 }
 
 if (p < 0.0) p = 0.0;
}

double EoSCMFeaux::p(double e, double nb) {
 if (e < EPS_CUTOFF) return 0.0;
 const double de = (emax - emin) / (ne - 1);
 const double dn = (nmax - nmin) / (nn - 1);
 int ie = (int)((e - emin) / de);
 int in = (int)((nb - nmin) / dn);
 if (ie < 0) ie = 0;
 if (in < 0) in = 0;
 if (ie > ne - 2) ie = ne - 2;
 if (in > nn - 2) in = nn - 2;
 const double em = e - emin - ie * de;
 const double nm = nb - nmin - in * dn;

 double we[2] = {1. - em / de, em / de};
 double wn[2] = {1. - nm / dn, nm / dn};

 double p = 0.0;
 for (int je = 0; je < 2; je++)
  for (int jn = 0; jn < 2; jn++) p += we[je] * wn[jn] * ptab[ie + je][in + jn];

 if(isnan(p)) {
   if(e> 0.8*EPS_0  || nb> 40.*N_0){
    p = 0.2964 * e;
   }else{
    p = 0.;
   }
 }
 if (p < 0.0) p = 0.0;
 return p;
}

EoSCMFe::EoSCMFe() {
 eosbig = new EoSCMFeaux("eos/cmfe.dat", 4001, 401);
 eossmall = new EoSCMFeaux("eos/cmfe_small.dat", 1501, 501);
}

EoSCMFe::~EoSCMFe() {
 delete eosbig;
 delete eossmall;
}

void EoSCMFe::eos(double e, double nb, double nq, double ns, double& T,
                    double& mub, double& muq, double& mus, double& p) {
 if (e < 15*EPS_0/1000.0 && nb < 5*N_0)
  eossmall->get(e, nb, p, T, mub, mus);
 else if (e < 0.8*EPS_0 && nb < 40.*N_0)
  eosbig->get(e, nb, p, T, mub, mus);
 else {
  p = 0.2964 * e;
  T = 0.15120476935 * pow(e, 0.25);
  mub = mus = 0.0;
 }

 muq = 0.0;  
}

double EoSCMFe::p(double e, double nb, double nq, double ns) {
 if (e < 15*EPS_0/1000.0 && nb < 5*N_0)
  return eossmall->p(e, nb);
 else if (e < 0.8*EPS_0 && nb < 40.*N_0)
  return eosbig->p(e, nb);
 else
  return 0.2964 * e;
}
