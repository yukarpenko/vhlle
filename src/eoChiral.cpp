//#include <TError.h>
//#include <TApplication.h>
//#include <TGraph.h>
//#include <TCanvas.h>
//#include <TMath.h>
//#include <TGraph.h>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <cstdlib>

#include <iostream>
#include <fstream>
//#include <TF1.h>

#include "eos.h"
#include "eoChiral.h"

using namespace std;

// ---- auxiliary EoS class. Two instances (objects) of this class will be
// created
// ---- to store two EoS tables: chiraleos and chiralsmall
// ---- Not performance-wise, but the code is elegant :)
class EoSaux {
  double emax, nmax, emin, nmin;
  int ne, nn;
  double** ptab, **Ttab, **mubtab, **mustab, **stab;

 public:
  EoSaux(char* filename, int Ne, int Nn);
  ~EoSaux();
  void get(double e, double nb, double& p, double& T, double& mub, double& mus);
  double p(double e, double nb);
};

EoSaux::EoSaux(char* filename, int Ne, int Nn) {
  ne = Ne;
  nn = Nn;
  ptab = new double* [ne];
  Ttab = new double* [ne];
  mubtab = new double* [ne];
  mustab = new double* [ne];
  stab = new double* [ne];
  for (int i = 0; i < ne; i++) {
    ptab[i] = new double[nn];
    Ttab[i] = new double[nn];
    mubtab[i] = new double[nn];
    mustab[i] = new double[nn];
    stab[i] = new double[nn];
  }
  double* e = new double[ne];
  double* n = new double[nn];
  ifstream fin(filename);
  double a;  // dummy variable
  if (!fin.good()) {
    cout << "I/O error with " << filename << endl;
    exit(1);
  }
  for (int in = 0; in < nn; in++)
    for (int ie = 0; ie < ne; ie++) {
      // T [MeV], mu_q [MeV], energy density [e_0], pressure [e_0], baryon
      // density [n_0], entropy density [n_0], mu_S [MeV], placeholder (not
      // important).
      fin >> Ttab[ie][in] >> mubtab[ie][in] >> e[ie] >> ptab[ie][in] >> n[in] >>
          stab[ie][in] >> mustab[ie][in] >> a;
      Ttab[ie][in] /= 1000.0;    // --> T[GeV]
      mubtab[ie][in] /= 1000.0;  // --> mub[GeV]
      mustab[ie][in] /= 1000.0;  // --> mus[GeV]
      ptab[ie][in] *= 0.146;     // --> p[GeV/fm3]
      stab[ie][in] *= 0.15;      // --> s[1/fm3]
    }
  emin = e[0] * 0.146;
  emax = e[ne - 1] * 0.146;
  nmin = n[0] * 0.15;
  nmax = n[nn - 1] * 0.15;
  cout << "EoSaux: table " << filename
       << " read, [emin,emax,nmin,nmax] = " << emin << "  " << emax << "  "
       << nmin << "  " << nmax << endl;
  delete[] e;
  delete[] n;
}

EoSaux::~EoSaux() {
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

void EoSaux::get(double e, double nb, double& p, double& T, double& mub,
                 double& mus) {
  if (e < 0.) {
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
  if (p < 0.0) p = 0.0;
  // cout <<  e <<" "<< nb <<" "<< nq <<" "<< ns <<" "<< _T<<endl;
}

double EoSaux::p(double e, double nb) {
  if (e < 0.) return 0.0;
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
    for (int jn = 0; jn < 2; jn++)
      p += we[je] * wn[jn] * ptab[ie + je][in + jn];

  if (p < 0.0) p = 0.0;
  return p;
  // cout <<  e <<" "<< nb <<" "<< nq <<" "<< ns <<" "<< _T<<endl;
}

EoSChiral::EoSChiral() {
  eosbig = new EoSaux("eos/chiraleos.dat", 2001, 401);
  eossmall = new EoSaux("eos/chiralsmall.dat", 201, 201);
}

EoSChiral::~EoSChiral() {
  delete eosbig;
  delete eossmall;
}

void EoSChiral::eos(double e, double nb, double nq, double ns, double& T,
                    double& mub, double& muq, double& mus, double& p, double tau) {
  if (e < 1.46 && nb < 0.3)
    eossmall->get(e, nb, p, T, mub, mus);
  else if (e < 146. && nb < 6.)
    eosbig->get(e, nb, p, T, mub, mus);
  else {
    p = 0.2964 * e;
    T = 0.15120476935 * pow(e, 0.25);
    mub = mus = 0.0;
  }
  muq = 0.0;  // generally it's not zero, but...but
}

double EoSChiral::p(double e, double nb, double nq, double ns, double tau) {
  if (e < 1.46 && nb < 0.3)
    return eossmall->p(e, nb);
  else if (e < 146. && nb < 6.)
    return eosbig->p(e, nb);
  else
    return 0.2964 * e;
}
