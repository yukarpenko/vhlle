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
#include "eoSmash.h"

using namespace std;

EoSSmash::EoSSmash(char* filename, int Ne, int Nn) {
 ne = Ne;
 nn = Nn;
 ptab = new double*[ne];
 Ttab = new double*[ne];
 mubtab = new double*[ne];
 mustab = new double*[ne];
 for (int i = 0; i < ne; i++) {
  ptab[i] = new double[nn];
  Ttab[i] = new double[nn];
  mubtab[i] = new double[nn];
  mustab[i] = new double[nn];
 }
 double* e = new double[ne];
 double* n = new double[nn];
 ifstream fin(filename);
 if (!fin.good()) {
  cout << "I/O error with " << filename << endl;
  exit(1);
 }
 string header_str;
 getline(fin, header_str); // skip the header line
 for (int ie = 0; ie < ne; ie++)
  for (int in = 0; in < nn; in++) {
   // T [GeV], mu [GeV], energy density [GeV/fm3], baryon density [1/fm3]
   fin >> e[ie] >> n[in] >> Ttab[ie][in] >> mubtab[ie][in] >> mustab[ie][in];
   ptab[ie][in] = 0.0;
  }
 emin = e[0];
 emax = e[ne - 1];
 nmin = n[0];
 nmax = n[nn - 1];
 cout << "EoSaux: table " << filename
      << " read, [emin,emax,nmin,nmax] = " << emin << "  " << emax << "  "
      << nmin << "  " << nmax << endl;
 delete[] e;
 delete[] n;
}

EoSSmash::~EoSSmash() {
 for (int i = 0; i < ne; i++) {
  delete[] ptab[i];
  delete[] Ttab[i];
  delete[] mubtab[i];
  delete[] mustab[i];
 }
 delete ptab;
 delete Ttab;
 delete mubtab;
 delete mustab;
}


void EoSSmash::eos(double e, double nb, double nq, double ns, double& T,
                    double& mub, double& muq, double& mus, double& p) {
 if (e < 0.) {
  T = mub = muq = mus = p = 0.;
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
  muq = 0.;
 if (p < 0.0) p = 0.0;
}

double EoSSmash::p(double e, double nb, double nq, double ns) {
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
  for (int jn = 0; jn < 2; jn++) p += we[je] * wn[jn] * ptab[ie + je][in + jn];

 if (p < 0.0) p = 0.0;
 return p;
}
