#include <math.h>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "eoAZH.h"
#include "eos.h"

using namespace std;

// ---- auxiliary EoS class. Six instances (objects) of this class will be
// created
// ---- to store six EoS tables: p, T, mu in two energy/density regions
class EoSAZHaux {
 double dn, de, emin, nmin;
 int ne, nn;
 double **tab;

public:
 EoSAZHaux(const char *filename);
 ~EoSAZHaux();
 inline double emin_() { return emin; }
 inline double emax() { return emin + de * (ne - 1); }
 inline double nmax() { return nmin + dn * (nn - 1); }
 double get(double e, double nb);
 double getLow(double e, double nb);
};

EoSAZHaux::EoSAZHaux(const char *filename) {
 ifstream fin(filename);
 if (!fin.good()) {
  cout << "I/O error with " << filename << endl;
  exit(1);
 }
 fin >> nmin >> emin;
 fin >> dn >> de >> nn >> ne;
 nn++;
 ne++;
 tab = new double *[ne];
 for (int i = 0; i < ne; i++) {
  tab[i] = new double[nn];
 }
 for (int ie = 0; ie < ne; ie++)
  for (int in = 0; in < nn; in++) {
   fin >> tab[ie][in];
  }
 for (int ie = 0; ie < 15; ie++) cout << "[" << tab[ie][0] << "]";
 cout << endl;
 cout << "EoSaux: table " << filename << " read, [ne nb] = " << ne << "  " << nn
      << "  done\n";
}

EoSAZHaux::~EoSAZHaux() {
 for (int i = 0; i < ne; i++) delete[] tab[i];
 delete tab;
}

double EoSAZHaux::get(double e, double nb) {
 if (e <= 0.) {
  return 0.0;
 }
 int ie = (int)((e - emin) / de);
 int in = (int)((nb - nmin) / dn);
 if (ie < 0.) ie = 0;
 if (in < 0.) in = 0;
 if (ie > ne - 2) ie = ne - 2;
 if (in > nn - 2) in = nn - 2;
 const double em = (e - emin - ie * de) / de;
 const double nm = (nb - nmin - in * dn) / dn;

 double we[2] = {1. - em, em};
 double wn[2] = {1. - nm, nm};

 double value = 0.0;
 for (int je = 0; je < 2; je++)
  for (int jn = 0; jn < 2; jn++)
   value += we[je] * wn[jn] * tab[ie + je][in + jn];
 //	if(value<0.0) value = 0.0 ;
 if (value != value)
  cout << "EoSAZH: " << e << "  " << nb << "  " << ie << "  " << in << endl;
 return value;
}

double EoSAZHaux::getLow(double e, double nb) {
 if (e <= 0.) {
  return 0.0;
 }
 if ((nb - nmin) / dn > 10000) {
  cout << "getLow: " << e << "  " << nb * e / emin << endl;
  exit(1);
 }
 int in = (int)((nb - nmin) / dn);
 if (in < 0.) in = 0;
 if (in > nn - 2) in = nn - 2;
 const double nm = (nb - nmin - in * dn) / dn;
 double wn[2] = {1. - nm, nm};
 double value = wn[0] * tab[0][in] + wn[1] * tab[0][in + 1];
 if (value != value)
  cout << "EoSAZH,getLow: " << e << "  " << nb << "  " << in << endl;
 return value * e / emin;
}

EoSAZH::EoSAZH() {
 p1 = new EoSAZHaux("eos/azhydro0p2/aa1_p.dat");
 p2 = new EoSAZHaux("eos/azhydro0p2/aa2_p.dat");
 T1 = new EoSAZHaux("eos/azhydro0p2/aa1_t.dat");
 T2 = new EoSAZHaux("eos/azhydro0p2/aa2_t.dat");
 mu1 = new EoSAZHaux("eos/azhydro0p2/aa1_mb.dat");
 mu2 = new EoSAZHaux("eos/azhydro0p2/aa2_mb.dat");
}

EoSAZH::~EoSAZH() {
 delete p1;
 delete p2;
 delete T1;
 delete T2;
 delete mu1;
 delete mu2;
}

void EoSAZH::eos(double e, double nb, double nq, double ns, double &T,
                 double &mub, double &muq, double &mus, double &p) {
 if (e > p1->emax()) {
  p = p2->get(e, fabs(nb));
  T = T2->get(e, fabs(nb));
  mub = mu2->get(e, nb);
 } else {
  if (e > p1->emin_()) {
   p = p1->get(e, fabs(nb));
   T = T1->get(e, fabs(nb));
   mub = mu1->get(e, nb);
  } else {
   p = p1->getLow(e, fabs(nb));
   T = T1->getLow(e, fabs(nb));
   mub = mu1->getLow(e, nb);
  }
 }
 muq = mus = 0.0;
}

double EoSAZH::p(double e, double nb, double nq, double ns) {
 if (e > p1->emax())
  return p2->get(e, fabs(nb));
 else if (e > p1->emin_())
  return p1->get(e, fabs(nb));
 else
  return p1->getLow(e, fabs(nb));
}
