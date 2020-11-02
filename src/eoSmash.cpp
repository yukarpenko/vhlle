#include <iostream>
#include <fstream>

#include "eos.h"
#include "eoSmash.h"

using namespace std;

EoSSmash::EoSSmash(char* filename, int Ne, int Nnb, int Nq) {
  ifstream fin(filename);
  if (!fin.good()) {
   cout << "I/O error with " << filename << endl;
   exit(1);
  }
 // Set members
 ne = Ne;
 nnb = Nnb;
 nnq = Nq;

 // Allocate arrays to later contain tabularised values of T, p, muB, muQ, muS
 ptab = new double[ne * nnb * nnq];
 Ttab = new double[ne * nnb * nnq];
 mubtab = new double[ne * nnb * nnq];
 muqtab = new double[ne * nnb * nnq];
 mustab = new double[ne * nnb * nnq];

 double *e = new double[ne];
 double *nb = new double[nnb];
 double *nq = new double[nnq];

 // Skip header
 string header_str;
 getline(fin, header_str);

 // Read table from file
 for (int ie = 0; ie < ne; ie++)
  for (int inb = 0; inb < nnb; inb++)
   for (int inq = 0; inq < nnq; inq++) {
    // e nb nq T p mub mus muq
    fin >> e[ie] >> nb[inb] >> nq[inq] >> Ttab[index3(ie, inb, inq)] >>
        ptab[index3(ie, inb, inq)] >> mubtab[index3(ie, inb, inq)] >>
        mustab[index3(ie, inb, inq)] >> muqtab[index3(ie, inb, inq)];
}

// Find upper and lower bounds of e, nB, nQ
emin = e[0];
emax = e[ne - 1];
nbmin = nb[0];
nbmax = nb[nnb - 1];
nqmin = nq[0];
nqmax = nq[nnq - 1];

cout << "EoSSMASH: table " << filename
      << " read, [emin,emax,nbmin,nbmax,qmin,qmax] = " << emin << "  " << emax << "  "
      << nbmin << "  " << nbmax << " " << nqmin << "  " << nqmax << endl;

// e, nb and nq are not needed anymore
delete[] e;
delete[] nb;
delete[] nq;
}

EoSSmash::~EoSSmash() {
 delete[] ptab;
 delete[] Ttab;
 delete[] mubtab;
 delete[] muqtab;
 delete[] mustab;
}


void EoSSmash::eos(double e, double nb, double nq, double ns, double &T,
                    double &mub, double &muq, double &mus, double &p) {

 // e < 0 is physically impossible
 if (e <= 0.) {
  T = mub = muq = mus = p = 0.;
  return;
 }

 // Linearly scale down lowest tabularised value
 if (e < emin) {
  p = e / emin * ptab[index3(0, nnb / 2, nnq / 2)];
  T = e / emin * Ttab[index3(0, nnb / 2, nnq / 2)];
  mub = muq = mus = 0.;
  return;
 }

 // Find width of steps in e, nB, nQ direction
 const double de = (emax - emin) / (ne - 1);
 const double dnb = (nbmax - nbmin) / (nnb - 1);
 const double dnq = (nqmax - nqmin) / (nnq - 1);

 // Find appropriate indices
 int ie = (int)((e - emin) / de);
 int inb = (int)((nb - nbmin) / dnb);
 int inq = (int)((nq - nqmin) / dnq);
 if (ie < 0) ie = 0;
 if (inb < 0) inb = 0;
 if (inq < 0) inq = 0;
 if (ie > ne - 2) ie = ne - 2;
 if (inb > nnb - 2) inb = nnb - 2;
 if (inq > nnq - 2) inq = nnq - 2;

 // Find fraction of step in e, nB and nQ direction
 const double em = e - emin - ie * de;
 const double nbm = nb - nbmin - inb * dnb;
 const double nqm = nq - nqmin - inq * dnq;

 // Assign weights for interpolation
 double we[2] = {1. - em / de, em / de};
 double wnb[2] = {1. - nbm / dnb, nbm / dnb};
 double wnq[2] = {1. - nqm / dnq, nqm / dnq};

 // Compute T, p, muB, muQ, muS by interpolating between neighbouring
 // tabularised points
 T = mub = muq = mus = p = 0.0;
 for (int je = 0; je < 2; je++)
  for (int jnb = 0; jnb < 2; jnb++)
   for (int jnq = 0; jnq < 2; jnq++) {
    p += we[je] * wnb[jnb] * wnq[jnq] *
         ptab[index3(ie + je, inb + jnb, inq + jnq)];
    T += we[je] * wnb[jnb] * wnq[jnq] *
         Ttab[index3(ie + je, inb + jnb, inq + jnq)];
    mub += we[je] * wnb[jnb] * wnq[jnq] *
           mubtab[index3(ie + je, inb + jnb, inq + jnq)];
    muq += we[je] * wnb[jnb] * wnq[jnq] *
           muqtab[index3(ie + je, inb + jnb, inq + jnq)];
    mus += we[je] * wnb[jnb] * wnq[jnq] *
           mustab[index3(ie + je, inb + jnb, inq + jnq)];
   }
 if (p < 0.0) p = 0.0;
}

double EoSSmash::p(double e, double nb, double nq, double ns) {
  if (e <= 0.) return 0.0;
  if (e < emin) return e / emin * ptab[index3(0, nnb / 2, nnq / 2)];

  // Find width of steps in e, nB, nQ direction
  const double de = (emax - emin) / (ne - 1);
  const double dnb = (nbmax - nbmin) / (nnb - 1);
  const double dnq = (nqmax - nqmin) / (nnq - 1);

  // Find appropriate indices
  int ie = (int)((e - emin) / de);
  int inb = (int)((nb - nbmin) / dnb);
  int inq = (int)((nq - nqmin) / dnq);
  if (ie < 0) ie = 0;
  if (inb < 0) inb = 0;
  if (inq < 0) inq = 0;
  if (ie > ne - 2) ie = ne - 2;
  if (inb > nnb - 2) inb = nnb - 2;
  if (inq > nnq - 2) inq = nnq - 2;

  // Find fraction of step in e, nB and nQ direction
  const double em = e - emin - ie * de;
  const double nbm = nb - nbmin - inb * dnb;
  const double nqm = nq - nqmin - inq * dnq;

  // Assign weights for interpolation
  double we[2] = {1. - em / de, em / de};
  double wnb[2] = {1. - nbm / dnb, nbm / dnb};
  double wnq[2] = {1. - nqm / dnq, nqm / dnq};

  // Compute p
  double p = 0.0;
  for (int je = 0; je < 2; je++)
   for (int jnb = 0; jnb < 2; jnb++)
    for (int jnq = 0; jnq < 2; jnq++) {
     p += we[je] * wnb[jnb] * wnq[jnq] *
          ptab[index3(ie + je, inb + jnb, inq + jnq)];
    }

  if (p < 0.0) p = 0.0;

  return p;
}
