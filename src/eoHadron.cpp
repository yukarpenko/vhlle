#include <cmath>
#include <iomanip>
#include <ctime>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include "eos.h"
#include "eoHadron.h"

using namespace std;

EoSHadron::EoSHadron(char *filename) {
  ifstream fin(filename);
  if (!fin.good()) {
    cout << "I/O error with " << filename << endl;
    exit(1);
  }
  fin >> ne >> nnb >> nnq >> e0 >> n0 >> logemin >> logemax >> lognmax;

  ptab = new double[ne * nnb * nnq];
  Ttab = new double[ne * nnb * nnq];
  mubtab = new double[ne * nnb * nnq];
  muqtab = new double[ne * nnb * nnq];
  mustab = new double[ne * nnb * nnq];
  statustab = new int[ne * nnb * nnq];

  double *e = new double[ne];
  double *nb = new double[nnb];
  double *nq = new double[nnq];

  for (int ie = 0; ie < ne; ie++)
    for (int inb = 0; inb < nnb; inb++)
      for (int inq = 0; inq < nnq; inq++) {
        // e  nb  p  T  mub  muq  mus  status
        fin >> e[ie] >> nb[inb] >> nq[inq] >> ptab[index3(ie, inb, inq)] >>
            Ttab[index3(ie, inb, inq)] >> mubtab[index3(ie, inb, inq)] >>
            muqtab[index3(ie, inb, inq)] >> mustab[index3(ie, inb, inq)] >>
            statustab[index3(ie, inb, inq)];
      }
  double emin = e[0];
  double emax = e[ne - 1];
  double nbmin = nb[0];
  double nbmax = nb[nnb - 1];
  double nqmin = nq[0];
  double nqmax = nq[nnq - 1];

  const double dxnb = (2.0 * lognmax) / (nnb - 1);
  const double dxnq = (2.0 * lognmax) / (nnq - 1);
  nb_abs_min = n0 * exp(2.0 * dxnb - lognmax);
  nq_abs_min = n0 * exp(2.0 * dxnq - lognmax);
  e_min = e0 * exp(logemin);

  if (fabs((e0 * exp(logemin) - emin) / emin) > 1e-5 ||
      fabs((e0 * exp(logemax) - emax) / emax) > 1e-5 ||
      fabs(n0 * exp(lognmax) - nbmax) > 1e-4) {
    cout << "EoSHadron: wrong eps or nb range: " << setw(14) << emin << setw(14)
         << emax << setw(14) << nbmin << setw(14) << nbmax << endl;
    exit(1);
  }
  cout << "EoHadron: table " << filename
       << " read, [emin,emax,nmin,nmax] = " << emin << "  " << emax << "  "
       << nbmin << "  " << nbmax << "  " << nqmin << "  " << nqmax << endl;
  delete[] e;
  delete[] nb;
  delete[] nq;
}

EoSHadron::~EoSHadron() {
  delete[] ptab;
  delete[] Ttab;
  delete[] mubtab;
  delete[] muqtab;
  delete[] mustab;
  delete[] statustab;
}

void EoSHadron::eos(double e, double nb, double nq, double ns, double &T,
                    double &mub, double &muq, double &mus, double &p, double tau) {
  if (e <= 0.) {
    T = mub = muq = mus = p = 0.;
    return;
  }
  if (e < e_min) {
    p = e / e_min * ptab[index3(0, nnb / 2, nnq / 2)];
    T = e / e_min * Ttab[index3(0, nnb / 2, nnq / 2)];
    mub = muq = mus = 0.;
    return;
  }
  const double xe = log(e / e0);
  const double dxe = (logemax - logemin) / (ne - 1);
  const double dxnb = (2.0 * lognmax) / (nnb - 1);
  const double dxnq = (2.0 * lognmax) / (nnq - 1);
  double xnb, xnq;  // position of the point in the table
  // some white magic with the inverse log mapping function for nb, nq
  if (nb != 0.0) {
    const double znb = log(fabs(nb) / n0);
    if (znb >= -lognmax + 2.0 * dxnb)
      xnb = nb / fabs(nb) * 0.5 * (znb + lognmax);
    else
      xnb = nb / nb_abs_min * dxnb;
  } else
    xnb = 0.;
  if (nq != 0.0) {
    const double znq = log(fabs(nq) / n0);
    if (znq >= -lognmax + 2.0 * dxnq)
      xnq = nq / fabs(nq) * 0.5 * (znq + lognmax);
    else
      xnq = nq / nq_abs_min * dxnq;
  } else
    xnq = 0.;

  int ie = (int)((xe - logemin) / dxe);
  int inb = (int)((xnb + lognmax) / dxnb);
  int inq = (int)((xnq + lognmax) / dxnq);
  if (ie < 0) ie = 0;
  if (inb < 0) inb = 0;
  if (inq < 0) inq = 0;
  if (ie > ne - 2) ie = ne - 2;
  if (inb > nnb - 2) inb = nnb - 2;
  if (inq > nnq - 2) inq = nnq - 2;
  const double em = xe - logemin - ie * dxe;
  const double nbm = xnb + lognmax - inb * dxnb;
  const double nqm = xnq + lognmax - inq * dxnq;

  if (statustab[index3(ie, inb, inq)] == 1) {
    T = mub = muq = mus = p = 0.;
    return;
  }

  double we[2] = {1. - em / dxe, em / dxe};
  double wnb[2] = {1. - nbm / dxnb, nbm / dxnb};
  double wnq[2] = {1. - nqm / dxnq, nqm / dxnq};

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
  // cout <<  e <<" "<< nb <<" "<< nq <<" "<< ns <<" "<< _T<<endl;
}

double EoSHadron::p(double e, double nb, double nq, double ns) {
  if (e <= 0.) return 0.0;
  if (e < e_min) return e / e_min * ptab[index3(0, nnb / 2, nnq / 2)];
  const double xe = log(e / e0);
  const double dxe = (logemax - logemin) / (ne - 1);
  const double dxnb = (2.0 * lognmax) / (nnb - 1);
  const double dxnq = (2.0 * lognmax) / (nnq - 1);
  double xnb, xnq;  // position of the point in the table
  // some white magic with the inverse log mapping function for nb, nq
  if (nb != 0.0) {
    const double znb = log(fabs(nb) / n0);
    if (znb >= -lognmax + 2.0 * dxnb)
      xnb = nb / fabs(nb) * 0.5 * (znb + lognmax);
    else
      xnb = nb / nb_abs_min * dxnb;
  } else
    xnb = 0.;
  if (nq != 0.0) {
    const double znq = log(fabs(nq) / n0);
    if (znq >= -lognmax + 2.0 * dxnq)
      xnq = nq / fabs(nq) * 0.5 * (znq + lognmax);
    else
      xnq = nq / nq_abs_min * dxnq;
  } else
    xnq = 0.;

  int ie = (int)((xe - logemin) / dxe);
  int inb = (int)((xnb + lognmax) / dxnb);
  int inq = (int)((xnq + lognmax) / dxnq);
  if (ie < 0) ie = 0;
  if (inb < 0) inb = 0;
  if (inq < 0) inq = 0;
  if (ie > ne - 2) ie = ne - 2;
  if (inb > nnb - 2) inb = nnb - 2;
  if (inq > nnq - 2) inq = nnq - 2;
  const double em = xe - logemin - ie * dxe;
  const double nbm = xnb + lognmax - inb * dxnb;
  const double nqm = xnq + lognmax - inq * dxnq;

  if (statustab[index3(ie, inb, inq)] == 1) return 0.0;

  double we[2] = {1. - em / dxe, em / dxe};
  double wnb[2] = {1. - nbm / dxnb, nbm / dxnb};
  double wnq[2] = {1. - nqm / dxnq, nqm / dxnq};

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
