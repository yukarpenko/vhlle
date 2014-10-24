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
#include "eo1.h"

using namespace std;

EoS1f::EoS1f(char *filename) {
  ifstream fin(filename);

  fin >> emax >> e0 >> nmax >> n0 >> ne >> nn;

  egrid = new double[ne];
  ngrid = new double[nn];
  lgegrid = new double[ne];
  lgngrid = new double[nn];

  for (int ixe = 0; ixe < ne; ixe++) fin >> egrid[ixe];

  for (int ixnb = 0; ixnb < nn; ixnb++) fin >> ngrid[ixnb];

  T = new double[ne * nn];
  pre = new double[ne * nn];
  mub = new double[ne * nn];

  for (int ie = 0; ie < ne; ie++)
    for (int inb = 0; inb < nn; inb++) {
      fin >> pre[index(ie, inb)] >> T[index(ie, inb)] >> mub[index(ie, inb)];
      if (pre[index(ie, inb)] < 1e-18) pre[index(ie, inb)] = 0.;
    }
  fin.close();
}

EoS1f::~EoS1f() {
  delete[] T;
  delete[] pre;
  delete[] mub;
}

void EoS1f::eosranges(double &_emax, double &_e0, double &_nmax, double &_n0,
                      int &_ne, int &_nn) {
  _emax = emax;
  _e0 = e0;
  _nmax = nmax;
  _n0 = n0;
  _ne = ne;
  _nn = nn;
}

void EoS1f::getue(double e, int &ixe, double &ue) {
  int k, k1, k2;
  k1 = 0;
  k2 = ne - 1;
  for (; k2 - k1 > 1;) {
    k = (k1 + k2) / 2;
    if (e < egrid[k])
      k2 = k;
    else
      k1 = k;
  }
  ixe = k1;
  if (ixe > ne - 2) ixe = ne - 2;
  double dxe = egrid[k2] - egrid[k1];
  double xem = e - egrid[k1];
  ue = xem / dxe;
}

void EoS1f::getun(double n, int &ixn, double &un) {
  int k, k1, k2;
  k1 = 0;
  k2 = nn - 1;
  for (; k2 - k1 > 1;) {
    k = (k1 + k2) / 2;
    if (n < ngrid[k])
      k2 = k;
    else
      k1 = k;
  }
  ixn = k1;
  if (ixn > nn - 2) ixn = nn - 2;
  if (ixn < 0) ixn = 0;
  double dxn = ngrid[k2] - ngrid[k1];
  double xnm = n - ngrid[k1];
  un = xnm / dxn;
}

void EoS1f::eos(double e, double nb, double nq, double ns, double &_T,
                double &_mub, double &_muq, double &_mus, double &_p) {
  if (e < 0.) {
    _T = _mub = _muq = _mus = _p = 0.;
  }
  int ixe, ixnb;
  double ue, ub;

  getue(e, ixe, ue);
  getun(nb, ixnb, ub);

  double we[2] = {1. - ue, ue};
  double wnb[2] = {1. - ub, ub};

  _T = _mub = _muq = _mus = _p = 0.;
  for (int je = 0; je < 2; je++)
    for (int jnb = 0; jnb < 2; jnb++) {
      _p += we[je] * wnb[jnb] * pre[index(ixe + je, ixnb + jnb)];
      _T += we[je] * wnb[jnb] * T[index(ixe + je, ixnb + jnb)];
      _mub += we[je] * wnb[jnb] * mub[index(ixe + je, ixnb + jnb)];
    }
  if (_p < 0.) _p = 0.;
  // cout <<  e <<" "<< nb <<" "<< nq <<" "<< ns <<" "<< _T<<endl;
}
