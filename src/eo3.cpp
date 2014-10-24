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
#include "eo3.h"

using namespace std;

double sign(double x) {
  if (x > 0)
    return 1.;
  else if (x < 0.)
    return -1.;
  else
    return 0.;
}

EoS3f::EoS3f(char *filename, double _B, double _volex0, double _delta0,
             double _aaa, double _bbb) {
  ifstream fin(filename);

  fin >> B >> volex0 >> delta0 >> aaa >> bbb;
  if (B > 0) {
    if (fabs(B - _B) > 1e-5) {
      cout << "\n\n\n      B =      " << B << " instead of " << _B
           << "; exiting\n\n\n\n";
      exit(0);
    }
    if (fabs(volex0 - _volex0) > 1e-5) {
      cout << "\n\n\n      volex0 = " << volex0 << " instead of " << _volex0
           << "; exiting\n\n\n\n";
      exit(0);
    }
    if (fabs(delta0 - _delta0) > 1e-5) {
      cout << "\n\n\n      delta0 = " << delta0 << " instead of " << _delta0
           << "; exiting\n\n\n\n";
      exit(0);
    }
    if (fabs(aaa - _aaa) > 1e-5) {
      cout << "\n\n\n      aaa =    " << aaa << " instead of " << _aaa
           << "; exiting\n\n\n\n";
      exit(0);
    }
    if (fabs(bbb - _bbb) > 1e-5) {
      cout << "\n\n\n      bbb =    " << bbb << " instead of " << _bbb
           << "; exiting\n\n\n\n";
      exit(0);
    }
  }

  fin >> emax >> e0 >> nmax >> n0 >> ne >> nn;

  egrid = new double[ne];
  ngrid = new double[nn];
  lgegrid = new double[ne];
  lgngrid = new double[nn];

  for (int ixe = 0; ixe < ne; ixe++) fin >> egrid[ixe];

  for (int ixnb = 0; ixnb < nn; ixnb++) fin >> ngrid[ixnb];

  T = new double[ne * nn * nn * nn];
  pre = new double[ne * nn * nn * nn];
  mub = new double[ne * nn * nn * nn];
  muq = new double[ne * nn * nn * nn];
  mus = new double[ne * nn * nn * nn];

  for (int ie = 0; ie < ne; ie++)
    for (int inb = 0; inb < nn; inb++)
      for (int inq = 0; inq < nn; inq++)
        for (int ins = 0; ins < nn; ins++) {
          fin >> pre[index(ie, inb, inq, ins)] >> T[index(ie, inb, inq, ins)] >>
              mub[index(ie, inb, inq, ins)] >> muq[index(ie, inb, inq, ins)] >>
              mus[index(ie, inb, inq, ins)];
          if (pre[index(ie, inb, inq, ins)] < 1e-18)
            pre[index(ie, inb, inq, ins)] = 0.;
        }
  fin.close();
}

EoS3f::~EoS3f() {
  delete[] T;
  delete[] pre;
  delete[] mub;
  delete[] muq;
  delete[] mus;
}

void EoS3f::eosranges(double &_emax, double &_e0, double &_nmax, double &_n0,
                      int &_ne, int &_nn) {
  _emax = emax;
  _e0 = e0;
  _nmax = nmax;
  _n0 = n0;
  _ne = ne;
  _nn = nn;
}

void EoS3f::getue(double e, int &ixe, double &ue) {
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

void EoS3f::getun(double n, int &ixn, double &un) {
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

void EoS3f::eos(double e, double nb, double nq, double ns, double &_T,
                double &_mub, double &_muq, double &_mus, double &_p) {
  if (e < 0.) {
    _T = _mub = _muq = _mus = _p = 0.;
  }
  int ixe, ixnb, ixns, ixnq;
  double ue, ub, uq, us;

  getue(e, ixe, ue);
  getun(nb, ixnb, ub);
  getun(nq, ixnq, uq);
  getun(ns, ixns, us);

  double we[2] = {1. - ue, ue};
  double wnb[2] = {1. - ub, ub};
  double wnq[2] = {1. - uq, uq};
  double wns[2] = {1. - us, us};

  _T = _mub = _muq = _mus = _p = 0.;
  for (int je = 0; je < 2; je++)
    for (int jnb = 0; jnb < 2; jnb++)
      for (int jnq = 0; jnq < 2; jnq++)
        for (int jns = 0; jns < 2; jns++) {
          _p += we[je] * wnb[jnb] * wnq[jnq] * wns[jns] *
                pre[index(ixe + je, ixnb + jnb, ixnq + jnq, ixns + jns)];
          _T += we[je] * wnb[jnb] * wnq[jnq] * wns[jns] *
                T[index(ixe + je, ixnb + jnb, ixnq + jnq, ixns + jns)];
          _mub += we[je] * wnb[jnb] * wnq[jnq] * wns[jns] *
                  mub[index(ixe + je, ixnb + jnb, ixnq + jnq, ixns + jns)];
          _muq += we[je] * wnb[jnb] * wnq[jnq] * wns[jns] *
                  muq[index(ixe + je, ixnb + jnb, ixnq + jnq, ixns + jns)];
          _mus += we[je] * wnb[jnb] * wnq[jnq] * wns[jns] *
                  mus[index(ixe + je, ixnb + jnb, ixnq + jnq, ixns + jns)];
        }
  if (_p < 0.) _p = 0.;
  // cout <<  e <<" "<< nb <<" "<< nq <<" "<< ns <<" "<< _T<<endl;
}
