/******************************************************************************
*                                                                             *
*            vHLLE : a 3D viscous hydrodynamic code                           *
*            version 1.1,            October 2014                            *
*            by Iurii Karpenko                                                *
*  contact:  yu.karpenko@gmail.com                                            *
*  For the detailed description please refer to:                              *
*  http://arxiv.org/abs/1312.4160                                             *
*                                                                             *
*  This code can be freely used and redistributed, provided that this         *
*  copyright appear in all the copies. If you decide to make modifications    *
*  to the code, please contact the authors, especially if you plan to publish *
* the results obtained with such modified code. Any publication of results    *
* obtained using this code must include the reference to                      *
* arXiv:1312.4160 [nucl-th] or the published version of it, when available.   *
*                                                                             *
*******************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "inc.h"
#include "rmn.h"
#include "fld.h"
#include "cll.h"
#include "eos.h"
#include "trancoeff.h"
#include "cornelius.h"

// pi^{mu nu} and Pi output
//#define OUTPI

// change to hadron EoS (e.g. Laine) to calculate v,T,mu at the surface
#define SWAP_EOS

using namespace std;

// returns the velocities in cartesian coordinates, fireball rest frame.
// Y=longitudinal rapidity of fluid
void Fluid::getCMFvariables(Cell *c, double tau, double &e, double &nb,
                            double &nq, double &ns, double &vx, double &vy,
                            double &Y) {
  double p, vz;
  c->getPrimVar(eos, tau, e, p, nb, nq, ns, vx, vy, vz);
  double eta = getZ(c->getZ());
  //	Y = eta + TMath::ATanH(vz) ;
  Y = eta + 1. / 2. * log((1. + vz) / (1. - vz));
  vx = vx * cosh(Y - eta) / cosh(Y);
  vy = vy * cosh(Y - eta) / cosh(Y);
}

Fluid::Fluid(EoS *_eos, EoS *_eosH, TransportCoeff *_trcoeff, int _nx, int _ny,
             int _nz, double _minx, double _maxx, double _miny, double _maxy,
             double _minz, double _maxz, double _dt, double eCrit) {
  eos = _eos;
  eosH = _eosH;
  trcoeff = _trcoeff;
  nx = _nx;
  ny = _ny;
  nz = _nz;
  minx = _minx;
  maxx = _maxx;
  miny = _miny;
  maxy = _maxy;
  minz = _minz;
  maxz = _maxz;
  dx = (maxx - minx) / (nx - 1);
  dy = (maxy - miny) / (ny - 1);
  dz = (maxz - minz) / (nz - 1);
  dt = _dt;

  cell = new Cell[nx * ny * nz];

  cell0 = new Cell;
  cell0->setPrimVar(eos, 1.0, 0., 0., 0., 0., 0., 0.,
                    0.);  // tau is not important here, since *0
  cell0->setAllM(0.);

  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++)
      for (int iz = 0; iz < nz; iz++) {
        getCell(ix, iy, iz)->setPrev(X_, getCell(ix - 1, iy, iz));
        getCell(ix, iy, iz)->setNext(X_, getCell(ix + 1, iy, iz));
        getCell(ix, iy, iz)->setPrev(Y_, getCell(ix, iy - 1, iz));
        getCell(ix, iy, iz)->setNext(Y_, getCell(ix, iy + 1, iz));
        getCell(ix, iy, iz)->setPrev(Z_, getCell(ix, iy, iz - 1));
        getCell(ix, iy, iz)->setNext(Z_, getCell(ix, iy, iz + 1));
        getCell(ix, iy, iz)->setPos(ix, iy, iz);
      }

  output_nt = 0;
  output_nx = 0;
  output_ny = 0;

  //---- Cornelius init
  double arrayDx[4] = {dt, dx, dy, dz};
  cornelius = new Cornelius;
  cornelius->init(4, eCrit, arrayDx);
  ecrit = eCrit;
  vEff = 0.;
  EtotSurf = 0.0;
}

Fluid::~Fluid() {
  delete[] cell;
  delete cell0;
}

void Fluid::initOutput(char *dir, int maxstep, double tau0, int cmpr2dOut) {
  //    directory = dir ;
  compress2dOut = cmpr2dOut;
  cout << "maxstep=" << maxstep << endl;
  char command[255];
  sprintf(command, "mkdir -p %s", dir);
  int return_mkdir = system(command);
  cout << "mkdir returns: " << return_mkdir << endl;
  string outx = dir;
  outx.append("/outx.dat");
  string outxvisc = dir;
  outxvisc.append("/outx.visc.dat");
  string outyvisc = dir;
  outyvisc.append("/outy.visc.dat");
  string outdiagvisc = dir;
  outdiagvisc.append("/diag.visc.dat");
  string outy = dir;
  outy.append("/outy.dat");
  string outdiag = dir;
  outdiag.append("/outdiag.dat");
  string outz = dir;
  outz.append("/outz.dat");
  string outaniz = dir;
  outaniz.append("/out.aniz.dat");
  string out2d = dir;
  out2d.append("/out2D.dat");
  string out3d = dir;
  out3d.append("/3D.dat");
  string outfreeze = dir;
  outfreeze.append("/freezeout.dat");
  foutx.open(outx.c_str());
  fouty.open(outy.c_str());
  foutz.open(outz.c_str());
  foutdiag.open(outdiag.c_str());
  fout2d.open(out2d.c_str());
  fout3d.open(out3d.c_str());
  foutxvisc.open(outxvisc.c_str());
  foutyvisc.open(outyvisc.c_str());
  foutdiagvisc.open(outdiagvisc.c_str());
  fout_aniz.open(outaniz.c_str());
  ffreeze.open(outfreeze.c_str());
  //################################################################
  // important remark. for correct diagonal output, nx=ny must hold.
  //################################################################
  foutxvisc << maxstep + 1 << "  " << getNX() << endl;
  foutyvisc << maxstep + 1 << "  " << getNY() << endl;
  foutdiagvisc << maxstep + 1 << "  " << getNX() << endl;
  foutx << "# " << maxstep + 1 << "  " << getNX() << endl;
  fouty << "# " << maxstep + 1 << "  " << getNY() << endl;
  foutdiag << "# " << maxstep + 1 << "  " << getNX() << endl;
  foutz << "# " << maxstep + 1 << "  " << getNZ() << endl;
  fout2d << " " << maxstep + 1 << "  " << (getNX() - 5) + 1 << "  "
         << (getNY() - 5) + 1 << endl;
  fout2d << tau0 << "  " << tau0 + 0.05 * maxstep << "  " << getX(2) << "  "
         << getX(getNX() - 3) << "  " << getY(2) << "  " << getY(getNY() - 3)
         << endl;
  outputGnuplot(tau0);
  fout_aniz << "#  tau  <<v_T>>  e_p  e'_p  (to compare with SongHeinz)\n";
}

void Fluid::correctImagCells(void) {
  double Q[7];
  // Z
  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++) {
      // left boundary
      getCell(ix, iy, 2)->getQ(Q);
      getCell(ix, iy, 1)->setQ(Q);
      getCell(ix, iy, 0)->setQ(Q);
      // right boundary
      getCell(ix, iy, nz - 3)->getQ(Q);
      getCell(ix, iy, nz - 2)->setQ(Q);
      getCell(ix, iy, nz - 1)->setQ(Q);
    }
  // Y
  for (int ix = 0; ix < nx; ix++)
    for (int iz = 0; iz < nz; iz++) {
      // left boundary
      getCell(ix, 2, iz)->getQ(Q);
      getCell(ix, 1, iz)->setQ(Q);
      getCell(ix, 0, iz)->setQ(Q);
      // right boundary
      getCell(ix, ny - 3, iz)->getQ(Q);
      getCell(ix, ny - 2, iz)->setQ(Q);
      getCell(ix, ny - 1, iz)->setQ(Q);
    }
  // X
  for (int iy = 0; iy < ny; iy++)
    for (int iz = 0; iz < nz; iz++) {
      // left boundary
      getCell(2, iy, iz)->getQ(Q);
      getCell(1, iy, iz)->setQ(Q);
      getCell(0, iy, iz)->setQ(Q);
      // right boundary
      getCell(nx - 3, iy, iz)->getQ(Q);
      getCell(nx - 2, iy, iz)->setQ(Q);
      getCell(nx - 1, iy, iz)->setQ(Q);
    }
}

void Fluid::correctImagCellsFull(void) {
  double Q[7], _pi[4][4], _Pi;
  // Z
  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++) {
      // left boundary
      getCell(ix, iy, 2)->getQ(Q);
      getCell(ix, iy, 1)->setQ(Q);
      getCell(ix, iy, 0)->setQ(Q);
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) _pi[i][j] = getCell(ix, iy, 2)->getpi(i, j);
      _Pi = getCell(ix, iy, 2)->getPi();

      for (int i = 0; i < 4; i++)
        for (int j = 0; j <= i; j++) {
          getCell(ix, iy, 0)->setpi(i, j, _pi[i][j]);
          getCell(ix, iy, 1)->setpi(i, j, _pi[i][j]);
        }
      getCell(ix, iy, 0)->setPi(_Pi);
      getCell(ix, iy, 1)->setPi(_Pi);
      // right boundary
      getCell(ix, iy, nz - 3)->getQ(Q);
      getCell(ix, iy, nz - 2)->setQ(Q);
      getCell(ix, iy, nz - 1)->setQ(Q);
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
          _pi[i][j] = getCell(ix, iy, nz - 3)->getpi(i, j);
      _Pi = getCell(ix, iy, nz - 3)->getPi();

      for (int i = 0; i < 4; i++)
        for (int j = 0; j <= i; j++) {
          getCell(ix, iy, nz - 2)->setpi(i, j, _pi[i][j]);
          getCell(ix, iy, nz - 1)->setpi(i, j, _pi[i][j]);
        }
      getCell(ix, iy, nz - 2)->setPi(_Pi);
      getCell(ix, iy, nz - 1)->setPi(_Pi);
    }
  // Y
  for (int ix = 0; ix < nx; ix++)
    for (int iz = 0; iz < nz; iz++) {
      // left boundary
      getCell(ix, 2, iz)->getQ(Q);
      getCell(ix, 1, iz)->setQ(Q);
      getCell(ix, 0, iz)->setQ(Q);
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) _pi[i][j] = getCell(ix, 2, iz)->getpi(i, j);
      _Pi = getCell(ix, 2, iz)->getPi();

      for (int i = 0; i < 4; i++)
        for (int j = 0; j <= i; j++) {
          getCell(ix, 0, iz)->setpi(i, j, _pi[i][j]);
          getCell(ix, 1, iz)->setpi(i, j, _pi[i][j]);
        }
      getCell(ix, 0, iz)->setPi(_Pi);
      getCell(ix, 1, iz)->setPi(_Pi);
      // right boundary
      getCell(ix, ny - 3, iz)->getQ(Q);
      getCell(ix, ny - 2, iz)->setQ(Q);
      getCell(ix, ny - 1, iz)->setQ(Q);
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
          _pi[i][j] = getCell(ix, ny - 3, iz)->getpi(i, j);
      _Pi = getCell(ix, ny - 3, iz)->getPi();

      for (int i = 0; i < 4; i++)
        for (int j = 0; j <= i; j++) {
          getCell(ix, ny - 2, iz)->setpi(i, j, _pi[i][j]);
          getCell(ix, ny - 1, iz)->setpi(i, j, _pi[i][j]);
        }
      getCell(ix, ny - 2, iz)->setPi(_Pi);
      getCell(ix, ny - 1, iz)->setPi(_Pi);
    }
  // X
  for (int iy = 0; iy < ny; iy++)
    for (int iz = 0; iz < nz; iz++) {
      // left boundary
      getCell(2, iy, iz)->getQ(Q);
      getCell(1, iy, iz)->setQ(Q);
      getCell(0, iy, iz)->setQ(Q);
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) _pi[i][j] = getCell(2, iy, iz)->getpi(i, j);
      _Pi = getCell(2, iy, iz)->getPi();

      for (int i = 0; i < 4; i++)
        for (int j = 0; j <= i; j++) {
          getCell(0, iy, iz)->setpi(i, j, _pi[i][j]);
          getCell(1, iy, iz)->setpi(i, j, _pi[i][j]);
        }
      getCell(0, iy, iz)->setPi(_Pi);
      getCell(1, iy, iz)->setPi(_Pi);
      // right boundary
      getCell(nx - 3, iy, iz)->getQ(Q);
      getCell(nx - 2, iy, iz)->setQ(Q);
      getCell(nx - 1, iy, iz)->setQ(Q);
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
          _pi[i][j] = getCell(nx - 3, iy, iz)->getpi(i, j);
      _Pi = getCell(nx - 3, iy, iz)->getPi();

      for (int i = 0; i < 4; i++)
        for (int j = 0; j <= i; j++) {
          getCell(nx - 2, iy, iz)->setpi(i, j, _pi[i][j]);
          getCell(nx - 1, iy, iz)->setpi(i, j, _pi[i][j]);
        }
      getCell(nx - 2, iy, iz)->setPi(_Pi);
      getCell(nx - 1, iy, iz)->setPi(_Pi);
    }
}

void Fluid::updateM(double tau, double dt) {
  for (int ix = 0; ix < getNX(); ix++)
    for (int iy = 0; iy < getNY(); iy++)
      for (int iz = 0; iz < getNZ(); iz++) {
        Cell *c = getCell(ix, iy, iz);
        c->setDM(X_, 0.);
        c->setDM(Y_, 0.);
        c->setDM(Z_, 0.);
        if (getCell(ix, iy, iz)->getMaxM() < 1.) {
          if (getCell(ix + 1, iy, iz)->getM(X_) >= 1. ||
              getCell(ix - 1, iy, iz)->getM(X_) >= 1.)
            c->setDM(X_, dt / dx);
          if (getCell(ix, iy + 1, iz)->getM(Y_) >= 1. ||
              getCell(ix, iy - 1, iz)->getM(Y_) >= 1.)
            c->setDM(Y_, dt / dy);
          if (getCell(ix, iy, iz + 1)->getM(Z_) >= 1. ||
              getCell(ix, iy, iz - 1)->getM(Z_) >= 1.)
            c->setDM(Z_, dt / dz / tau);

          if (c->getDM(X_) == 0. && c->getDM(Y_) == 0.) {
            if (getCell(ix + 1, iy + 1, iz)->getMaxM() >= 1. ||
                getCell(ix + 1, iy - 1, iz)->getMaxM() >= 1. ||
                getCell(ix - 1, iy + 1, iz)->getMaxM() >= 1. ||
                getCell(ix - 1, iy - 1, iz)->getMaxM() >= 1.) {
              c->setDM(X_, 0.707 * dt / dx);
              c->setDM(Y_, 0.707 * dt / dy);
            }
          }
        }  // if
      }

  for (int ix = 0; ix < getNX(); ix++)
    for (int iy = 0; iy < getNY(); iy++)
      for (int iz = 0; iz < getNZ(); iz++) {
        Cell *c = getCell(ix, iy, iz);
        c->addM(X_, c->getDM(X_));
        c->addM(Y_, c->getDM(Y_));
        c->addM(Z_, c->getDM(Z_));
      }
}


void Fluid::outputGnuplot(double tau) {
  double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz;

  // X direction
  for (int ix = 0; ix < nx; ix++) {
    double x = getX(ix);
    Cell *c = getCell(ix, ny / 2, nz / 2);
    getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
    eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
    foutx << setw(14) << tau << setw(14) << x << setw(14) << vx << setw(14)
          << vy << setw(14) << e << setw(14) << nb << setw(14) << t << setw(14)
          << mub;
    for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
     foutx << setw(14) << c->getDbeta(i,j);
    foutx << endl;
  }
  foutx << endl;

  // Y direction
  for (int iy = 0; iy < ny; iy++) {
    double y = getY(iy);
    Cell *c = getCell(nx / 2, iy, nz / 2);
    getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
    eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
    fouty << setw(14) << tau << setw(14) << y << setw(14) << vy << setw(14)
          << vx << setw(14) << e << setw(14) << nb << setw(14) << t << setw(14)
          << mub;
    for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
     fouty << setw(14) << c->getDbeta(i,j);
    fouty << endl;
  }
  fouty << endl;

  // diagonal
  for (int ix = 0; ix < nx; ix++) {
    double x = getY(ix);
    Cell *c = getCell(ix, ix, nz / 2);
    getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
    eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
    foutdiag << setw(14) << tau << setw(14) << sqrt(2.) * x << setw(14) << vx
             << setw(14) << vy << setw(14) << e << setw(14) << nb << setw(14)
             << t << setw(14) << mub << endl;
    foutdiag << setw(14) << c->getpi(0, 0) << setw(14) << c->getpi(0, 1)
             << setw(14) << c->getpi(0, 2);
    foutdiag << setw(14) << c->getpi(0, 3) << setw(14) << c->getpi(1, 1)
             << setw(14) << c->getpi(1, 2);
    foutdiag << setw(14) << c->getpi(1, 3) << setw(14) << c->getpi(2, 2)
             << setw(14) << c->getpi(2, 3);
    foutdiag << setw(14) << c->getpi(3, 3) << setw(14) << c->getPi() << setw(14)
             << c->getViscCorrCutFlag() << endl;
  }
  foutdiag << endl;

  // Z direction
  for (int iz = 0; iz < nz; iz++) {
    double z = getZ(iz);
    Cell *c = getCell(nx / 2, ny / 2, iz);
    getCMFvariables(getCell(nx / 2, ny / 2, iz), tau, e, nb, nq, ns, vx, vy,
                    vz);
    eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
    foutz << setw(14) << tau << setw(14) << z << setw(14) << vz << setw(14)
          << vx << setw(14) << e << setw(14) << nb << setw(14) << t << setw(14)
          << mub;
    for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
     foutz << setw(14) << c->getDbeta(i,j);
    foutz << endl;
  }
  foutz << endl;
}

void Fluid::outputSnapshot(double tau) {
  double e, p, nb, nq, ns, T, mub, muq, mus, vx, vy, vz;

  // X direction
  for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
  for (int iz = 0; iz < nz; iz++) {
    double x = getX(ix);
    double y = getY(iy);
    double eta = getZ(iz);
    Cell *c = getCell(ix, iy, iz);
    getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
    eos->eos(e, nb, nq, ns, T, mub, muq, mus, p);
    if(T<0.01) T = 100.0;
    double gam = 1.0/sqrt(1.0 - vx*vx - vy*vy - tanh(vz)*tanh(vz));
    fout3d << setw(14) << tau << setw(14) << x << setw(14) << y << setw(14)
          << eta << setw(14) << e << setw(14) << vx << setw(14) << vy << setw(14)
          << vz;
    for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
     fout3d << setw(14) << c->getDbeta(i,j);
    fout3d << endl;
  }
  fout3d << endl;
}

// unput: geom. rapidity + velocities in Bjorken frame, --> output: velocities
// in lab.frame
void transformToLab(double eta, double &vx, double &vy, double &vz) {
  const double Y = eta + 1. / 2. * log((1. + vz) / (1. - vz));
  vx = vx * cosh(Y - eta) / cosh(Y);
  vy = vy * cosh(Y - eta) / cosh(Y);
  vz = tanh(Y);
}

// for the procedure below, the following approximations are used:
// dv/d(tau) = v^{t+dt}_ideal - v^{t}
// dv/dx_i ~ v^{x+dx}-v{x-dx},
// which makes sense after non-viscous step
void Fluid::NSquant(double tau, int ix, int iy, int iz, double pi[4][4],
                    double &Pi, double dmu[4][4], double dbeta[4][4],
                    double &du) {
  const double VMIN = 1e-2;
  const double UDIFF = 3.0;
  double e0, e1, p, nb, nq, ns, vx1, vy1, vz1, vx0, vy0, vz0, vxH, vyH, vzH;
  double T, T0, T1, mub, muq, mus;
  double ut0, ux0, uy0, uz0, ut1, ux1, uy1, uz1;
  //	double dmu [4][4] ; // \partial_\mu u^\nu matrix
  // coordinates: 0=tau, 1=x, 2=y, 3=eta
  double Z[4][4][4][4];  // Z[mu][nu][lambda][rho]
  double uuu[4];         // the 4-velocity
  double gmunu[4][4] = {{1, 0, 0, 0},
                        {0, -1, 0, 0},
                        {0, 0, -1, 0},
                        {0, 0, 0, -1}};  // omit 1/tau^2 in g^{eta,eta}
  Cell *c = getCell(ix, iy, iz);
  // check if the cell is next to vacuum from +-x, +-y side:
  if (c->getNext(X_)->getMaxM() <= 0.9 || c->getNext(Y_)->getMaxM() <= 0.9 ||
      c->getPrev(X_)->getMaxM() <= 0.9 || c->getPrev(Y_)->getMaxM() <= 0.9 ||
      getCell(ix + 1, iy + 1, iz)->getMaxM() <= 0.9 ||
      getCell(ix + 1, iy - 1, iz)->getMaxM() <= 0.9 ||
      getCell(ix - 1, iy + 1, iz)->getMaxM() <= 0.9 ||
      getCell(ix - 1, iy - 1, iz)->getMaxM() <= 0.9) {
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        pi[i][j] = 0.;
        dmu[i][j] = 0.;
        dbeta[i][j] = 0.;
      }
    Pi = du = 0.;
    return;
  }
  // calculation of \partial_\mu u^\nu matrix
  // mu=first index, nu=second index
  // centered differences with respect to the values at (it+1/2, ix, iy, iz)
  // d_tau u^\mu
  c->getPrimVarPrev(eos, tau - dt, e0, p, nb, nq, ns, vx0, vy0, vz0);
  eos->eos(e0, nb, nq, ns, T0, mub, muq, mus, p);
  c->getPrimVar(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1);
  eos->eos(e1, nb, nq, ns, T1, mub, muq, mus, p);
  c->getPrimVarHCenter(eos, tau - 0.5 * dt, e1, p, nb, nq, ns, vxH, vyH, vzH);
  //############## get transport coefficients
  double etaS, zetaS;
  double s = eos->s(e1, nb, nq, ns);  // entropy density in the current cell
  eos->eos(e1, nb, nq, ns, T, mub, muq, mus, p);
  trcoeff->getEta(e1, T, etaS, zetaS);
  //##############
  // if(e1<0.00004) s=0. ; // negative pressure due to pi^zz for small e
  ut0 = 1.0 / sqrt(1.0 - vx0 * vx0 - vy0 * vy0 - vz0 * vz0);
  ux0 = ut0 * vx0;
  uy0 = ut0 * vy0;
  uz0 = ut0 * vz0;
  ut1 = 1.0 / sqrt(1.0 - vx1 * vx1 - vy1 * vy1 - vz1 * vz1);
  ux1 = ut1 * vx1;
  uy1 = ut1 * vy1;
  uz1 = ut1 * vz1;
  uuu[0] = 1.0 / sqrt(1.0 - vxH * vxH - vyH * vyH - vzH * vzH);
  uuu[1] = uuu[0] * vxH;
  uuu[2] = uuu[0] * vyH;
  uuu[3] = uuu[0] * vzH;

  dmu[0][0] = (ut1 * ut1 - ut0 * ut0) / 2. / uuu[0] / dt;
  dmu[0][1] = (ux1 * ux1 - ux0 * ux0) / 2. / uuu[1] / dt;
  dmu[0][2] = (uy1 * uy1 - uy0 * uy0) / 2. / uuu[2] / dt;
  dmu[0][3] = (uz1 * uz1 - uz0 * uz0) / 2. / uuu[3] / dt;
  dbeta[0][0] = (ut1 / T1 - ut0 / T0) / dt;
  dbeta[0][1] = (ux1 / T1 - ux0 / T0) / dt;
  dbeta[0][2] = (uy1 / T1 - uy0 / T0) / dt;
  dbeta[0][3] = (uz1 / T1 - uz0 / T0) / dt;
  if (fabs(0.5 * (ut1 + ut0) / ut1) > UDIFF) dmu[0][0] = (ut1 - ut0) / dt;
  if (fabs(uuu[1]) < VMIN || fabs(0.5 * (ux1 + ux0) / ux1) > UDIFF)
    dmu[0][1] = (ux1 - ux0) / dt;
  if (fabs(uuu[2]) < VMIN || fabs(0.5 * (uy1 + uy0) / uy1) > UDIFF)
    dmu[0][2] = (uy1 - uy0) / dt;
  if (fabs(uuu[3]) < VMIN || fabs(0.5 * (uz1 + uz0) / uz1) > UDIFF)
    dmu[0][3] = (uz1 - uz0) / dt;
  if (e1 <= 0. || e0 <= 0.) {  // matter-vacuum
    dmu[0][0] = dmu[0][1] = dmu[0][2] = dmu[0][3] = 0.;
  }
  if(e1 <= 0. || e0 <= 0. || T1<=0. || T0<=0.) {
    dbeta[0][0] = dbeta[0][1] = dbeta[0][2] = dbeta[0][3] = 0.;
  }
  // d_x u^\mu
  getCell(ix + 1, iy, iz)
      ->getPrimVarHCenter(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1);
  eos->eos(e1, nb, nq, ns, T1, mub, muq, mus, p);
  getCell(ix - 1, iy, iz)
      ->getPrimVarHCenter(eos, tau, e0, p, nb, nq, ns, vx0, vy0, vz0);
  eos->eos(e0, nb, nq, ns, T0, mub, muq, mus, p);
  if (e1 > 0. && e0 > 0.) {
    ut0 = 1.0 / sqrt(1.0 - vx0 * vx0 - vy0 * vy0 - vz0 * vz0);
    ux0 = ut0 * vx0;
    uy0 = ut0 * vy0;
    uz0 = ut0 * vz0;
    ut1 = 1.0 / sqrt(1.0 - vx1 * vx1 - vy1 * vy1 - vz1 * vz1);
    ux1 = ut1 * vx1;
    uy1 = ut1 * vy1;
    uz1 = ut1 * vz1;
    dmu[1][0] = 0.25 * (ut1 * ut1 - ut0 * ut0) / uuu[0] / dx;
    dmu[1][1] = 0.25 * (ux1 * ux1 - ux0 * ux0) / uuu[1] / dx;
    dmu[1][2] = 0.25 * (uy1 * uy1 - uy0 * uy0) / uuu[2] / dx;
    dmu[1][3] = 0.25 * (uz1 * uz1 - uz0 * uz0) / uuu[3] / dx;
    dbeta[1][0] = 0.5 * (ut1 / T1 - ut0 / T0) / dx;
    dbeta[1][1] = 0.5 * (ux1 / T1 - ux0 / T0) / dx;
    dbeta[1][2] = 0.5 * (uy1 / T1 - uy0 / T0) / dx;
    dbeta[1][3] = 0.5 * (uz1 / T1 - uz0 / T0) / dx;
    if (fabs(0.5 * (ut1 + ut0) / uuu[0]) > UDIFF)
      dmu[1][0] = 0.5 * (ut1 - ut0) / dx;
    if (fabs(uuu[1]) < VMIN || fabs(0.5 * (ux1 + ux0) / uuu[1]) > UDIFF)
      dmu[1][1] = 0.5 * (ux1 - ux0) / dx;
    if (fabs(uuu[2]) < VMIN || fabs(0.5 * (uy1 + uy0) / uuu[2]) > UDIFF)
      dmu[1][2] = 0.5 * (uy1 - uy0) / dx;
    if (fabs(uuu[3]) < VMIN || fabs(0.5 * (uz1 + uz0) / uuu[3]) > UDIFF)
      dmu[1][3] = 0.5 * (uz1 - uz0) / dx;
  } else {  // matter-vacuum
    dmu[1][0] = dmu[1][1] = dmu[1][2] = dmu[1][3] = 0.;
  }
  if(e1 <= 0. || e0 <= 0. || T1<=0. || T0<=0.) {
    dbeta[1][0] = dbeta[1][1] = dbeta[1][2] = dbeta[1][3] = 0.;
  }
  if (fabs(dmu[1][3]) > 1e+10)
    cout << "dmu[1][3]:  " << uz1 << "  " << uz0 << "  " << uuu[3] << endl;
  // d_y u^\mu
  getCell(ix, iy + 1, iz)
      ->getPrimVarHCenter(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1);
  eos->eos(e1, nb, nq, ns, T1, mub, muq, mus, p);
  getCell(ix, iy - 1, iz)
      ->getPrimVarHCenter(eos, tau, e0, p, nb, nq, ns, vx0, vy0, vz0);
  eos->eos(e0, nb, nq, ns, T0, mub, muq, mus, p);
  if (e1 > 0. && e0 > 0.) {
    ut0 = 1.0 / sqrt(1.0 - vx0 * vx0 - vy0 * vy0 - vz0 * vz0);
    ux0 = ut0 * vx0;
    uy0 = ut0 * vy0;
    uz0 = ut0 * vz0;
    ut1 = 1.0 / sqrt(1.0 - vx1 * vx1 - vy1 * vy1 - vz1 * vz1);
    ux1 = ut1 * vx1;
    uy1 = ut1 * vy1;
    uz1 = ut1 * vz1;
    dmu[2][0] = 0.25 * (ut1 * ut1 - ut0 * ut0) / uuu[0] / dy;
    dmu[2][1] = 0.25 * (ux1 * ux1 - ux0 * ux0) / uuu[1] / dy;
    dmu[2][2] = 0.25 * (uy1 * uy1 - uy0 * uy0) / uuu[2] / dy;
    dmu[2][3] = 0.25 * (uz1 * uz1 - uz0 * uz0) / uuu[3] / dy;
    dbeta[2][0] = 0.5 * (ut1 / T1 - ut0 / T0) / dy;
    dbeta[2][1] = 0.5 * (ux1 / T1 - ux0 / T0) / dy;
    dbeta[2][2] = 0.5 * (uy1 / T1 - uy0 / T0) / dy;
    dbeta[2][3] = 0.5 * (uz1 / T1 - uz0 / T0) / dy;
    if (fabs(0.5 * (ut1 + ut0) / uuu[0]) > UDIFF)
      dmu[2][0] = 0.5 * (ut1 - ut0) / dy;
    if (fabs(uuu[1]) < VMIN || fabs(0.5 * (ux1 + ux0) / uuu[1]) > UDIFF)
      dmu[2][1] = 0.5 * (ux1 - ux0) / dy;
    if (fabs(uuu[2]) < VMIN || fabs(0.5 * (uy1 + uy0) / uuu[2]) > UDIFF)
      dmu[2][2] = 0.5 * (uy1 - uy0) / dy;
    if (fabs(uuu[3]) < VMIN || fabs(0.5 * (uz1 + uz0) / uuu[3]) > UDIFF)
      dmu[2][3] = 0.5 * (uz1 - uz0) / dy;
  } else {  // matter-vacuum
    dmu[2][0] = dmu[2][1] = dmu[2][2] = dmu[2][3] = 0.;
  }
  if(e1 <= 0. || e0 <= 0. || T1<=0. || T0<=0.){
    dbeta[2][0] = dbeta[2][1] = dbeta[2][2] = dbeta[2][3] = 0.;
  }
  // d_z u^\mu
  getCell(ix, iy, iz + 1)
      ->getPrimVarHCenter(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1);
  eos->eos(e1, nb, nq, ns, T1, mub, muq, mus, p);
  getCell(ix, iy, iz - 1)
      ->getPrimVarHCenter(eos, tau, e0, p, nb, nq, ns, vx0, vy0, vz0);
  eos->eos(e0, nb, nq, ns, T0, mub, muq, mus, p);
  if (e1 > 0. && e0 > 0.) {
    ut0 = 1.0 / sqrt(1.0 - vx0 * vx0 - vy0 * vy0 - vz0 * vz0);
    ux0 = ut0 * vx0;
    uy0 = ut0 * vy0;
    uz0 = ut0 * vz0;
    ut1 = 1.0 / sqrt(1.0 - vx1 * vx1 - vy1 * vy1 - vz1 * vz1);
    ux1 = ut1 * vx1;
    uy1 = ut1 * vy1;
    uz1 = ut1 * vz1;
    dmu[3][0] = 0.25 * (ut1 * ut1 - ut0 * ut0) / uuu[0] / dz / (tau - 0.5 * dt);
    dmu[3][1] = 0.25 * (ux1 * ux1 - ux0 * ux0) / uuu[1] / dz / (tau - 0.5 * dt);
    dmu[3][2] = 0.25 * (uy1 * uy1 - uy0 * uy0) / uuu[2] / dz / (tau - 0.5 * dt);
    dmu[3][3] = 0.25 * (uz1 * uz1 - uz0 * uz0) / uuu[3] / dz / (tau - 0.5 * dt);
    dbeta[3][0] = 0.5 * (ut1 / T1 - ut0 / T0) / (dz * (tau - 0.5 * dt));
    dbeta[3][1] = 0.5 * (ux1 / T1 - ux0 / T0) / (dz * (tau - 0.5 * dt));
    dbeta[3][2] = 0.5 * (uy1 / T1 - uy0 / T0) / (dz * (tau - 0.5 * dt));
    dbeta[3][3] = 0.5 * (uz1 / T1 - uz0 / T0) / (dz * (tau - 0.5 * dt));
    if (fabs(0.5 * (ut1 + ut0) / uuu[0]) > UDIFF)
      dmu[3][0] = 0.5 * (ut1 - ut0) / dz / (tau + 0.5 * dt);
    if (fabs(uuu[1]) < VMIN || fabs(0.5 * (ux1 + ux0) / uuu[1]) > UDIFF)
      dmu[3][1] = 0.5 * (ux1 - ux0) / dz / (tau + 0.5 * dt);
    if (fabs(uuu[2]) < VMIN || fabs(0.5 * (uy1 + uy0) / uuu[2]) > UDIFF)
      dmu[3][2] = 0.5 * (uy1 - uy0) / dz / (tau + 0.5 * dt);
    if (fabs(uuu[3]) < VMIN || fabs(0.5 * (uz1 + uz0) / uuu[3]) > UDIFF)
      dmu[3][3] = 0.5 * (uz1 - uz0) / dz / (tau + 0.5 * dt);
  } else {  // matter-vacuum
    dmu[3][0] = dmu[3][1] = dmu[3][2] = dmu[3][3] = 0.;
  }
  if(e1 <= 0. || e0 <= 0. || T1<=0. || T0<=0.){
    dbeta[3][0] = dbeta[3][1] = dbeta[3][2] = dbeta[3][3] = 0.;
  }
  // additional terms from Christoffel symbols :)
  dmu[3][0] += uuu[3] / (tau - 0.5 * dt);
  dmu[3][3] += uuu[0] / (tau - 0.5 * dt);
  if(T>0.){
  dbeta[3][0] += uuu[3] / (T * (tau - 0.5 * dt));
  dbeta[3][3] += uuu[0] / (T * (tau - 0.5 * dt));
  }
  // calculation of Z[mu][nu][lambda][rho]
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      for (int k = 0; k < 4; k++)
        for (int l = 0; l < 4; l++) Z[i][j][k][l] = 0.0;
  // filling Z matrix
  for (int mu = 0; mu < 4; mu++)
    for (int nu = 0; nu < 4; nu++)
      for (int lam = 0; lam < 4; lam++)
        for (int rho = 0; rho < 4; rho++) {
          if (nu == rho)
            Z[mu][nu][lam][rho] += 0.5 * (gmunu[mu][lam] - uuu[mu] * uuu[lam]);
          if (mu == rho)
            Z[mu][nu][lam][rho] += 0.5 * (gmunu[nu][lam] - uuu[nu] * uuu[lam]);
          if (lam == rho)
            Z[mu][nu][lam][rho] -= (gmunu[mu][nu] - uuu[mu] * uuu[nu]) / 3.0;
        }
  // calculating sigma[mu][nu]
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) {
      pi[i][j] = 0.0;
      for (int k = 0; k < 4; k++)
        for (int l = 0; l < 4; l++) {
          pi[i][j] += Z[i][j][k][l] * dmu[k][l] * 2.0 * etaS * s / 5.068;
        }
    }
  Pi = -zetaS * s * (dmu[0][0] + dmu[1][1] + dmu[2][2] + dmu[3][3]) /
       5.068;  // fm^{-4} --> GeV/fm^3
  du = dmu[0][0] + dmu[1][1] + dmu[2][2] + dmu[3][3];
  //--------- debug part: NaN/inf check, trace check, diag check, transversality
  //check
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) {
      if (pi[i][j] != 0. && fabs(1.0 - pi[j][i] / pi[i][j]) > 1e-10)
        cout << "non-diag: " << pi[i][j] << "  " << pi[j][i] << endl;
      if (isinf(pi[i][j]) || isnan(pi[i][j])) {
        cout << "hydro:NSquant: inf/nan i " << i << " j " << j << endl;
        exit(1);
      }
    }
}

void Fluid::outputSurface(double tau) {
  static double nbSurf = 0.0;
  double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz, Q[7];
  double E = 0., Efull = 0., S = 0., Px = 0., vt_num = 0., vt_den = 0.,
         vxvy_num = 0., vxvy_den = 0., pi0x_num = 0., pi0x_den = 0.,
         txxyy_num = 0., txxyy_den = 0., Nb1 = 0., Nb2 = 0.;
  double eta = 0;
  int nelements = 0, nsusp = 0;  // all surface emenents and suspicious ones
  int nCoreCells = 0,
      nCoreCutCells = 0;  // cells with e>eCrit and cells with cut visc.corr.
                          //-- Cornelius: allocating memory for corner points
  double ****ccube = new double ***[2];
  for (int i1 = 0; i1 < 2; i1++) {
    ccube[i1] = new double **[2];
    for (int i2 = 0; i2 < 2; i2++) {
      ccube[i1][i2] = new double *[2];
      for (int i3 = 0; i3 < 2; i3++) {
        ccube[i1][i2][i3] = new double[2];
      }
    }
  }
//----end Cornelius
#ifdef SWAP_EOS
  swap(eos, eosH);
#endif
  fout2d << endl;
  for (int ix = 2; ix < nx - 2; ix++)
    for (int iy = 2; iy < ny - 2; iy++)
      for (int iz = 2; iz < nz - 2; iz++) {
        Cell *c = getCell(ix, iy, iz);
        getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
        c->getQ(Q);
        eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
        double s = eos->s(e, nb, nq, ns);
        eta = getZ(iz);
        const double cosh_int =
            (sinh(eta + 0.5 * dz) - sinh(eta - 0.5 * dz)) / dz;
        const double sinh_int =
            (cosh(eta + 0.5 * dz) - cosh(eta - 0.5 * dz)) / dz;
        E += tau * (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                 (cosh_int - tanh(vz) * sinh_int) -
             tau * p * cosh_int;
        Nb1 += Q[NB_];
        Nb2 += tau * nb * (cosh_int - tanh(vz) * sinh_int) /
               sqrt(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
        //---- inf check
        if (isinf(E)) {
          cout << "EEinf" << setw(14) << e << setw(14) << p << setw(14) << vx
               << setw(14) << vy << setw(14) << vz << endl;
          exit(1);
        }
        //--------------
        Efull += tau * (e + p) /
                     (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                     (cosh(eta) - tanh(vz) * sinh(eta)) -
                 tau * p * cosh(eta);
        if (trcoeff->isViscous())
          Efull += tau * c->getpi(0, 0) * cosh(eta) +
                   tau * c->getpi(0, 3) * sinh(eta);
        if (e > ecrit) {
          nCoreCells++;
          if (c->getViscCorrCutFlag() < 0.9) nCoreCutCells++;
        }
        // -- noneq. corrections to entropy flux
        const double gmumu[4] = {1., -1., -1., -1.};
        double deltas = 0.;
        if (trcoeff->isViscous())
          for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
              deltas += pow(c->getpi(i, j), 2) * gmumu[i] * gmumu[j];
        if (t > 0.02) {
          s += 1.5 * deltas / ((e + p) * t);
          S += tau * s * (cosh_int - tanh(vz) * sinh_int) /
               sqrt(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
        }
        Px +=
            tau * (e + p) * vx / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
        vt_num += e / sqrt(1. - vx * vx - vy * vy) * sqrt(vx * vx + vy * vy);
        vt_den += e / sqrt(1. - vx * vx - vy * vy);
        vxvy_num += e * (fabs(vx) - fabs(vy));
        vxvy_den += e;
        txxyy_num += (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                     (vx * vx - vy * vy);
        txxyy_den += (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                         (vx * vx + vy * vy) +
                     2. * p;
        pi0x_num += e / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                    fabs(c->getpi(0, 1));
        pi0x_den += e / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
        //----- Cornelius stuff
        double QCube[2][2][2][2][7];
        double piSquare[2][2][2][10], PiSquare[2][2][2];
        double dbetaSq [2][2][2][4][4];
        for (int jx = 0; jx < 2; jx++)
          for (int jy = 0; jy < 2; jy++)
            for (int jz = 0; jz < 2; jz++) {
              double _p, _nb, _nq, _ns, _vx, _vy, _vz;
              Cell *cc = getCell(ix + jx, iy + jy, iz + jz);
              cc->getPrimVar(eos, tau, e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
              cc->getQ(QCube[1][jx][jy][jz]);
              ccube[1][jx][jy][jz] = e;
              cc->getPrimVarPrev(eos, tau - dt, e, _p, _nb, _nq, _ns, _vx, _vy,
                                 _vz);
              cc->getQprev(QCube[0][jx][jy][jz]);
              ccube[0][jx][jy][jz] = e;
              // ---- get viscous tensor
              for (int ii = 0; ii < 4; ii++)
                for (int jj = 0; jj <= ii; jj++)
                  piSquare[jx][jy][jz][index44(ii, jj)] = cc->getpi(ii, jj);
              for (int i=0; i<4; i++)
              for (int j=0; j<4; j++)
               dbetaSq[jx][jy][jz][i][j] = cc->getDbeta(i,j);
              PiSquare[jx][jy][jz] = cc->getPi();
            }
        cornelius->find_surface_4d(ccube);
        const int Nsegm = cornelius->get_Nelements();
        for (int isegm = 0; isegm < Nsegm; isegm++) {
          nelements++;
          ffreeze.precision(15);
          ffreeze << setw(24) << tau + cornelius->get_centroid_elem(isegm, 0)
                  << setw(24)
                  << getX(ix) + cornelius->get_centroid_elem(isegm, 1)
                  << setw(24)
                  << getY(iy) + cornelius->get_centroid_elem(isegm, 2)
                  << setw(24)
                  << getZ(iz) + cornelius->get_centroid_elem(isegm, 3);
          // ---- interpolation procedure
          double vxC = 0., vyC = 0., vzC = 0., TC = 0., mubC = 0., muqC = 0.,
                 musC = 0., piC[10], PiC = 0., nbC = 0.,
                 nqC = 0.;  // values at the centre, to be interpolated
          double QC[7] = {0., 0., 0., 0., 0., 0., 0.};
          double eC = 0., pC = 0.;
          for (int ii = 0; ii < 10; ii++) piC[ii] = 0.0;
          double wCenT[2] = {1. - cornelius->get_centroid_elem(isegm, 0) / dt,
                             cornelius->get_centroid_elem(isegm, 0) / dt};
          double wCenX[2] = {1. - cornelius->get_centroid_elem(isegm, 1) / dx,
                             cornelius->get_centroid_elem(isegm, 1) / dx};
          double wCenY[2] = {1. - cornelius->get_centroid_elem(isegm, 2) / dy,
                             cornelius->get_centroid_elem(isegm, 2) / dy};
          double wCenZ[2] = {1. - cornelius->get_centroid_elem(isegm, 3) / dz,
                             cornelius->get_centroid_elem(isegm, 3) / dz};
          for (int jt = 0; jt < 2; jt++)
            for (int jx = 0; jx < 2; jx++)
              for (int jy = 0; jy < 2; jy++)
                for (int jz = 0; jz < 2; jz++)
                  for (int i = 0; i < 7; i++) {
                    QC[i] += QCube[jt][jx][jy][jz][i] * wCenT[jt] * wCenX[jx] *
                             wCenY[jy] * wCenZ[jz];
                  }
          for (int i = 0; i < 7; i++)
            QC[i] = QC[i] / (tau + cornelius->get_centroid_elem(isegm, 0));
          double _ns = 0.0;
          transformPV(eos, QC, eC, pC, nbC, nqC, _ns, vxC, vyC, vzC);
          eos->eos(eC, nbC, nqC, _ns, TC, mubC, muqC, musC, pC);
          if (TC > 0.4 || fabs(mubC) > 0.85) {
            cout << "#### Error (surface): high T/mu_b ####\n";
          }
          if (eC > ecrit * 2.0 || eC < ecrit * 0.5) nsusp++;
          double dbetaC [4][4] = {0.};
          for (int jx = 0; jx < 2; jx++)
            for (int jy = 0; jy < 2; jy++)
              for (int jz = 0; jz < 2; jz++) {
                for (int ii = 0; ii < 10; ii++)
                  piC[ii] += piSquare[jx][jy][jz][ii] * wCenX[jx] * wCenY[jy] *
                             wCenZ[jz];
                PiC += PiSquare[jx][jy][jz] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
                for(int i=0; i<4; i++)
                for(int j=0; j<4; j++)
                 dbetaC[i][j] += dbetaSq[jx][jy][jz][i][j]*wCenX[jx]*wCenY[jy]*wCenZ[jz];
              }
          double v2C = vxC * vxC + vyC * vyC + vzC * vzC;
          if (v2C > 1.) {
            vxC *= sqrt(0.99 / v2C);
            vyC *= sqrt(0.99 / v2C);
            vzC *= sqrt(0.99 / v2C);
            v2C = 0.99;
          }
          double etaC = getZ(iz) + cornelius->get_centroid_elem(isegm, 3);
          transformToLab(etaC, vxC, vyC, vzC);  // viC is now in lab.frame!
          double gammaC = 1. / sqrt(1. - vxC * vxC - vyC * vyC - vzC * vzC);

          double uC[4] = {gammaC, gammaC * vxC, gammaC * vyC, gammaC * vzC};
          const double tauC = tau + cornelius->get_centroid_elem(isegm, 0);
          double dsigma[4];
          // ---- transform dsigma to lab.frame :
          const double ch = cosh(etaC);
          const double sh = sinh(etaC);
          dsigma[0] = tauC * (ch * cornelius->get_normal_elem(0, 0) -
                              sh / tauC * cornelius->get_normal_elem(0, 3));
          dsigma[3] = tauC * (-sh * cornelius->get_normal_elem(0, 0) +
                              ch / tauC * cornelius->get_normal_elem(0, 3));
          dsigma[1] = tauC * cornelius->get_normal_elem(0, 1);
          dsigma[2] = tauC * cornelius->get_normal_elem(0, 2);
          double dVEff = 0.0;
          for (int ii = 0; ii < 4; ii++)
            dVEff += dsigma[ii] * uC[ii];  // normalize for Delta eta=1
          vEff += dVEff;
          for (int ii = 0; ii < 4; ii++) ffreeze << setw(24) << dsigma[ii];
          for (int ii = 0; ii < 4; ii++) ffreeze << setw(24) << uC[ii];
          ffreeze << setw(24) << TC << setw(24) << mubC << setw(24) << muqC
                  << setw(24) << musC;
          double picart[10];
          /*pi00*/ picart[index44(0, 0)] = ch * ch * piC[index44(0, 0)] +
                                           2. * ch * sh * piC[index44(0, 3)] +
                                           sh * sh * piC[index44(3, 3)];
          /*pi01*/ picart[index44(0, 1)] =
              ch * piC[index44(0, 1)] + sh * piC[index44(3, 1)];
          /*pi02*/ picart[index44(0, 2)] =
              ch * piC[index44(0, 2)] + sh * piC[index44(3, 2)];
          /*pi03*/ picart[index44(0, 3)] =
              ch * sh * (piC[index44(0, 0)] + piC[index44(3, 3)]) +
              (ch * ch + sh * sh) * piC[index44(0, 3)];
          /*pi11*/ picart[index44(1, 1)] = piC[index44(1, 1)];
          /*pi12*/ picart[index44(1, 2)] = piC[index44(1, 2)];
          /*pi13*/ picart[index44(1, 3)] =
              sh * piC[index44(0, 1)] + ch * piC[index44(3, 1)];
          /*pi22*/ picart[index44(2, 2)] = piC[index44(2, 2)];
          /*pi23*/ picart[index44(2, 3)] =
              sh * piC[index44(0, 2)] + ch * piC[index44(3, 2)];
          /*pi33*/ picart[index44(3, 3)] = sh * sh * piC[index44(0, 0)] +
                                           ch * ch * piC[index44(3, 3)] +
                                           2. * sh * ch * piC[index44(0, 3)];
#ifdef OUTPI
          for (int ii = 0; ii < 10; ii++) ffreeze << setw(24) << picart[ii];
          ffreeze << setw(24) << PiC ;
#endif
          const double jacob [4][4] =
          {{ch, 0., 0., -sh}, {0., 1., 0., 0.}, {0., 0., 1., 0.},
           {-sh, 0., 0., ch}}; // Jacobian to transform covariant (lower index)
           // vector from Milne to Cartesian coordinate system
          double dbetaCart [4][4] = {0};
          for(int i=0; i<4; i++)
          for(int j=0; j<4; j++)
          for(int k=0; k<4; k++)
          for(int l=0; l<4; l++)
           dbetaCart [i][j] += jacob[i][k] * jacob[j][l] * dbetaC[k][l]
                               * gmumu[l]; // which equals to d_i beta_j
          for(int i=0; i<4; i++)
          for(int j=0; j<4; j++)
           ffreeze << setw(24) << dbetaCart[i][j];
          ffreeze << endl;
          double dEsurfVisc = 0.;
          for (int i = 0; i < 4; i++)
            dEsurfVisc += picart[index44(0, i)] * dsigma[i];
          EtotSurf += (eC + pC) * uC[0] * dVEff - pC * dsigma[0] + dEsurfVisc;
          nbSurf += nbC * dVEff;
        }
        // if(cornelius->get_Nelements()>1) cout << "oops, Nelements>1\n" ;
        //----- end Cornelius
      }
  E = E * dx * dy * dz;
  Efull = Efull * dx * dy * dz;
  S = S * dx * dy * dz;
  Nb1 *= dx * dy * dz;
  Nb2 *= dx * dy * dz;
  fout_aniz << setw(12) << tau << setw(14) << vt_num / vt_den << setw(14)
            << vxvy_num / vxvy_den << setw(14) << pi0x_num / pi0x_den << endl;
  cout << setw(10) << tau << setw(13) << E << setw(13) << Efull << setw(13)
       << nbSurf << setw(13) << S << setw(10) << nelements << setw(10) << nsusp
       << setw(13) << (float)(nCoreCutCells) / (float)(nCoreCells) << endl;
  //-- Cornelius: all done, let's free memory
  for (int i1 = 0; i1 < 2; i1++) {
    for (int i2 = 0; i2 < 2; i2++) {
      for (int i3 = 0; i3 < 2; i3++) {
        delete[] ccube[i1][i2][i3];
      }
      delete[] ccube[i1][i2];
    }
    delete[] ccube[i1];
  }
  delete[] ccube;
//----end Cornelius
#ifdef SWAP_EOS
  swap(eos, eosH);
#endif
  if (nelements == 0) exit(0);
}

void Fluid::outputCorona(double tau) {
  static double nbSurf = 0.0;
  double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz, Q[7];
  double E = 0., Efull = 0., S = 0., Px = 0., vt_num = 0., vt_den = 0.,
         vxvy_num = 0., vxvy_den = 0., pi0x_num = 0., pi0x_den = 0.,
         txxyy_num = 0., txxyy_den = 0., Nb1 = 0., Nb2 = 0.;
  double eta = 0;
  int nelements = 0;

#ifdef SWAP_EOS
  swap(eos, eosH);
#endif
  fout2d << endl;
  for (int ix = 2; ix < nx - 2; ix++)
    for (int iy = 2; iy < ny - 2; iy++)
      for (int iz = 2; iz < nz - 2; iz++) {
        Cell *c = getCell(ix, iy, iz);
        getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
        c->getQ(Q);
        eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
        double s = eos->s(e, nb, nq, ns);
        eta = getZ(iz);
        const double cosh_int =
            (sinh(eta + 0.5 * dz) - sinh(eta - 0.5 * dz)) / dz;
        const double sinh_int =
            (cosh(eta + 0.5 * dz) - cosh(eta - 0.5 * dz)) / dz;
        E += tau * (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                 (cosh_int - tanh(vz) * sinh_int) -
             tau * p * cosh_int;
        Nb1 += Q[NB_];
        Nb2 += tau * nb * (cosh_int - tanh(vz) * sinh_int) /
               sqrt(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
        //---- inf check
        if (isinf(E)) {
          cout << "EEinf" << setw(14) << e << setw(14) << p << setw(14) << vx
               << setw(14) << vy << setw(14) << vz << endl;
          exit(1);
        }
        //--------------
        Efull += tau * (e + p) /
                     (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                     (cosh(eta) - tanh(vz) * sinh(eta)) -
                 tau * p * cosh(eta);
        if (trcoeff->isViscous())
          Efull += tau * c->getpi(0, 0) * cosh(eta) +
                   tau * c->getpi(0, 3) * sinh(eta);
        // -- noneq. corrections to entropy flux
        const double gmumu[4] = {1., -1., -1., -1.};
        double deltas = 0.;
        if (trcoeff->isViscous())
          for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
              deltas += pow(c->getpi(i, j), 2) * gmumu[i] * gmumu[j];
        if (t > 0.02) {
          s += 1.5 * deltas / ((e + p) * t);
          S += tau * s * (cosh_int - tanh(vz) * sinh_int) /
               sqrt(1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
        }
        Px +=
            tau * (e + p) * vx / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));
        vt_num += e / sqrt(1. - vx * vx - vy * vy) * sqrt(vx * vx + vy * vy);
        vt_den += e / sqrt(1. - vx * vx - vy * vy);
        vxvy_num += e * (fabs(vx) - fabs(vy));
        vxvy_den += e;
        txxyy_num += (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                     (vx * vx - vy * vy);
        txxyy_den += (e + p) / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                         (vx * vx + vy * vy) +
                     2. * p;
        pi0x_num += e / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz)) *
                    fabs(c->getpi(0, 1));
        pi0x_den += e / (1. - vx * vx - vy * vy - tanh(vz) * tanh(vz));

        //----- Cornelius stuff
        bool isCorona = true, isTail = true;
        double QCube[2][2][2][7];
        double piSquare[2][2][2][10], PiSquare[2][2][2];
        for (int jx = 0; jx < 2; jx++)
          for (int jy = 0; jy < 2; jy++)
            for (int jz = 0; jz < 2; jz++) {
              double _p, _nb, _nq, _ns, _vx, _vy, _vz;
              Cell *cc = getCell(ix + jx, iy + jy, iz + jz);
              cc->getPrimVar(eos, tau, e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
              cc->getQ(QCube[jx][jy][jz]);
              if (e > ecrit) isCorona = false;
              if (e > 0.01) isTail = false;
              // ---- get viscous tensor
              for (int ii = 0; ii < 4; ii++)
                for (int jj = 0; jj <= ii; jj++)
                  piSquare[jx][jy][jz][index44(ii, jj)] = cc->getpi(ii, jj);
              PiSquare[jx][jy][jz] = cc->getPi();
            }
        if (isCorona && !isTail) {
          nelements++;
          ffreeze.precision(15);
          ffreeze << setw(24) << tau << setw(24) << getX(ix) + 0.5 * dx
                  << setw(24) << getY(iy) + 0.5 * dy << setw(24)
                  << getZ(iz) + 0.5 * dz;
          // ---- interpolation procedure
          double vxC = 0., vyC = 0., vzC = 0., TC = 0., mubC = 0., muqC = 0.,
                 musC = 0., piC[10], PiC = 0., nbC = 0.,
                 nqC = 0.;  // values at the centre, to be interpolated
          double QC[7] = {0., 0., 0., 0., 0., 0., 0.};
          double eC = 0., pC = 0.;
          for (int ii = 0; ii < 10; ii++) piC[ii] = 0.0;
          for (int jx = 0; jx < 2; jx++)
            for (int jy = 0; jy < 2; jy++)
              for (int jz = 0; jz < 2; jz++)
                for (int i = 0; i < 7; i++) {
                  QC[i] += QCube[jx][jy][jz][i] * 0.125;
                }
          for (int i = 0; i < 7; i++) QC[i] = QC[i] / tau;
          double _ns = 0.0;
          transformPV(eos, QC, eC, pC, nbC, nqC, _ns, vxC, vyC, vzC);
          eos->eos(eC, nbC, nqC, _ns, TC, mubC, muqC, musC, pC);
          if (TC > 0.4 || fabs(mubC) > 0.85) {
            cout << "#### Error (surface): high T/mu_b ####\n";
          }
          for (int jx = 0; jx < 2; jx++)
            for (int jy = 0; jy < 2; jy++)
              for (int jz = 0; jz < 2; jz++) {
                for (int ii = 0; ii < 10; ii++)
                  piC[ii] += piSquare[jx][jy][jz][ii] * 0.125;
                PiC += PiSquare[jx][jy][jz] * 0.125;
              }
          double v2C = vxC * vxC + vyC * vyC + vzC * vzC;
          if (v2C > 1.) {
            vxC *= sqrt(0.99 / v2C);
            vyC *= sqrt(0.99 / v2C);
            vzC *= sqrt(0.99 / v2C);
            v2C = 0.99;
          }
          double etaC = getZ(iz) + 0.5 * dz;
          transformToLab(etaC, vxC, vyC, vzC);  // viC is now in lab.frame!
          double gammaC = 1. / sqrt(1. - vxC * vxC - vyC * vyC - vzC * vzC);

          double uC[4] = {gammaC, gammaC * vxC, gammaC * vyC, gammaC * vzC};
          const double tauC = tau;
          double dsigma[4];
          // ---- transform dsigma to lab.frame :
          const double ch = cosh(etaC);
          const double sh = sinh(etaC);
          dsigma[0] = tauC * (ch * dx * dy * dz);
          dsigma[3] = tauC * (-sh * dx * dy * dz);
          dsigma[1] = 0.0;
          dsigma[2] = 0.0;
          double dVEff = 0.0;
          for (int ii = 0; ii < 4; ii++)
            dVEff += dsigma[ii] * uC[ii];  // normalize for Delta eta=1
          vEff += dVEff;
          for (int ii = 0; ii < 4; ii++) ffreeze << setw(24) << dsigma[ii];
          for (int ii = 0; ii < 4; ii++) ffreeze << setw(24) << uC[ii];
          ffreeze << setw(24) << TC << setw(24) << mubC << setw(24) << muqC
                  << setw(24) << musC;
          double picart[10];
          /*pi00*/ picart[index44(0, 0)] = ch * ch * piC[index44(0, 0)] +
                                           2. * ch * sh * piC[index44(0, 3)] +
                                           sh * sh * piC[index44(3, 3)];
          /*pi01*/ picart[index44(0, 1)] =
              ch * piC[index44(0, 1)] + sh * piC[index44(3, 1)];
          /*pi02*/ picart[index44(0, 2)] =
              ch * piC[index44(0, 2)] + sh * piC[index44(3, 2)];
          /*pi03*/ picart[index44(0, 3)] =
              ch * sh * (piC[index44(0, 0)] + piC[index44(3, 3)]) +
              (ch * ch + sh * sh) * piC[index44(0, 3)];
          /*pi11*/ picart[index44(1, 1)] = piC[index44(1, 1)];
          /*pi12*/ picart[index44(1, 2)] = piC[index44(1, 2)];
          /*pi13*/ picart[index44(1, 3)] =
              sh * piC[index44(0, 1)] + ch * piC[index44(3, 1)];
          /*pi22*/ picart[index44(2, 2)] = piC[index44(2, 2)];
          /*pi23*/ picart[index44(2, 3)] =
              sh * piC[index44(0, 2)] + ch * piC[index44(3, 2)];
          /*pi33*/ picart[index44(3, 3)] = sh * sh * piC[index44(0, 0)] +
                                           ch * ch * piC[index44(3, 3)] +
                                           2. * sh * ch * piC[index44(0, 3)];
#ifdef OUTPI
          for (int ii = 0; ii < 10; ii++) ffreeze << setw(24) << picart[ii];
          ffreeze << setw(24) << PiC;
#endif
          for(int i=0; i<4; i++)
          for(int j=0; j<4; j++)
           ffreeze << setw(24) << 0.0;
          ffreeze << endl;
          double dEsurfVisc = 0.;
          for (int i = 0; i < 4; i++)
            dEsurfVisc += picart[index44(0, i)] * dsigma[i];
          EtotSurf += (eC + pC) * uC[0] * dVEff - pC * dsigma[0] + dEsurfVisc;
          nbSurf += nbC * dVEff;
        }
        //----- end Cornelius
      }
  E = E * dx * dy * dz;
  Efull = Efull * dx * dy * dz;
  S = S * dx * dy * dz;
  Nb1 *= dx * dy * dz;
  Nb2 *= dx * dy * dz;
  fout_aniz << setw(12) << tau << setw(14) << vt_num / vt_den << setw(14)
            << vxvy_num / vxvy_den << setw(14) << pi0x_num / pi0x_den << endl;
  cout << setw(10) << "tau" << setw(13) << "E" << setw(13) << "Efull"
       << setw(13) << "Nb" << setw(13) << "Sfull" << setw(10) << "elements"
       << setw(10) << "susp." << setw(13) << "\%cut" << endl;
  cout << setw(10) << tau << setw(13) << E << setw(13) << Efull << setw(13)
       << nbSurf << setw(13) << S << endl;
#ifdef SWAP_EOS
  swap(eos, eosH);
#endif
  cout << "corona elements : " << nelements << endl;
}


void Fluid::outputOmegaXZ()
{
 double omegaXZ = 0.0;
 int ncell = 0;
 for (int ix = 2; ix < nx - 2; ix++)
  for (int iy = 2; iy < ny - 2; iy++)
   for (int iz = nz/2 - 1; iz < nz/2 + 2; iz++) {
    Cell *c = getCell(ix, iy, iz);
    double dOmegaXZ = 0.5* (c->getDbeta(1,3) - c->getDbeta(3,1));
    omegaXZ += dOmegaXZ;
    if(fabs(dOmegaXZ)>1e-10) ncell++;
   }
  cout << "initial_av_omegaXEta: " << omegaXZ/ncell << endl;
}
