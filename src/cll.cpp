#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "rmn.h"
#include "cll.h"

using namespace std;

// slope limiter; chooses minimal abs of the neighbouring slopes

double minmod(double a, double b) {
 if (a * b <= 0.) return 0.;
 //	else return (a*a*b+a*b*b)/(a*a+b*b) ;
 if (fabs(a) > fabs(b))
  return b;
 else
  return a;
}

// index44: returns an index of pi^{mu nu} mu,nu component in a plain 1D array
int index44(const int &i, const int &j) {
 if (i > 3 || j > 3 || i < 0 || j < 0) {
  std::cout << "index44: i j " << i << " " << j << endl;
  exit(1);
 }
 if (j < i)
  return (i * (i + 1)) / 2 + j;
 else
  return (j * (j + 1)) / 2 + i;
}

Cell::Cell() {
 for (int i = 0; i < 7; i++) {
  Q[i] = 0.;
  Qh[i] = 0.;
  Qprev[i] = 0.;
  flux[i] = 0.;
 }
 viscCorrCut = 1.;
 for (int i = 0; i < 10; i++) {
  pi[i] = 0.0;
  piH[i] = 0.0;
  pi0[i] = 0.0;
  piH0[i] = 0.0;
 }
 Pi = 0.0;
 PiH = 0.0;
 Pi0 = 0.0;
 PiH0 = 0.0;
 setAllM(0.);
}

void Cell::importVars(Cell* c) {
 c->getQ(Q);
 c->getQh(Qh);
 c->getQprev(Qprev); // not necessary but for consistency
 Pi = c->getPi();
 PiH = c->getPiH();
 for(int i=0; i<4; i++)
 for(int j=0; j<4; j++) {
  pi[index44(i, j)] = c->getpi(i, j);
  piH[index44(i, j)] = c->getpiH(i, j); // not necessary either
 }
 for(int i=X_; i<=Z_; i++)
  m[i-1] = c->getM(i);
 viscCorrCut = c->getViscCorrCutFlag();
}

void Cell::updateByFlux() {
 if(Q[0]+flux[0]<0.)
  return;
 for (int i = 0; i < 7; i++) Q[i] += flux[i];
}

void Cell::updateByViscFlux() {
 if(fabs(flux[0]) <= 0.5*Q[0]) {
  for (int i = 0; i < 7; i++) Q[i] += flux[i];
 } else if (flux[0]!=0.){
  double fac;
  fac = fabs(0.5*Q[0]/flux[0]);
  for (int i = 0; i < 7; i++) Q[i] += fac*flux[i];
 }
}

void Cell::updateQtoQhByFlux() {
 for (int i = 0; i < 7; i++) Qh[i] = Q[i] + flux[i];
}

void Cell::getPrimVar(EoS *eos, double tau, double &_e, double &_p, double &_nb,
                      double &_nq, double &_ns, double &_vx, double &_vy,
                      double &_vz) {
 double _Q[7];
 for (int i = 0; i < 7; i++) _Q[i] = Q[i] / tau;
 transformPV(eos, _Q, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
 //-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVar:\n";
  Dump(tau);
 }
 #endif
 //------------------------------------------
}

void Cell::getPrimVarLeft(EoS *eos, double tau, double &_e, double &_p,
                          double &_nb, double &_nq, double &_ns, double &_vx,
                          double &_vy, double &_vz, int dir) {
 double Qr[7], Ql[7], dQ[7];

 next[dir - 1]->getQ(Qr);
 prev[dir - 1]->getQ(Ql);

 for (int i = 0; i < 7; i++)
  dQ[i] = minmod((Qr[i] - Q[i]) / 2., (Q[i] - Ql[i]) / 2.);

 for (int i = 0; i < 7; i++) Ql[i] = (Q[i] - dQ[i]) / tau;
 transformPV(eos, Ql, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
 //-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVarLeft:\n";
  Dump(tau);
  next[dir - 1]->Dump(tau);
  prev[dir - 1]->Dump(tau);
 }
 #endif
 //------------------------------------------
}

void Cell::getPrimVarRight(EoS *eos, double tau, double &_e, double &_p,
                           double &_nb, double &_nq, double &_ns, double &_vx,
                           double &_vy, double &_vz, int dir) {
 double Qr[7], Ql[7], dQ[7];

 next[dir - 1]->getQ(Qr);
 prev[dir - 1]->getQ(Ql);

 for (int i = 0; i < 7; i++)
  dQ[i] = minmod((Qr[i] - Q[i]) / 2., (Q[i] - Ql[i]) / 2.);

 for (int i = 0; i < 7; i++) Qr[i] = (Q[i] + dQ[i]) / tau;
 transformPV(eos, Qr, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
 //-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVarRight:\n";
  Dump(tau);
  next[dir - 1]->Dump(tau);
  prev[dir - 1]->Dump(tau);
 }
 #endif
 //------------------------------------------
}

void Cell::getPrimVarHLeft(EoS *eos, double tau, double &_e, double &_p,
                           double &_nb, double &_nq, double &_ns, double &_vx,
                           double &_vy, double &_vz, int dir) {
 double Qr[7], Ql[7], dQ[7];

 next[dir - 1]->getQh(Qr);
 prev[dir - 1]->getQh(Ql);

 for (int i = 0; i < 7; i++)
  dQ[i] = minmod((Qr[i] - Qh[i]) / 2., (Qh[i] - Ql[i]) / 2.);

 for (int i = 0; i < 7; i++) Ql[i] = (Qh[i] - dQ[i]) / tau;
 transformPV(eos, Ql, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
 //-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVarHLeft:\n";
  Dump(tau);
  next[dir - 1]->Dump(tau);
  prev[dir - 1]->Dump(tau);
 }
 #endif
 //------------------------------------------
}

void Cell::getPrimVarHRight(EoS *eos, double tau, double &_e, double &_p,
                            double &_nb, double &_nq, double &_ns, double &_vx,
                            double &_vy, double &_vz, int dir) {
 double Qr[7], Ql[7], dQ[7];

 next[dir - 1]->getQh(Qr);
 prev[dir - 1]->getQh(Ql);

 for (int i = 0; i < 7; i++)
  dQ[i] = minmod((Qr[i] - Qh[i]) / 2., (Qh[i] - Ql[i]) / 2.);

 for (int i = 0; i < 7; i++) Qr[i] = (Qh[i] + dQ[i]) / tau;
 transformPV(eos, Qr, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
 //-------------------- debug ---------------
 #ifdef NAN_DEBUG
 if (std::isinf(_nb) or std::isnan(_nb)) {
  cout << "---error in getPrimVarHRight:\n";
  Dump(tau);
  next[dir - 1]->Dump(tau);
  prev[dir - 1]->Dump(tau);
 }
 #endif
 //------------------------------------------
}

void Cell::getPrimVarHCenter(EoS *eos, double tau, double &_e, double &_p,
                             double &_nb, double &_nq, double &_ns, double &_vx,
                             double &_vy, double &_vz) {
 double _Q[7];
 for (int i = 0; i < 7; i++) _Q[i] = Qh[i] / tau;
 transformPV(eos, _Q, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
}

void Cell::getPrimVarPrev(EoS *eos, double tau, double &_e, double &_p,
                          double &_nb, double &_nq, double &_ns, double &_vx,
                          double &_vy, double &_vz) {
 double _Q[7];
 for (int i = 0; i < 7; i++) _Q[i] = Qprev[i] / tau;
 transformPV(eos, _Q, _e, _p, _nb, _nq, _ns, _vx, _vy, _vz);
}

void Cell::setPrimVar(EoS *eos, double tau, double _e, double _nb, double _nq,
                      double _ns, double _vx, double _vy, double _vz) {
 const double gamma2 = 1. / (1 - _vx * _vx - _vy * _vy - _vz * _vz);
 const double p = eos->p(_e, _nb, _nq, _ns);
 Q[T_] = tau * (_e + p * (_vx * _vx + _vy * _vy + _vz * _vz)) * gamma2;
 Q[X_] = tau * (_e + p) * _vx * gamma2;
 Q[Y_] = tau * (_e + p) * _vy * gamma2;
 Q[Z_] = tau * (_e + p) * _vz * gamma2;
 Q[NB_] = tau * _nb * sqrt(gamma2);
 Q[NQ_] = tau * _nq * sqrt(gamma2);
 Q[NS_] = tau * _ns * sqrt(gamma2);
 if (std::isinf(Q[NB_]) or std::isnan(Q[NB_])) {
  cout << "init error!\n";
  eos->p(_e, _nb, _nq, _ns);
  cout << "e = " << _e << " p = " << p << " vx = " << _vx << " vy = " << _vy
       << " vz = " << _vz << " gamma2 = " << gamma2 << endl;
  //		exit(1) ;
  return;
 }
}

void Cell::Dump(double tau) {
 cout << "---------cell values dump-------\n";
 cout << setw(5) << ix << setw(5) << iy << setw(5) << iz << endl;
 cout << setw(14) << Q[0] / tau << setw(14) << Q[1] / tau << setw(14)
      << Q[2] / tau << setw(14) << Q[3] / tau << endl;
 cout << setw(14) << Q[4] / tau << setw(14) << Q[5] / tau << setw(14)
      << Q[6] / tau << endl;
 cout << setw(14) << Qh[0] / tau << setw(14) << Qh[1] / tau << setw(14)
      << Qh[2] / tau << setw(14) << Qh[3] / tau << endl;
 cout << setw(14) << Qh[4] / tau << setw(14) << Qh[5] / tau << setw(14)
      << Qh[6] / tau << endl;

 cout << "--------------------------------\n";
}
