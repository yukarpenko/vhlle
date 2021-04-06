#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cfloat>

#include "fld.h"
#include "eos.h"
#include "ic3F.h"
#include "rmn.h"
#include "inc.h"

using namespace std;

IC3F::IC3F(Fluid *f_p, Fluid *f_t, double tau, int _nevents, double _snn, int _projA, int _targA, int _projZ, int _targZ, double _Rg) {
 // random generators
 random_device rd;
 mt19937 gen(rd());
 uniform_real_distribution<> dis(0.0, C_PI);
 uniform_real_distribution<> dis2(0.0, 2*C_PI);

 nx = f_p->getNX();
 ny = f_p->getNY();
 nz = f_p->getNZ();
 dx = f_p->getDx();
 dy = f_p->getDy();
 dz = f_p->getDz();
 xmin = f_p->getX(0);
 xmax = f_p->getX(nx - 1);
 ymin = f_p->getY(0);
 ymax = f_p->getY(ny - 1);
 zmin = f_p->getZ(0);
 zmax = f_p->getZ(nz - 1);

 snn = _snn;
 nevents = _nevents;
 projA = _projA;
 targA = _targA;
 projZ = _projZ;
 targZ = _targZ;

 const double vcoll = sqrt(1-pow(2*nucleon_mass/snn, 2)); // velocity of the beam
 const double gamma = 1/sqrt(1 - vcoll * vcoll); // gamma factor of nuclei
 const double rap_beam = 0.5 * log((1 + vcoll) / (1-vcoll)); // beam rapidity
 const double Rproj = 1.25 * pow(projA, 0.333); // projectile nucleus radius
 const double Rtarg = 1.25 * pow(targA, 0.333); // target nucleus radius
 const double z0_proj = - Rproj / gamma;
 const double z0_targ = Rtarg / gamma;
 tau0 = 2 * max(-z0_proj, z0_targ);
 Rgx = _Rg;
 Rgy = _Rg;
 Rgz = _Rg;
 nsmoothx = (int)(6.0 * Rgx / dx);
 nsmoothy = (int)(6.0 * Rgy / dy);
 nsmoothz = (int)(6.0 * Rgz / dz);

 // allocate tensor field for projectile
 T00_p = new double**[nx];
 T0z_p = new double**[nx];
 QB_p = new double**[nx];
 QE_p = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  T00_p[ix] = new double*[ny];
  T0z_p[ix] = new double*[ny];
  QB_p[ix] = new double*[ny];
  QE_p[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   T00_p[ix][iy] = new double[nz];
   T0z_p[ix][iy] = new double[nz];
   QB_p[ix][iy] = new double[nz];
   QE_p[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    T00_p[ix][iy][iz] = 0.0;
    T0z_p[ix][iy][iz] = 0.0;
    QB_p[ix][iy][iz] = 0.0;
    QE_p[ix][iy][iz] = 0.0;
   }
  }
 }

  // allocate tensor field for target
 T00_t = new double**[nx];
 T0z_t = new double**[nx];
 QB_t = new double**[nx];
 QE_t = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  T00_t[ix] = new double*[ny];
  T0z_t[ix] = new double*[ny];
  QB_t[ix] = new double*[ny];
  QE_t[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   T00_t[ix][iy] = new double[nz];
   T0z_t[ix][iy] = new double[nz];
   QB_t[ix][iy] = new double[nz];
   QE_t[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    T00_t[ix][iy][iz] = 0.0;
    T0z_t[ix][iy][iz] = 0.0;
    QB_t[ix][iy][iz] = 0.0;
    QE_t[ix][iy][iz] = 0.0;
   }
  }
 }

 // generating nucleons of projectile nucleus
 for (int i = 0; i < projA * nevents; i++) {
  double r = GenerateRadius(Rproj);
  double phi = dis2(gen);
  double theta = dis(gen);
  double x = r * sin(theta) * cos(phi);
  double y = r * sin(theta) * sin(phi);
  double z = (r * cos(theta)) / gamma - z0_proj;
  double eta = asinh(z * cosh(rap_beam) / tau0 - sinh(rap_beam)) + rap_beam;
  int charge = i < projZ ? 1 : 0;
  makeSmoothPart(x, y, eta, charge, rap_beam, true);
 }

 // generating nucleons of target nucleus
 for (int i = 0; i < targA * nevents; i++) {
  double r = GenerateRadius(Rtarg);
  double phi = dis2(gen);
  double theta = dis(gen);
  double x = r * sin(theta) * cos(phi);
  double y = r * sin(theta) * sin(phi);
  double z = (r * cos(theta)) / gamma - z0_targ;
  double eta = asinh(z * cosh(rap_beam) / tau0 + sinh(rap_beam)) - rap_beam;
  int charge = i < targZ ? 1 : 0;
  makeSmoothPart(x, y, eta, charge, rap_beam, false);
 }
}

IC3F::~IC3F() {
 for (int ix = 0; ix < nx; ix++) {
  for (int iy = 0; iy < ny; iy++) {
   delete[] T00_p[ix][iy];
   delete[] T0z_p[ix][iy];
   delete[] QB_p[ix][iy];
   delete[] QE_p[ix][iy];
   delete[] T00_t[ix][iy];
   delete[] T0z_t[ix][iy];
   delete[] QB_t[ix][iy];
   delete[] QE_t[ix][iy];
  }
  delete[] T00_p[ix];
  delete[] T0z_p[ix];
  delete[] QB_p[ix];
  delete[] QE_p[ix];
  delete[] T00_t[ix];
  delete[] T0z_t[ix];
  delete[] QB_t[ix];
  delete[] QE_t[ix];
 }
 delete[] T00_p;
 delete[] T0z_p;
 delete[] QB_p;
 delete[] QE_p;
 delete[] T00_t;
 delete[] T0z_t;
 delete[] QB_t;
 delete[] QE_t;
}

double IC3F::GenerateRadius(double R) {
 bool generated = false;
 random_device rd;
 mt19937 gen(rd());
 uniform_real_distribution<> dis(0.0, R + 4.0);
 uniform_real_distribution<> dis2(0.0, 1.0);
 double r;
 while (!generated) {
  r = dis(gen);
  if (1/(1 + exp((r - R) / WSdelta)) > dis2(gen)) generated = true;
 }
 return r;
}

void IC3F::makeSmoothPart(double x, double y, double eta, int Charge, double rap, bool isProjectile) {
 int ixc = (int)round((x - xmin) / dx);
 int iyc = (int)round((y - ymin) / dy);
 int izc = (int)round((eta - zmin) / dz);
 double m = nucleon_mass;
 double norm_gauss = 0.0;
 for (int ix = ixc - nsmoothx; ix < ixc + nsmoothx + 1; ix++)
  for (int iy = iyc - nsmoothy; iy < iyc + nsmoothy + 1; iy++)
   for (int iz = izc - nsmoothz; iz < izc + nsmoothz + 1; iz++)
    if (ix > 0 && ix < nx && iy > 0 && iy < ny && iz > 0 && iz < nz) {
     const double xdiff = x - (xmin + ix * dx);
     const double ydiff = y - (ymin + iy * dy);
     const double zdiff = eta - (zmin + iz * dz);
     norm_gauss +=
         exp(-xdiff * xdiff / 2 / Rgx / Rgx - ydiff * ydiff / 2 / Rgy / Rgy -
             zdiff * zdiff / 2 / Rgz / Rgz * tau0 * tau0 * cosh(eta) * cosh(eta) *
                 cosh(rap) * cosh(rap));
    }

 for (int ix = ixc - nsmoothx; ix < ixc + nsmoothx + 1; ix++)
  for (int iy = iyc - nsmoothy; iy < iyc + nsmoothy + 1; iy++)
   for (int iz = izc - nsmoothz; iz < izc + nsmoothz + 1; iz++)
    if (ix > 0 && ix < nx && iy > 0 && iy < ny && iz > 0 && iz < nz) {
     const double xdiff = x - (xmin + ix * dx);
     const double ydiff = y - (ymin + iy * dy);
     const double zdiff = eta - (zmin + iz * dz);
     double weight;
      weight = 1.0 / norm_gauss *
          exp(-xdiff * xdiff / 2 / Rgx / Rgx - ydiff * ydiff / 2 / Rgy / Rgy -
             zdiff * zdiff / 2 / Rgz / Rgz * tau0 * tau0 * cosh(eta) * cosh(eta) *
                 cosh(rap) * cosh(rap));
     if (weight != weight || fabs(weight) > DBL_MAX) {
      weight = 0.0;
     }
     if (isProjectile) {
      T00_p[ix][iy][iz] += weight * m * cosh(rap - eta + zdiff);
      T0z_p[ix][iy][iz] += weight * m * sinh(rap - eta + zdiff);
      QB_p[ix][iy][iz] += weight;
      QE_p[ix][iy][iz] += Charge * weight;
     } else {
      T00_t[ix][iy][iz] += weight * m * cosh(rap - eta + zdiff);
      T0z_t[ix][iy][iz] += weight * m * sinh(rap - eta + zdiff);
      QB_t[ix][iy][iz] += weight;
      QE_t[ix][iy][iz] += Charge * weight;
     }
     //}
    }
}

void IC3F::setIC(Fluid* f_p, Fluid* f_t, EoS* eos) {
 double E = 0.0, Px = 0.0, Py = 0.0, Pz = 0.0, Nb = 0.0, S = 0.0;
 double Q_p[7], e_p, p_p, nb_p, nq_p, ns_p, vx_p, vy_p, vz_p;
 double Q_t[7], e_t, p_t, nb_t, nq_t, ns_t, vx_t, vy_t, vz_t;
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    Q_p[T_] = T00_p[ix][iy][iz] / nevents / dx / dy / dz / tau0;  // /tau for Milne
    Q_p[X_] = 0.0;
    Q_p[Y_] = 0.0;
    Q_p[Z_] = T0z_p[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    Q_p[NB_] = QB_p[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    Q_p[NQ_] = QE_p[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    Q_p[NS_] = 0.0;

    Q_t[T_] = T00_t[ix][iy][iz] / nevents / dx / dy / dz / tau0;  // /tau for Milne
    Q_t[X_] = 0.0;
    Q_t[Y_] = 0.0;
    Q_t[Z_] = T0z_t[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    Q_t[NB_] = QB_t[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    Q_t[NQ_] = QE_t[ix][iy][iz] / nevents / dx / dy / dz / tau0;
    Q_t[NS_] = 0.0;

    transformPV(eos, Q_p, e_p, p_p, nb_p, nq_p, ns_p, vx_p, vy_p, vz_p);
    transformPV(eos, Q_t, e_t, p_t, nb_t, nq_t, ns_t, vx_t, vy_t, vz_t);
    if (e_p < 1e-7 || fabs(f_p->getX(ix)) > 10. || fabs(f_p->getY(iy)) > 10. ||
        fabs(f_p->getZ(iz)) > 5.) {
     e_p = nb_p = nq_p = 0.0;
     vx_p = vy_p = vz_p = 0.0;
    }
    if (e_t < 1e-7 || fabs(f_t->getX(ix)) > 10. || fabs(f_t->getY(iy)) > 10. ||
        fabs(f_t->getZ(iz)) > 5.) {
     e_t = nb_t = nq_t = 0.0;
     vx_t = vy_t = vz_t = 0.0;
    }
    Cell* c_p = f_p->getCell(ix, iy, iz);
    c_p->setPrimVar(eos, tau0, e_p, nb_p, nq_p, ns_p, vx_p, vy_p, vz_p);

    Cell* c_t = f_t->getCell(ix, iy, iz);
    c_t->setPrimVar(eos, tau0, e_t, nb_t, nq_t, ns_t, vx_t, vy_t, vz_t);
    if (e_p > 0.) c_p->setAllM(1.);
    if (e_t > 0.) c_t->setAllM(1.);
    /*const double gamma = 1.0 / sqrt(1.0 - vx * vx - vy * vy - vz * vz);
    double u[4] = {gamma, gamma * vx, gamma * vy, gamma * vz};
    double eta = zmin + iz * dz;
    E += tau0 * ((e + p) * u[0] * (u[0] * cosh(eta) + u[3] * sinh(eta)) -
                 p * cosh(eta)) *
         dx * dy * dz;
    // if(zmin+iz*dz>0.)
    Pz += tau0 * ((e + p) * u[0] * (u[0] * sinh(eta) + u[3] * cosh(eta)) -
                  p * sinh(eta)) *
          dx * dy * dz;
    Px += tau0 * (e + p) * u[1] * u[0] * dx * dy * dz;
    Py += tau0 * (e + p) * u[2] * u[0] * dx * dy * dz;
    Nb += tau0 * nb * u[0] * dx * dy * dz;
    S += tau0 * eos->s(e, nb, nq, ns) * u[0] * dx * dy * dz;

    if (iz == (int)nz/2 and iy == (int)ny/2) {
      cout << "x " << xmin + ix * dx << "  " << vx << "  " << e << endl;
    }
    if (iz == (int)nz/2 and ix == (int)ny/2) {
      cout << "y " << ymin + iy * dy << "  " << vy << "  " << e << endl;
    }
   }
 cout << "hydrodynamic E = " << E << "  Pz = " << Pz << "  Nbar = " << Nb
      << endl
      << "  Px = " << Px << "  Py = " << Py << endl;
 cout << "initial_entropy S_ini = " << S << endl;*/
   }
 exit(1);
}
