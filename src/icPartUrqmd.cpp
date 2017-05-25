#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cfloat>

#include "eos.h"
#include "eoChiral.h"
#include "rmn.h"
#include "fld.h"
#include "icPartUrqmd.h"

using namespace std;

IcPartUrqmd::IcPartUrqmd(Fluid* f, char* filename, double _Rgt, double _Rgz,
                         double _tau0) {
  nx = f->getNX();
  ny = f->getNY();
  nz = f->getNZ();
  dx = f->getDx();
  dy = f->getDy();
  dz = f->getDz();
  xmin = f->getX(0);
  xmax = f->getX(nx - 1);
  ymin = f->getY(0);
  ymax = f->getY(ny - 1);
  zmin = f->getZ(0);
  zmax = f->getZ(nz - 1);

  tau0 = _tau0;
  Rgx = _Rgt;
  Rgy = _Rgt;
  Rgz = _Rgz;
  nsmoothx =
      (int)(6.0 * Rgx / dx);  // smoothly distribute to +- this many cells
  nsmoothy = (int)(6.0 * Rgy / dy);
  nsmoothz = (int)(1.5 * Rgz / dz);

  T00 = new double** [nx];
  T0x = new double** [nx];
  T0y = new double** [nx];
  T0z = new double** [nx];
  QB = new double** [nx];
  QE = new double** [nx];
  for (int ix = 0; ix < nx; ix++) {
    T00[ix] = new double* [ny];
    T0x[ix] = new double* [ny];
    T0y[ix] = new double* [ny];
    T0z[ix] = new double* [ny];
    QB[ix] = new double* [ny];
    QE[ix] = new double* [ny];
    for (int iy = 0; iy < ny; iy++) {
      T00[ix][iy] = new double[nz];
      T0x[ix][iy] = new double[nz];
      T0y[ix][iy] = new double[nz];
      T0z[ix][iy] = new double[nz];
      QB[ix][iy] = new double[nz];
      QE[ix][iy] = new double[nz];
      for (int iz = 0; iz < nz; iz++) {
        T00[ix][iy][iz] = 0.0;
        T0x[ix][iy][iz] = 0.0;
        T0y[ix][iy][iz] = 0.0;
        T0z[ix][iy][iz] = 0.0;
        QB[ix][iy][iz] = 0.0;
        QE[ix][iy][iz] = 0.0;
      }
    }
  }
#ifdef TSHIFT
  tau0 += tshift;
#endif
  // ---- read the events
  nevents = 0;
  ifstream fin(filename);
  if (!fin.good()) {
    cout << "I/O error with " << filename << endl;
    exit(1);
  }
  int np = 0;  // particle counter
  string line;
  istringstream instream;
  while (!fin.eof()) {
    getline(fin, line);
    instream.str(line);
    instream.seekg(0);
    instream.clear();
    instream >> Tau[np] >> X[np] >> Y[np] >> Eta[np] >> Mt[np] >> Px[np] >>
        Py[np] >> Rap[np] >> Id[np] >> Charge[np];
#ifdef TSHIFT
    Eta[np] = TMath::ATanH(Tau[np] * sinh(Eta[np]) /
                           (Tau[np] * cosh(Eta[np]) + tshift));
    Tau[np] += tshift;
#endif
    if (!instream.fail())
      np++;
    else if (np > 0) {
      // cout<<"readF14:instream: failure reading data\n" ;
      // cout<<"stream = "<<instream.str()<<endl ;
      if (nevents % 100 == 0) {
        cout << "event = " << nevents << "  np = " << np << "\r";
        cout << flush;
      }
      makeSmoothTable(np);
      np = 0;
      nevents++;
      // if(nevents>10000) return ;
    }
    if (np > NP - 1) cout << "ERROR: increase NP constant\n";
  }
  if (nevents > 1)
    cout << "++ Warning: loaded " << nevents << "  initial UrQMD events\n";
}

IcPartUrqmd::~IcPartUrqmd() {
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      delete[] T00[ix][iy];
      delete[] T0x[ix][iy];
      delete[] T0y[ix][iy];
      delete[] T0z[ix][iy];
      delete[] QB[ix][iy];
      delete[] QE[ix][iy];
    }
    delete[] T00[ix];
    delete[] T0x[ix];
    delete[] T0y[ix];
    delete[] T0z[ix];
    delete[] QB[ix];
    delete[] QE[ix];
  }
  delete[] T00;
  delete[] T0x;
  delete[] T0y;
  delete[] T0z;
  delete[] QB;
  delete[] QE;
}

void IcPartUrqmd::makeSmoothTable(int npart) {
  for (int ip = 0; ip < npart; ip++) {  // particle loop
    int ixc = (int)((X[ip] - xmin) / dx);
    int iyc = (int)((Y[ip] - ymin) / dy);
    int izc = (int)((Eta[ip] - zmin) / dz);
    // finding the norm
    double norm_gauss = 0.0;
    const double gammaz = cosh(Rap[ip] - Eta[ip]);
    for (int ix = ixc - nsmoothx; ix < ixc + nsmoothx + 1; ix++)
      for (int iy = iyc - nsmoothy; iy < iyc + nsmoothy + 1; iy++)
        for (int iz = izc - nsmoothz; iz < izc + nsmoothz + 1; iz++)
          if (ix > 0 && ix < nx && iy > 0 && iy < ny && iz > 0 && iz < nz) {
            const double xdiff = X[ip] - (xmin + ix * dx);
            const double ydiff = Y[ip] - (ymin + iy * dy);
            const double zdiff = Eta[ip] - (zmin + iz * dz);
            if (fabs(Rgz) < 1e-5)
              norm_gauss +=
                  exp(-xdiff * xdiff / Rgx / Rgx - ydiff * ydiff / Rgy / Rgy);
            else
              norm_gauss +=
                  exp(-xdiff * xdiff / Rgx / Rgx - ydiff * ydiff / Rgy / Rgy -
                      zdiff * zdiff / Rgz / Rgz * gammaz * gammaz * tau0 *
                          tau0);  // this term may become big, >800 !
          }

    for (int ix = ixc - nsmoothx; ix < ixc + nsmoothx + 1; ix++)
      for (int iy = iyc - nsmoothy; iy < iyc + nsmoothy + 1; iy++)
        for (int iz = izc - nsmoothz; iz < izc + nsmoothz + 1; iz++)
          if (ix > 0 && ix < nx && iy > 0 && iy < ny && iz > 0 && iz < nz) {
            const double xdiff = X[ip] - (xmin + ix * dx);
            const double ydiff = Y[ip] - (ymin + iy * dy);
            const double zdiff = Eta[ip] - (zmin + iz * dz);
            double weight;
            if (fabs(Rgz) < 1e-5)  // jun19
              weight = 1.0 / norm_gauss * exp(-xdiff * xdiff / Rgx / Rgx -
                                              ydiff * ydiff / Rgy / Rgy);
            else
              weight =
                  1.0 / norm_gauss *
                  exp(-xdiff * xdiff / Rgx / Rgx - ydiff * ydiff / Rgy / Rgy -
                      zdiff * zdiff / Rgz / Rgz * gammaz * gammaz * tau0 *
                          tau0);  // -zdiff*zdiff*...*tau*tau (Jun 16)
            if (weight != weight || fabs(weight) > DBL_MAX) {
              weight = 0.0;
            }
            T00[ix][iy][iz] +=
                weight * Mt[ip] * cosh(Rap[ip] - Eta[ip] + zdiff);
            T0x[ix][iy][iz] += weight * Px[ip];
            T0y[ix][iy][iz] += weight * Py[ip];
            T0z[ix][iy][iz] +=
                weight * Mt[ip] * sinh(Rap[ip] - Eta[ip] + zdiff);
            // if(ix==50 && iy==50 && iz==17){ // debug output
            // foutd<<setw(14)<<Mt[ip]<<setw(14)<<Rap[ip]<<setw(14)<<Eta[ip]-zdiff<<setw(14)<<weight<<endl
            // ;
            //}
            if (Id[ip] > 0 && Id[ip] < 56)
              QB[ix][iy][iz] += weight;  // correct expression (only weight)!
            if (Id[ip] < 0 && Id[ip] > -56) QB[ix][iy][iz] -= weight;
            QE[ix][iy][iz] += Charge[ip] * weight;
            //}
          }
  }  // end particle loop
}

void IcPartUrqmd::setIC(Fluid* f, EoS* eos) {
  double E = 0.0, Px = 0.0, Py = 0.0, Pz = 0.0, Nb = 0.0, S = 0.0;
  double Jx0 = 0., Jy0 = 0., Jint1 = 0., Jint2 = 0., Jint3 = 0.;
  double E_midrap = 0.0; // same quantity at midrapidity
  double Xcm = 0., Ycm = 0., Zcm = 0., Tcm = 0.0;
  double Q[7], e, p, nb, nq, ns, vx, vy, vz;
  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++)
      for (int iz = 0; iz < nz; iz++) {
        Q[T_] =
            T00[ix][iy][iz] / nevents / dx / dy / dz / tau0;  // /tau for Milne
        Q[X_] = T0x[ix][iy][iz] / nevents / dx / dy / dz / tau0;
        Q[Y_] = T0y[ix][iy][iz] / nevents / dx / dy / dz / tau0;
        Q[Z_] = T0z[ix][iy][iz] / nevents / dx / dy / dz / tau0;
        Q[NB_] = QB[ix][iy][iz] / nevents / dx / dy / dz / tau0;
        Q[NQ_] = QE[ix][iy][iz] / nevents / dx / dy / dz / tau0;
        Q[NS_] = 0.0;
        if (ix == nx / 2 && iy == ny / 2 && iz == nz / 2)
          cout << "IcUrqmd, center: " << xmin + ix* dx << "  " << zmin + iz* dz
               << "  " << Q[T_] << "  " << Q[Z_] << endl;
        transformPV(eos, Q, e, p, nb, nq, ns, vx, vy, vz);
        if (e < 1e-7 || fabs(f->getX(ix)) > 10. || fabs(f->getY(iy)) > 10. ||
            fabs(f->getZ(iz)) > 5.) {
          e = nb = nq = 0.0;
          vx = vy = vz = 0.0;
        }
        Cell* c = f->getCell(ix, iy, iz);
        c->setPrimVar(eos, tau0, e, nb, nq, ns, vx, vy, vz);
        if (e > 0.) c->setAllM(1.);
        const double gamma = 1.0 / sqrt(1.0 - vx * vx - vy * vy - vz * vz);
        double u[4] = {gamma, gamma * vx, gamma * vy, gamma * vz};
        double eta = zmin + iz * dz;
        double coshEta = cosh(eta);
        double sinhEta = sinh(eta);
        double u0lab = u[0] * coshEta + u[3] * sinhEta;
        double uzlab = u[0] * sinhEta + u[3] * coshEta;
        double dE = tau0 * ((e + p) * u[0] * u0lab - p * coshEta) * dx * dy * dz;
        E += dE;
        // if(zmin+iz*dz>0.)
        Pz += tau0 * ((e + p) * u[0] * uzlab - p * sinhEta) * dx * dy * dz;
        Px += tau0 * (e + p) * u[1] * u[0] * dx * dy * dz;
        Py += tau0 * (e + p) * u[2] * u[0] * dx * dy * dz;
        Nb += tau0 * nb * u[0] * dx * dy * dz;
        S += tau0 * eos->s(e, nb, nq, ns) * u[0] * dx * dy * dz;
        // angular momentum calculation
        const double t = tau0 * coshEta;
        const double z = tau0 * sinhEta;
        const double x = xmin + ix * dx;
        const double y = ymin + iy * dy;
        if(iz>nz/2 - 2 && iz<nz/2 + 2){
         Xcm += x * dE;
         Ycm += y * dE;
         Zcm += z * dE;
         Tcm += t * dE;
         Jx0 += tau0 * (e + p) * u[0] * (z * u[2] - y * uzlab) *
               dx * dy * dz * gevtofm;
         Jy0 += tau0 * (e + p) * u[0] * (z * u[1] - x * uzlab) *
               dx * dy * dz * gevtofm;
         Jint1 += tau0 * (e + p) * u[0] * u[1] * dx * dy * dz * gevtofm;
         Jint2 += tau0 * (e + p) * u[0] * u[2] * dx * dy * dz * gevtofm;
         Jint3 += tau0 * (e + p) * u[0] * uzlab * dx * dy * dz * gevtofm;
         E_midrap += dE;
	       }
      }
  Xcm = Xcm / E_midrap;
  Ycm = Ycm / E_midrap;
  Zcm = Zcm / E_midrap;
  Tcm = Tcm / E_midrap;
  double Jy = Jy0 - Zcm * Jint1 + Xcm * Jint3;
  double Jx = Jx0 - Zcm * Jint2 + Ycm * Jint3;
  cout << "hydrodynamic E = " << E << "  Pz = " << Pz << "  Nbar = " << Nb
       << endl << "  Px = " << Px << "  Py = " << Py << endl;
  cout << "initial_entropy S_ini = " << S << endl;
  cout << "Xcm: " << sqrt(Tcm*Tcm - Zcm*Zcm) << "  " << Xcm << "  " << Ycm <<
       "  " << 0.5*log((Tcm+Zcm)/(Tcm-Zcm)) << endl;
  cout << "initial_J_midrap  " << Jx0 << " " << Jy0 << endl;
  cout << "corrected_J_midrap  " << Jx << " " << Jy << endl;
  cout << "J_to_analyze " << setw(14) << E << setw(14) << Nb << setw(14) <<
       Jy << endl;
  cout << "midrapidity_E: " << E_midrap << endl;
//  exit(1);
}
