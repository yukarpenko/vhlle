#include <cfloat>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstring>

#include "eos.h"
#include "eoChiral.h"
#include "fld.h"
#include "icGlissando.h"
#include "rmn.h"
#include "s95p.h"

using namespace std;

IcGlissando::IcGlissando(Fluid* f, const char* filename, double _tau0, const char* setup) {
 cout << "loading GLISSANDO IC\n";
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
 
 if(strcmp(setup,"LHC276")==0) {
  eta0 = 2.3; // midrapidity plateau
  sigEta = 1.4; // diffuseness of rapidity profile
  ybeam = 7.98; // beam rapidity, for 200 GeV RHIC
  alphaMix = 0.15; // WN/binary mixing
  Rg = 0.4; // Gaussian smearing in transverse dir
  sNorm = 0.96; // normalization of initial entropy profile
  A = 3.6e-5; // initial shear flow
  cout << "IcGlissando: setup for 2.76 TeV LHC\n";
 } else if(strcmp(setup,"RHIC200")==0) {
  eta0 = 1.5; // midrapidity plateau
  sigEta = 1.4; // diffuseness of rapidity profile
  ybeam = 5.36; // beam rapidity, for 200 GeV RHIC
  alphaMix = 0.125; // WN/binary mixing
  Rg = 0.4; // Gaussian smearing in transverse dir
  sNorm = 0.56; // normalization of initial entropy profile
  A = 5e-4; // initial shear flow
  cout << "IcGlissando: setup for 200 GeV RHIC\n";
 } else {
  cout << "IcGlissando: optional parameter LHC276 or RHIC200 is expected\n";
  exit(0);
 }

 nsmoothx = (int)(3.0 * Rg / dx);  // smoothly distribute to +- this many cells
 nsmoothy = nsmoothx;

 rho = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  rho[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   rho[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    rho[ix][iy][iz] = 0.0;
   }
  }
 }
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
  instream >> X[np] >> Y[np] >> C[np] >> W[np];
  if (!instream.fail() && C[np]!=0) // binary scattering entries are ignored!
   np++;
  if (instream.fail() && np > 0) {
   // cout<<"readF14:instream: failure reading data\n" ;
   // cout<<"stream = "<<instream.str()<<endl ;
   if (nevents % 1 == 0) {
    cout << "event = " << nevents << "  np = " << np << "\n";
    //cout << flush;
   }
   makeSmoothTable(np);
   np = 0;
   nevents++;
   //if(nevents>2) return ;
  }
  if (np > NP - 1) cout << "ERROR: increase NP constant\n";
 }
 if (nevents > 1)
  cout << "IcGlissando: loaded " << nevents << "  initial UrQMD events\n";
}

IcGlissando::~IcGlissando() {
 for (int ix = 0; ix < nx; ix++) {
  for (int iy = 0; iy < ny; iy++) {
   delete[] rho[ix][iy];
  }
  delete[] rho[ix];
 }
 delete[] rho;
}

void IcGlissando::makeSmoothTable(int npart) {
 for (int ip = 0; ip < npart; ip++) {  // particle loop
  int ixc = (int)((X[ip] - xmin) / dx);
  int iyc = (int)((Y[ip] - ymin) / dy);
  // finding the norm
  const double norm_gauss = 2.0*C_PI*Rg*Rg;

  for (int ix = ixc - nsmoothx; ix < ixc + nsmoothx + 1; ix++)
   for (int iy = iyc - nsmoothy; iy < iyc + nsmoothy + 1; iy++)
    if (ix > 0 && ix < nx && iy > 0 && iy < ny) {
     const double xdiff = X[ip] - (xmin + ix * dx);
     const double ydiff = Y[ip] - (ymin + iy * dy);
     double trSmear;
     trSmear = 1.0 / norm_gauss *
              exp(-xdiff * xdiff / (2. * Rg * Rg) - ydiff * ydiff / (2. * Rg * Rg));
     if (trSmear != trSmear || fabs(trSmear) > DBL_MAX) {
      trSmear = 0.0;
     }
     for (int iz = 0; iz < nz; iz++) {
      // longidudinal profile here
      const double eta = zmin + iz * dz;
      const double tilt = C[ip]>0 ? 1. + eta/ybeam : 1. - eta/ybeam;
      double fEta = 0.;
      if(fabs(eta)<eta0) fEta = 1.0;
      else if (fabs(eta)<ybeam) fEta = exp(-0.5*pow((fabs(eta)-eta0)/sigEta,2));
      rho[ix][iy][iz] += trSmear * fEta * tilt
       * ((1.0 - alphaMix) + abs(C[ip])*alphaMix);
     } // Z(eta) loop
    }
 }  // end particle loop
}

void IcGlissando::setIC(Fluid* f, EoS* eos) {
 double E = 0.0, Px = 0.0, Py = 0.0, Pz = 0.0, Nb = 0.0, S = 0.0;
 double Jy0 = 0.0, Jint1 = 0.0, Jint3 = 0.0, Xcm = 0.0, Ycm = 0.0, Zcm = 0.0;
 double E_midrap = 0.0, Jy0_midrap = 0.0;  // same quantity at midrapidity
 double Tcm = 0.0;
 double e, p;
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    e = s95p::s95p_e(sNorm * rho[ix][iy][iz] / nevents / dx / dy);
    p = eos->p(e, 0., 0., 0.);
    Cell* c = f->getCell(ix, iy, iz);
    const double ueta = tanh(A*f->getX(ix))*sinh(ybeam-fabs(f->getZ(iz)));
    double u[4] = {sqrt(1.0+ueta*ueta), 0., 0., ueta};
    c->setPrimVar(eos, tau0, e, 0., 0., 0., 0., 0., u[3]/u[0]);
    if (e > 0.) c->setAllM(1.);
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
    Nb += 0.;
    S += tau0 * eos->s(e, 0., 0., 0.) * u[0] * dx * dy * dz;
    // angular momentum calculation
    const double t = tau0 * coshEta;
    const double z = tau0 * sinhEta;
    const double x = xmin + ix * dx;
    const double y = ymin + iy * dy;
    Xcm += x * dE;
    Ycm += y * dE;
    Zcm += z * dE;
    Tcm += t * dE;
    Jy0 +=
        tau0 * (e + p) * u[0] * (z * u[1] - x * uzlab) * dx * dy * dz * gevtofm;
    Jint1 += tau0 * (e + p) * u[0] * u[1] * dx * dy * dz * gevtofm;
    Jint3 += tau0 * (e + p) * u[0] * uzlab * dx * dy * dz * gevtofm;
    if (iz > nz / 2 - 2 && iz < nz / 2 + 2) {
     E_midrap += dE;
     Jy0_midrap += tau0 * (e + p) * u[0] * (z * u[1] - x * uzlab) * dx * dy *
                   dz * gevtofm;
    }
   }
 Xcm = Xcm / E;
 Ycm = Ycm / E;
 Zcm = Zcm / E;
 Tcm = Tcm / E;
 double Jy = Jy0 - Zcm * Jint1 + Xcm * Jint3;
 cout << "hydrodynamic E = " << E << "  Pz = " << Pz << "  Nbar = " << Nb
      << endl
      << "  Px = " << Px << "  Py = " << Py << endl;
 cout << "initial_entropy S_ini = " << S << endl;
 cout << "Xcm: " << sqrt(Tcm * Tcm - Zcm * Zcm) << "  " << Xcm << "  " << Ycm
      << "  " << 0.5 * log((Tcm + Zcm) / (Tcm - Zcm)) << endl;
 cout << "initial/corrected J_y  " << Jy0 << " " << Jy << endl;
 cout << "J_to_analyze " << setw(14) << E << setw(14) << Nb << setw(14) << Jy
      << endl;
 cout << "1/tau*dE/dy_ini: : " << E_midrap/(3.0*dz*tau0) << endl;
 cout << "1/tau*dJ/dy_ini: " << Jy0_midrap/(3.0*dz*tau0) << endl;
 //exit(1);
}
