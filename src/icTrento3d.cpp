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
#include "icTrento3d.h"
#include "rmn.h"
#include "s95p.h"

using namespace std;

IcTrento3d::IcTrento3d(Fluid* f, const char* filename, double _tau0, const char* setup) {
 cout << "loading TRENTO3d IC\n";
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
 double norm=2.5;

if(strcmp(setup,"LHC276")==0) {
  sNN = 2760;
  cout << "IcTrento3d: setup for 2.76 TeV LHC\n";
 } else if(strcmp(setup,"RHIC200")==0) {
  sNN = 200;
  cout << "IcTrento3d: setup for 200 GeV RHIC\n";
 } else if(strcmp(setup,"LHC5020")==0) {
  sNN = 5020;
  cout << "IcTrento3d: setup for 5.02 TeV LHC\n";
 } else if(strcmp(setup,"RHIC62")==0) {
  sNN = 62.4;
  cout << "IcTrento3d: setup for 62.4 GeV RHIC\n";
 } else if(strcmp(setup,"RHIC27")==0) {
  sNN = 27;
  cout << "IcTrento3d: setup for 27 GeV RHIC\n";
 } else {
  cout << "IcTrento3d: optional parameter LHC276 or RHIC200 is expected\n";
  exit(0);
 }

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
 int npart = 0;  // number of participants
 string line;
 istringstream instream;
 while (!fin.eof()) {
  for (int i = 0; i < 3; i++) { // read first 3 lines - the 3rd line contains npart
   getline(fin, line);
  }
  if (line == "") break;
  instream.str(line);
  instream.seekg(9);
  instream >> npart; // read npart
  for (int i = 0; i < 5; i++) { // read the rest 5 lines from header
   getline(fin, line);
  }
  getline(fin, line); // read grid-step
  getline(fin, line); // read grid-nsteps
  instream.str(line);
  instream.seekg(15);
  instream >> n_grid;

  getline(fin, line); // read grid-max
  instream.str(line);
  instream.seekg(12);
  instream >> xmaxG;

  getline(fin, line); 
  getline(fin, line); // read grid-nsteps-eta
  instream.str(line);
  instream.seekg(19);
  instream >> n_grid_eta;

  getline(fin, line); // read grid-max-eta
  instream.str(line);
  instream.seekg(16);
  instream >> etamaxG;
 
    // allocate source array 3D - [x,y,eta]
  source = new double**[n_grid];
  for (int ix = 0; ix < n_grid; ix++) {
      source[ix] = new double*[n_grid];
    for (int iy = 0; iy < n_grid; iy++) {
      source[ix][iy] = new double[n_grid_eta];
      for (int eta=0;eta<n_grid_eta;eta++){
        source[ix][iy][eta] = 0.0;
      }  
   }
  }

  ymaxG = xmaxG;
  xminG = -xmaxG;
  yminG = -xmaxG;
  etaminG=-etamaxG;
  
  if (nevents == 0) cout << "Trento IS grid: x,y,etamaxG = " << xmaxG <<" "<< etamaxG << "  n_grid = " << n_grid << " n_grid_eta= "<<n_grid_eta<< endl;
  for (int iy = 0; iy < n_grid; iy++) {
   for (int ix = 0; ix < n_grid; ix++) {
    for (int eta=0;eta<n_grid_eta;eta++){
      fin >> source[ix][iy][eta];
    }
   }
  }

  nevents += 1;
  cout << "event " << nevents << "  npart = " << npart << "\n";
  makeSmoothTable(npart);

  // delete source array
  for (int ix = 0; ix < n_grid; ix++) {
   delete[] source[ix];
  }
  delete[] source;
  getline(fin, line);
 }

 // autocalculation of sNorm and nNorm
 sNorm = 1.0;
 double old_sNorm = 0.0;
 do {
   old_sNorm = sNorm;
   sNorm = pow(setNormalization(npart), 0.75)*old_sNorm;
 } while (abs(sNorm-old_sNorm) > 0.0001);
 cout << "sNorm set to " << sNorm << endl;
 double old_nNorm = 0.0;

}

IcTrento3d::~IcTrento3d() {
 for (int ix = 0; ix < nx; ix++) {
  for (int iy = 0; iy < ny; iy++) {
   delete[] rho[ix][iy];
   
  }
  delete[] rho[ix];
  
 }
 delete[] rho;
 
}

double IcTrento3d::interpolateGrid(double x, double y,double eta) {
 const double dxG = (xmaxG - xminG) / (n_grid - 1);
 const double dyG = (ymaxG - yminG) / (n_grid - 1);
 const double detaG = (etamaxG - etaminG) / (n_grid_eta - 1);

 int ix = (int)((x - xminG) / dxG);
 int iy = (int)((y - yminG) / dyG);
 int ieta = (int)((eta - etaminG) / detaG);

 if (ix < 0) ix = 0;
 if (iy < 0) iy = 0;
 if (ieta < 0) ieta = 0;

 if (ix > n_grid - 2) ix = n_grid - 2;
 if (iy > n_grid - 2) iy = n_grid - 2;
 if (ieta > n_grid_eta - 2) ieta = n_grid_eta - 2;
 
 const double xm = x - xminG - ix * dxG;
 const double ym = y - yminG - iy * dyG;
 const double etam = eta - etaminG - ieta * detaG;
 
 double wx[2] = {1. - xm / dxG, xm / dxG};
 double wy[2] = {1. - ym / dyG, ym / dyG};
 double weta[2] = {1. - etam / detaG, etam / detaG};
 
 double return_val = 0.;
 for (int jx = 0; jx < 2; jx++){
  for (int jy = 0; jy < 2; jy++) {
    for (int jeta = 0; jeta < 2; jeta++) {
   return_val += wx[jx] * wy[jy] * weta[jeta] *source[ix + jx][iy + jy][ieta + jeta];
           }
      }
  }
 return_val = std::max(return_val, 0.);
 return return_val;
}

void IcTrento3d::makeSmoothTable(int npart) {
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    // longidudinal profile here
    const double x = xmin + ix * dx;
    const double y = ymin + iy * dy;
    const double eta = zmin + iz * dz;
    double baryonGaussian;
    
    rho[ix][iy][iz] += interpolateGrid(x,y,eta);
 } // Z(eta) loop
}

void IcTrento3d::setIC(Fluid* f, EoS* eos) {
 double E = 0.0, Px = 0.0, Py = 0.0, Pz = 0.0, Nb = 0.0, S = 0.0;
 double Jy0 = 0.0, Jint1 = 0.0, Jint3 = 0.0, Xcm = 0.0, Ycm = 0.0, Zcm = 0.0;
 double E_midrap = 0.0, Jy0_midrap = 0.0;  // same quantity at midrapidity
 double Tcm = 0.0;
 double e, p, nb;
 double total_energy = 0.0;
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    e = norm*(s95p::s95p_e(sNorm * rho[ix][iy][iz] / nevents));
    p = eos->p(e, 0., 0., 0.);
    Cell* c = f->getCell(ix, iy, iz);
    const double ueta = tanh(A*f->getX(ix))*sinh(ybeam-fabs(f->getZ(iz)));
    double u[4] = {sqrt(1.0+ueta*ueta), 0., 0., ueta};
    c->setPrimVar(eos, tau0, e, nb, 0.4*nb, 0., 0., 0., u[3]/u[0]);
    if (e > 0.) c->setAllM(1.);
    double eta = zmin + iz * dz;
    double coshEta = cosh(eta);
    double sinhEta = sinh(eta);
    total_energy += tau0*e*dx*dy*dz*coshEta;
    double u0lab = u[0] * coshEta + u[3] * sinhEta;
    double uzlab = u[0] * sinhEta + u[3] * coshEta;
    double dE = tau0 * ((e + p) * u[0] * u0lab - p * coshEta) * dx * dy * dz;
    E += dE;
    // if(zmin+iz*dz>0.)
    Pz += tau0 * ((e + p) * u[0] * uzlab - p * sinhEta) * dx * dy * dz;
    Px += tau0 * (e + p) * u[1] * u[0] * dx * dy * dz;
    Py += tau0 * (e + p) * u[2] * u[0] * dx * dy * dz;
    Nb += nb*tau0*dx*dy*dz;
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
 
}

double IcTrento3d::setNormalization(int npart) {
double e;
double total_energy = 0.0;
for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    e = s95p::s95p_e(sNorm * rho[ix][iy][iz] / nevents);
    double eta = zmin + iz * dz;
    double coshEta = cosh(eta);
    total_energy += tau0*e*dx*dy*dz*coshEta;
   }
 return npart*0.5*sNN/total_energy;
}

