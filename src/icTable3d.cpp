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
#include "icTable3d.h"
#include "rmn.h"
#include "s95p.h"

using namespace std;

IcTable3d::IcTable3d(Fluid* f, const char* filename, double _tau0)
 : ed_grid(nx_grid, std::vector<std::vector<double>>(ny_grid, std::vector<double>(neta_grid))) {
 cout << "loading 3D table from: " << filename << std::endl;

 tau0 = _tau0;
 ed_grid.resize(nx_grid);
 
 for (int ix = 0; ix < nx_grid; ix++)
  for (int iy = 0; iy < ny_grid; iy++) 
   for (int ieta = 0; ieta < neta_grid; ieta++) {
    ed_grid[ix][iy][ieta] = 0.0;
   }

 // ---- read the events
 ifstream fin(filename);
 if (!fin.good()) {
  cout << "I/O error with " << filename << endl;
  exit(1);
 }
 
 string line;
 istringstream instream;
  for (int i = 0; i < 1; i++) { // read and discard the first line
   getline(fin, line);
  }

 double x[nx_grid], y[ny_grid], eta[neta_grid];
 for (int ieta=0; ieta<neta_grid; ieta++)
  for (int iy = 0; iy < ny_grid; iy++)
    for (int ix = 0; ix < nx_grid; ix++) {
     fin >> x[ix] >> y[iy] >> eta[ieta] >> ed_grid[ix][iy][ieta];
  }

}

IcTable3d::~IcTable3d() {
}

double IcTable3d::interpolateGrid(double x, double y, double eta) {
 const double dxG = (xmaxG - xminG) / (nx_grid - 1);
 const double dyG = (ymaxG - yminG) / (ny_grid - 1);
 const double detaG = (etamaxG - etaminG) / (neta_grid - 1);
 int ix = (int)((x - xminG) / dxG);
 int iy = (int)((y - yminG) / dyG);
 int ieta = (int)((eta - etaminG) / detaG);
 if (ix < 0) ix = 0;
 if (iy < 0) iy = 0;
 if (ieta < 0) ieta = 0;
 if (ix > nx_grid - 2) ix = nx_grid - 2;
 if (iy > ny_grid - 2) iy = ny_grid - 2;
 if (ieta > neta_grid - 2) ieta = neta_grid - 2;
 const double xm = x - xminG - ix * dxG;
 const double ym = y - yminG - iy * dyG;
 const double etam = eta - etaminG - ieta * detaG;
 double wx[2] = {1. - xm / dxG, xm / dxG};
 double wy[2] = {1. - ym / dyG, ym / dyG};
 double weta[2] = {1. - etam / detaG, etam / detaG};
 double return_val = 0.;
 for (int jx = 0; jx < 2; jx++)
  for (int jy = 0; jy < 2; jy++)
   for (int jeta = 0; jeta < 2; jeta++){
    return_val += wx[jx] * wy[jy] * weta[jeta] * ed_grid[ix + jx][iy + jy][ieta + jeta];
  }
 return_val = std::max(return_val, 0.);
 return return_val;
}

void IcTable3d::setIC(Fluid* f, EoS* eos) {
 double E = 0.0, Px = 0.0, Py = 0.0, Pz = 0.0, Nb = 0.0, S = 0.0;
 double Jy0 = 0.0, Jint1 = 0.0, Jint3 = 0.0, Xcm = 0.0, Ycm = 0.0, Zcm = 0.0;
 double E_midrap = 0.0, Jy0_midrap = 0.0;  // same quantity at midrapidity
 double Tcm = 0.0;
 double e, p, nb;
 double total_energy = 0.0;
 for (int ix = 0; ix < f->getNX(); ix++)
  for (int iy = 0; iy < f->getNY(); iy++)
   for (int iz = 0; iz < f->getNZ(); iz++) {
    double x = f->getX(ix);
    double y = f->getY(iy);
    double eta = f->getZ(iz);
    e = interpolateGrid(x, y, eta);
    nb = 0.;
    p = eos->p(e, 0., 0., 0.);
    Cell* c = f->getCell(ix, iy, iz);
    const double ueta = 0.; //tanh(A*f->getX(ix))*sinh(ybeam-fabs(f->getZ(iz)));
    double u[4] = {sqrt(1.0+ueta*ueta), 0., 0., ueta};
    c->setPrimVar(eos, tau0, e, nb, 0.4*nb, 0., 0., 0., u[3]/u[0]);
    if (e > 0.) c->setAllM(1.);
    double coshEta = cosh(eta);
    double sinhEta = sinh(eta);
    double dxdydeta = f->getDx()*f->getDy()*f->getDz();
    total_energy += tau0*e*dxdydeta*coshEta;
    double u0lab = u[0] * coshEta + u[3] * sinhEta;
    double uzlab = u[0] * sinhEta + u[3] * coshEta;
    double dE = tau0 * ((e + p) * u[0] * u0lab - p * coshEta) * dxdydeta;
    E += dE;
    // if(zmin+iz*dz>0.)
    Pz += tau0 * ((e + p) * u[0] * uzlab - p * sinhEta) * dxdydeta;
    Px += tau0 * (e + p) * u[1] * u[0] * dxdydeta;
    Py += tau0 * (e + p) * u[2] * u[0] * dxdydeta;
    Nb += nb*tau0*dxdydeta;
    S += tau0 * eos->s(e, 0., 0., 0.) * u[0] * dxdydeta;
    // angular momentum calculation
    const double t = tau0 * coshEta;
    const double z = tau0 * sinhEta;
    Xcm += x * dE;
    Ycm += y * dE;
    Zcm += z * dE;
    Tcm += t * dE;
    Jy0 +=
        tau0 * (e + p) * u[0] * (z * u[1] - x * uzlab) * dxdydeta * gevtofm;
    Jint1 += tau0 * (e + p) * u[0] * u[1] * dxdydeta * gevtofm;
    Jint3 += tau0 * (e + p) * u[0] * uzlab * dxdydeta * gevtofm;
    if (iz > f->getNZ() / 2 - 2 && iz < f->getNZ() / 2 + 2) {
     E_midrap += dE;
     Jy0_midrap += tau0 * (e + p) * u[0] * (z * u[1] - x * uzlab) * dxdydeta * gevtofm;
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
 cout << "1/tau*dE/dy_ini: : " << E_midrap/(3.0*f->getDz()*tau0) << endl;
 cout << "1/tau*dJ/dy_ini: " << Jy0_midrap/(3.0*f->getDz()*tau0) << endl;
 //exit(1);
}

