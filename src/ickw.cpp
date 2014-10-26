#include <iostream>
#include <fstream>
#include <iomanip>
#include "ickw.h"
#include "eos.h"
#include "fld.h"
#include "inc.h"

using namespace std;

IC_KW::IC_KW(const char* filename) {
  ifstream finput(filename);
  if (!finput) {
    cout << "error! ICs file is not opened\n";
    exit(1);
  }

  double f1, f2, f3;

  nheadlines = 6;
  header = new char* [nheadlines];
  for (int i = 0; i < nheadlines; i++) {
    header[i] = new char[100];
    finput.getline(header[i], 100);
  }

  finput >> nx >> ny >> nz;
  e = new double** [nx];
  vx = new double** [nx];
  vy = new double** [nx];
  vz = new double** [nx];
  nu = new double** [nx];
  nd = new double** [nx];
  ns = new double** [nx];
  for (int ix = 0; ix < nx; ix++) {
    e[ix] = new double* [ny];
    vx[ix] = new double* [ny];
    vy[ix] = new double* [ny];
    vz[ix] = new double* [ny];
    nu[ix] = new double* [ny];
    nd[ix] = new double* [ny];
    ns[ix] = new double* [ny];
    for (int iy = 0; iy < ny; iy++) {
      e[ix][iy] = new double[nz];
      vx[ix][iy] = new double[nz];
      vy[ix][iy] = new double[nz];
      vz[ix][iy] = new double[nz];
      nu[ix][iy] = new double[nz];
      nd[ix][iy] = new double[nz];
      ns[ix][iy] = new double[nz];
    }
  }
  finput >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;
  dx = (xmax - xmin) / (nx);
  dy = (ymax - ymin) / (ny);
  dz = (zmax - zmin) / (nz);
  for (int iz = 0; iz < nz; iz++)
    for (int iy = 0; iy < ny; iy++)
      for (int ix = 0; ix < nx; ix++) {
        finput >> e[ix][iy][iz];
        if (e[ix][iy][iz] < 0.) e[ix][iy][iz] = 0.;
      }

  for (int iz = 0; iz < nz; iz++)
    for (int iy = 0; iy < ny; iy++)
      for (int ix = 0; ix < nx; ix++) {
        finput >> vx[ix][iy][iz] >> vy[ix][iy][iz] >> vz[ix][iy][iz];
        vz[ix][iy][iz] = vz[ix][iy][iz] * 0.6;
      }

  for (int iz = 0; iz < nz; iz++)
    for (int iy = 0; iy < ny; iy++)
      for (int ix = 0; ix < nx; ix++) {
        finput >> nu[ix][iy][iz] >> nd[ix][iy][iz] >> ns[ix][iy][iz];
      }

  // printf(" IC file : %s\n",filename);
  //	for(int ix=0; ix<nx; ix++)
  //		cout << setw(5) << ix << setw(15) << e[ix][ny/2][nz/2] << endl ;
  //	for(int iy=0; iy<ny; iy++)
  //		cout << setw(5) << iy << setw(15) << e[nx/2][iy][nz/2] << endl ;

  //	char a; cin>>a;
}

IC_KW::~IC_KW(void) {
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      delete[] e[ix][iy];
      delete[] vx[ix][iy];
      delete[] vy[ix][iy];
      delete[] vz[ix][iy];
      delete[] nu[ix][iy];
      delete[] nd[ix][iy];
      delete[] ns[ix][iy];
    }
    delete[] e[ix];
    delete[] vx[ix];
    delete[] vy[ix];
    delete[] vz[ix];
    delete[] nu[ix];
    delete[] nd[ix];
    delete[] ns[ix];
  }
  delete[] e;
  delete[] vx;
  delete[] vy;
  delete[] vz;
  delete[] nu;
  delete[] nd;
  delete[] ns;

  for (int i = 0; i < nheadlines; i++) delete[] header[i];
  delete[] header;
}

void IC_KW::writeHeader(ofstream& fout) {
  for (int i = 0; i < nheadlines; i++) fout << header[i] << endl;
}

void IC_KW::getICs(double x, double y, double z, double& _e, double& _vx,
                   double& _vy, double& _vz, double& _nu, double& _nd,
                   double& _ns) {
  //================= test ===========
  //	_e = (sqrt(x*x+y*y)<8.) ? 30.*exp(-(x*x+y*y)/5.4/5.4) : 0. ;
  //	_e = (fabs(x)<2. && fabs(y)<2.) ? 30. : 0. ;
  //	_nu = _nd = _ns = 0. ;
  //	_vx = _vy = 0. ;
  //	_vz = 0. ;
  //	return ;
  //==================================
  //	z = fabs(z) ;
  int ix = (x - xmin - dx / 2.) / dx;
  int iy = (y - ymin - dy / 2.) / dy;
  int iz = (z - zmin - dz / 2.) / dz;

  //	cout << "--------------IC_KW::setIC\n" ;
  //	cout << " x " << x << " y " << y << " z " << z << endl ;
  //	cout << "ix " << ix << " iy " << iy << " iz " << iz << endl ;

  if (ix < 0 || ix > nx - 2 || iy < 0 || iy > ny - 2 || iz < 0 || iz > nz - 2) {
    _e = _nu = _nd = _ns = _vx = _vy = _vz = 0.;
    return;
  }
  double xm = x - ix * dx - xmin - dx / 2.;
  double ym = y - iy * dy - ymin - dy / 2.;
  double zm = z - iz * dz - zmin - dz / 2.;

  //	for(int di=0; di<2; di++)
  //	for(int dj=0; dj<2; dj++)
  //	for(int dk=0; dk<2; dk++)
  //	 cout << "di " << di << " dj " << dj << " dk" << dk << " e " <<
  //vz[ix+di][iy+dj][iz+dk] << endl ;

  //	cout << " e = "	<< e[ix][iy][iz] << endl ;
  _e = 0.;
  _vx = 0.;
  _vy = 0.;
  _vz = 0., _nu = 0.;
  _nd = 0.;
  _ns = 0.;
  double wx[2] = {1. - xm / dx, xm / dx};
  double wy[2] = {1. - ym / dy, ym / dy};
  double wz[2] = {1. - zm / dz, zm / dz};
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        _e += wx[i] * wy[j] * wz[k] * e[ix + i][iy + j][iz + k];
        _vx += wx[i] * wy[j] * wz[k] * vx[ix + i][iy + j][iz + k];
        _vy += wx[i] * wy[j] * wz[k] * vy[ix + i][iy + j][iz + k];
        _vz += wx[i] * wy[j] * wz[k] * vz[ix + i][iy + j][iz + k];
        _nu += wx[i] * wy[j] * wz[k] * nu[ix + i][iy + j][iz + k];
        _nd += wx[i] * wy[j] * wz[k] * nd[ix + i][iy + j][iz + k];
        _ns += wx[i] * wy[j] * wz[k] * ns[ix + i][iy + j][iz + k];
      }
}

void IC_KW::setIC(Fluid* f, EoS* eos, double tau) {
  cout<<"IC_KW has to be modified for Cartesian; exit\n"; exit(1);
  for (int ix = 0; ix < f->getNX(); ix++)
    for (int iy = 0; iy < f->getNY(); iy++)
      for (int iz = 0; iz < f->getNZ(); iz++) {
        Cell* c = f->getCell(ix, iy, iz);
        double x = f->getX(ix);
        double y = f->getY(iy);
        double eta = f->getZ(iz);
        double e, vx, vy, vz, nu, nd, ns;
        getICs(x, y, eta, e, vx, vy, vz, nu, nd, ns);
        if (e < 0.) {
          e = nu = nd = ns = vx = vy = vz = 0.;
        }
        if (vx * vx + vy * vy + vz * vz > 1.) {
          //	cout << "setIC : " << ix << "  " << iy << "  " << iz << "  e = "
          //<< e << "  v^2 = " << vx*vx+vy*vy+vz*vz << endl ;
          double factor = sqrt(vx * vx + vy * vy + vz * vz);
          vx = vx * 0.99 / factor;
          vy = vy * 0.99 / factor;
          vz = vz * 0.99 / factor;
        }
        //============= test ==============
        //	vx = vy = 0. ;
        //=================================
        c->setPrimVar(eos, e, 0., 0., 0., vx, vy, vz);
        c->saveQprev();
        //=====viscous set-up: zero pi, Pi
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++) {
            c->setpi(i, j, 0.0);
            c->setpiH(i, j, 0.0);
            c->setpi0(i, j, 0.0);
            c->setpiH0(i, j, 0.0);
          }
        c->setPi(0.0);
        c->setPiH(0.0);
        c->setPi0(0.0);
        c->setPiH0(0.0);
        //===================
        if (e == 0.)
          c->setAllM(0.);
        else
          c->setAllM(1.);
        // test
        /*    double  ee, pp, nbx, nqx, nsx, vxx, vyx, vzx ;
            c->getPrimVar(eos, ee, pp, nbx, nqx, nsx, vxx, vyx, vzx) ;
            if( fabs(e-ee) > 0.01 ) {
            cout << "+++++++++++++++ ERROR in IC_KW::setIC"<< endl;
            cout << "+++++++++++++++ e_in  = "<< e << endl;
            cout << "+++++++++++++++ e_out = "<< ee << endl;
            cout <<ix<<" "<<iy<<" "<<iz<<" "<<e<<" "<<nb<<" "<<nq<<" "<<ns<<"
           "<<vx<<" "<<vy<<" "<<vz<<endl;
            cout <<ee<<" "<<pp<<" "<<nbx<<" "<<nqx<<" "<<nsx<<" "<<vxx<<"
           "<<vyx<<" "<<vzx<<endl;
            exit(1);
            }
        */
      }
}
