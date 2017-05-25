#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

extern int glauberVariable;

namespace icurqmd {
using namespace std;

int NX, NY, NZ;  // N*N table
double ***e, ***nb, ***nq, ***vx, ***vy, ***vz, xmin, xmax, ymin, ymax, zmin,
    zmax;

void loadIC(char* filename) {
 ifstream fin(filename);
 if (!fin.good()) {
  cout << "I/O error with " << filename << endl;
  exit(1);
 }
 string str;
 fin >> str >> NX >> NY >> NZ;
 cout << "ICUrQMD: NX " << NX << " NY " << NY << " NZ " << NZ << endl;
 double* x = new double[NX];
 double* y = new double[NY];
 double* z = new double[NZ];
 e = new double**[NX];
 nb = new double**[NX];
 nq = new double**[NX];
 vx = new double**[NX];
 vy = new double**[NX];
 vz = new double**[NX];
 for (int ix = 0; ix < NX; ix++) {
  e[ix] = new double*[NY];
  nb[ix] = new double*[NY];
  nq[ix] = new double*[NY];
  vx[ix] = new double*[NY];
  vy[ix] = new double*[NY];
  vz[ix] = new double*[NY];
  for (int iy = 0; iy < NY; iy++) {
   e[ix][iy] = new double[NZ];
   nb[ix][iy] = new double[NZ];
   nq[ix][iy] = new double[NZ];
   vx[ix][iy] = new double[NZ];
   vy[ix][iy] = new double[NZ];
   vz[ix][iy] = new double[NZ];
  }
 }

 for (int ix = 0; ix < NX; ix++)
  for (int iy = 0; iy < NY; iy++)
   for (int iz = 0; iz < NZ; iz++) {
    fin >> x[ix] >> y[iy] >> z[iz] >> e[ix][iy][iz] >> nb[ix][iy][iz] >>
        nq[ix][iy][iz] >> vx[ix][iy][iz] >> vy[ix][iy][iz] >> vz[ix][iy][iz];
   }
 xmin = x[0];
 xmax = x[NX - 1];
 ymin = y[0];
 ymax = y[NY - 1];
 zmin = z[0];
 zmax = z[NZ - 1];
 cout << "ICUrQMD: table read, limits are " << xmin << "  " << xmax << "  "
      << ymin << "  " << ymax << "  " << zmin << "  " << zmax << endl;
 delete[] x;
 delete[] y;
 delete[] z;
}

void getIC(double x, double y, double eta, double& eout, double& nbout,
           double& nqout, double& vxout, double& vyout, double& vzout) {
 const double dx = (xmax - xmin) / (NX - 1);
 const double dy = (ymax - ymin) / (NY - 1);
 const double dz = (zmax - zmin) / (NZ - 1);
 const int ix = (int)((x - xmin) / dx);
 const int iy = (int)((y - ymin) / dy);
 const int iz = (int)((eta - zmin) / dz);
 if (ix < 0 || ix > NX - 2 || iy < 0 || iy > NY - 2 || iz < 0 || iz > NZ - 2) {
  eout = vxout = vyout = vzout = 0.0;
  return;
 }
 const double sx = x - xmin - ix * dx;
 const double sy = y - ymin - iy * dy;
 const double sz = eta - zmin - iz * dz;
 double wx[2] = {(1. - sx / dx), sx / dx};
 double wy[2] = {(1. - sy / dy), sy / dy};
 double wz[2] = {(1. - sz / dz), sz / dz};
 eout = nbout = nqout = vxout = vyout = vzout = 0.0;
 for (int jx = 0; jx < 2; jx++)
  for (int jy = 0; jy < 2; jy++)
   for (int jz = 0; jz < 2; jz++) {
    eout += wx[jx] * wy[jy] * wz[jz] * e[ix + jx][iy + jy][iz + jz];
    nbout += wx[jx] * wy[jy] * wz[jz] * nb[ix + jx][iy + jy][iz + jz];
    nqout += wx[jx] * wy[jy] * wz[jz] * nq[ix + jx][iy + jy][iz + jz];
    vxout += wx[jx] * wy[jy] * wz[jz] * vx[ix + jx][iy + jy][iz + jz];
    vyout += wx[jx] * wy[jy] * wz[jz] * vy[ix + jx][iy + jy][iz + jz];
    vzout += wx[jx] * wy[jy] * wz[jz] * vz[ix + jx][iy + jy][iz + jz];
   }
}

}  // end icurqmd
