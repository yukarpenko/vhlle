#include <cmath>
#include <fstream>
#include <iostream>
#include <TMatrixDSymEigen.h>
#include <TMatrixDSym.h>
#include <TLorentzVector.h>

#include "multiHydro.h"
#include "hdo.h"
#include "fld.h"
#include "eos.h"
#include "rmn.h"
#include "trancoeff.h"
#include "cll.h"
#include "xsect.h"
#include "cornelius.h"

#define OUTPI

using namespace std;

MultiHydro::MultiHydro(Fluid *_f_p, Fluid *_f_t, Fluid *_f_f, Hydro *_h_p,
 Hydro *_h_t, Hydro *_h_f, EoS *_eos, TransportCoeff *_trcoeff, double _dtau,
 double eCrit)
{
 f_p = _f_p;
 f_t = _f_t;
 f_f = _f_f;
 h_p = _h_p;
 h_t = _h_t;
 h_f = _h_f;
 eos = _eos;
 trcoeff = _trcoeff;
 xsect = new CrossSections;
 nx = f_p->getNX();
 ny = f_p->getNY();
 nz = f_p->getNZ();
 dx = f_p->getDx();
 dy = f_p->getDy();
 dz = f_p->getDz();
 dtau = _dtau;

 //---- Cornelius init
 double arrayDx[4] = {h_p->getDtau(), f_p->getDx(), f_p->getDy(), f_p->getDz()};
 cornelius = new Cornelius;
 cornelius->init(4, eCrit, arrayDx);
 ecrit = eCrit;
 vEff = 0.;

 // allocate field for oveall energy density
 MHeps = new double**[nx];
 MHepsPrev = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  MHeps[ix] = new double*[ny];
  MHepsPrev[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   MHeps[ix][iy] = new double[nz];
   MHepsPrev[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    MHeps[ix][iy][iz] = 0.0;
    MHepsPrev[ix][iy][iz] = 0.0;
   }
  }
 }
}

MultiHydro::~MultiHydro() {
 for (int ix = 0; ix < nx; ix++) {
  for (int iy = 0; iy < ny; iy++) {
   delete[] MHeps[ix][iy];
   delete[] MHepsPrev[ix][iy];
  }
  delete[] MHeps[ix];
  delete[] MHepsPrev[ix];
 }
 delete[] MHeps;
 delete[] MHepsPrev;
}

void MultiHydro::setFluids(Fluid *_f_p, Fluid *_f_t, Fluid *_f_f, Hydro *_h_p,
 Hydro *_h_t, Hydro *_h_f) {
 f_p = _f_p;
 f_t = _f_t;
 f_f = _f_f;
 h_p = _h_p;
 h_t = _h_t;
 h_f = _h_f;
 nx = f_p->getNX();
 ny = f_p->getNY();
 nz = f_p->getNZ();
 dx = f_p->getDx();
 dy = f_p->getDy();
 dz = f_p->getDz();
 dtau = h_p->getDtau();

 //---- Cornelius init
 double arrayDx[4] = {h_p->getDtau(), f_p->getDx(), f_p->getDy(), f_p->getDz()};
 cornelius = new Cornelius;
 cornelius->init(4, ecrit, arrayDx);

 resizeMHeps();
}

void MultiHydro::initOutput(const char *dir) {
 string outfreeze_p = dir, outfreeze_f = dir, outfreeze_t = dir;
 outfreeze_p.append("/freezeout_p.dat");
 outfreeze_t.append("/freezeout_t.dat");
 outfreeze_f.append("/freezeout_f.dat");
 fmhfreeze_p.open(outfreeze_p.c_str());
 fmhfreeze_t.open(outfreeze_t.c_str());
 fmhfreeze_f.open(outfreeze_f.c_str());
}

void MultiHydro::performStep()
{
 h_p->performStep();
 h_t->performStep();
 h_f->performStep();
  /*if (trcoeff->isViscous()) {
   h_p->performViscSubstep();
   h_t->performViscSubstep();
   h_f->performViscSubstep();
  }*/
  frictionSubstep();
}

void MultiHydro::frictionSubstep()
{
 const double mN = 0.94; // nucleon mass [GeV]
 const double mpi = 0.1396; // pion mass [GeV]
 // here it is assumed that projectile and target grids
 // have same dimensions and physical sizes
 for (int iy = 0; iy < f_p->getNY(); iy++)
  for (int iz = 0; iz < f_p->getNZ(); iz++)
   for (int ix = 0; ix < f_p->getNX(); ix++) {
    Cell *c_p = f_p->getCell(ix, iy, iz);
    Cell *c_t = f_t->getCell(ix, iy, iz);
    Cell *c_f = f_f->getCell(ix, iy, iz);
    double ep, pp, nbp, nqp, nsp, vxp, vyp, vzp;
    double et, pt, nbt, nqt, nst, vxt, vyt, vzt;
    double ef, pf, nbf, nqf, nsf, vxf, vyf, vzf;
    c_p->getPrimVar(eos, h_p->getTau(), ep, pp, nbp, nqp, nsp, vxp, vyp, vzp);
    c_t->getPrimVar(eos, h_t->getTau(), et, pt, nbt, nqt, nst, vxt, vyt, vzt);
    c_f->getPrimVar(eos, h_f->getTau(), ef, pf, nbf, nqf, nsf, vxf, vyf, vzf);
    // 4-velocities, u_p and u_t
    double gammap = 1.0/sqrt(1.0-vxp*vxp-vyp*vyp-vzp*vzp);
    double up [4] = {gammap,gammap*vxp,gammap*vyp,gammap*vzp};
    double gammat = 1.0/sqrt(1.0-vxt*vxt-vyt*vyt-vzt*vzt);
    double ut [4] = {gammat,gammat*vxt,gammat*vyt,gammat*vzt};
    double gammaf = 1.0/sqrt(1.0-vxf*vxf-vyf*vyf-vzf*vzf);
    double uf [4] = {gammaf,gammaf*vxf,gammaf*vyf,gammaf*vzf};
    TLorentzVector upLV(up[1], up[2], up[3], up[0]);
    TLorentzVector utLV(ut[1], ut[2], ut[3], ut[0]);
    TLorentzVector ufLV(uf[1], uf[2], uf[3], uf[0]);
    double flux_p [4] = {0.}, flux_t [4] = {0.};
     // 1. projectile-target friction
    if (ep>0. && et>0.) {
    // u_p^\mu u_t_\mu
    double uput = gammap*gammat*(1.0 - vxp*vxt - vyp*vyt - vzp*vzt);
    // s (Mandelstam) variable
    double s = 2.0*mN*mN*(1.0 + uput);
    double Ekin = s/(2.0*mN) - 2.0*mN;
    double sigmaT, sigmaE, sigmaP;
    xsect->NN(Ekin, sigmaT, sigmaE, sigmaP);
    // Moeller factor
    double Vrel = sqrt(uput*uput - 1.0);
    // friction coefficient
    double D_P = mN*Vrel*sigmaP;
    double D_E = mN*Vrel*sigmaE; // SAME cross section moment for testing
    upLV.Boost(-vxt, -vyt, -vzt);
    utLV.Boost(-vxp, -vyp, -vzp);
    for(int i=0; i<4; i++){
     //flux_p[i] += -nbp*nbt*(D_P*(up[i] - ut[i]) + D_E*(up[i] + ut[i]))*h_p->getDtau();
     //flux_t[i] += -nbp*nbt*(D_P*(ut[i] - up[i]) + D_E*(up[i] + ut[i]))*h_p->getDtau();
     flux_p[i] += -upLV[(i+3)%4]*sqrt(ep*et)*h_p->getDtau()/lambda;
     flux_t[i] += -utLV[(i+3)%4]*sqrt(ep*et)*h_p->getDtau()/lambda;
    }
   }
   // 2. projectile-fireball friction
   if(ep>0. && ef>0.) {
    double s = mpi*mpi + mN*mN + 2.*mpi*mN*gammaf*gammap*
     (1.0 - vxf*vxp - vyf*vyp - vzf*vzp);
    double Vrel = 0.5/(mN*mpi)*sqrt(pow(s - mN*mN - mpi*mpi,2)
     - 4.*mN*mN*mpi*mpi);
    double D = Vrel*xsect->piN(s);
    for(int i=0; i<4; i++){
     flux_p[i] += D*nbp*(ef + pf)*uf[i]*h_p->getDtau();
    }
    flux_p[0] += -D*nbp*pf/uf[0]*h_p->getDtau();
   }
   if(et>0. && ef>0.) { // target-fireball friction
    double s = mpi*mpi + mN*mN + 2.*mpi*mN*gammaf*gammat*
     (1.0 - vxf*vxt - vyf*vyt - vzf*vzt);
    double Vrel = 0.5/(mN*mpi)*sqrt(pow(s - mN*mN - mpi*mpi,2)
     - 4.*mN*mN*mpi*mpi);
    double D = Vrel*xsect->piN(s);
    for(int i=0; i<4; i++){
     flux_t[i] += D*nbt*(ef + pf)*uf[i]*h_p->getDtau();
    }
    flux_t[0] += -D*nbt*pf/uf[0]*h_p->getDtau();
   }
   double taup = h_p->getTau();
   double taut = h_t->getTau();
   double tauf = h_f->getTau();
   double _Q_p[7], _Q_t[7], _Q_f[7];
   c_p->getQ(_Q_p);
   c_t->getQ(_Q_t);
   c_f->getQ(_Q_f);
   if (_Q_p[0] + flux_p[0]*taup >= 0.2*_Q_p[0] && _Q_t[0] + flux_t[0]*taut >= 0.2*_Q_t[0] && _Q_f[0] + (-flux_p[0]-flux_t[0])*tauf >= 0) {
    c_p->addFlux(flux_p[0]*taup, flux_p[1]*taup, flux_p[2]*taup, flux_p[3]*taup, 0., 0., 0.);
    c_t->addFlux(flux_t[0]*taut, flux_t[1]*taut, flux_t[2]*taut, flux_t[3]*taut, 0., 0., 0.);
    c_f->addFlux((-flux_p[0]-flux_t[0])*tauf, (-flux_p[1]-flux_t[1])*tauf,
     (-flux_p[2]-flux_t[2])*tauf, (-flux_p[3]-flux_t[3])*tauf, 0., 0., 0.);
    c_p->updateByFlux();
    c_t->updateByFlux();
    c_f->updateByFlux();
    c_p->clearFlux();
    c_t->clearFlux();
    c_f->clearFlux();
   }
   if(-flux_p[0]-flux_t[0] > 0. && c_f->getMaxM()<0.01)
    c_f->setAllM(1.0);
   } // end cell loop
}

void MultiHydro::getEnergyMomentumTensor(double (&T)[4][4], double Q_p[7], double Q_f[7], double Q_t[7])
{
 const double delta[4][4] = {
     {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
 double ep, pp, nbp, nqp, nsp, vxp, vyp, vzp;
 double et, pt, nbt, nqt, nst, vxt, vyt, vzt;
 double ef, pf, nbf, nqf, nsf, vxf, vyf, vzf;
 transformPV(eos, Q_p, ep, pp, nbp, nqp, nsp, vxp, vyp, vzp);
 transformPV(eos, Q_t, et, pt, nbt, nqt, nst, vxt, vyt, vzt);
 transformPV(eos, Q_f, ef, pf, nbf, nqf, nsf, vxf, vyf, vzf);
 // 4-velocities, u_p and u_t
 double gammap = 1.0/sqrt(1.0-vxp*vxp-vyp*vyp-vzp*vzp);
 double up [4] = {gammap,gammap*vxp,gammap*vyp,gammap*vzp};
 double gammat = 1.0/sqrt(1.0-vxt*vxt-vyt*vyt-vzt*vzt);
 double ut [4] = {gammat,gammat*vxt,gammat*vyt,gammat*vzt};
 double gammaf = 1.0/sqrt(1.0-vxf*vxf-vyf*vyf-vzf*vzf);
 double uf [4] = {gammaf,gammaf*vxf,gammaf*vyf,gammaf*vzf};

 // calculation of the energy-momentum tensor
 for (int i=0; i<4; i++){
  for (int j=0; j<4; j++){
   T[i][j] = (ep + pp) * up[i] * up[j] - pp * gmunu[i][j]
    + (et + pt) * ut[i] * ut[j] - pt * gmunu[i][j]
    + (ef + pf) * uf[i] * uf[j] - pf * gmunu[i][j];
  }
 }
}

void MultiHydro::getEnergyDensity()
{
 double Q_p[7], Q_f[7], Q_t[7];
 double Ttemp[4][4];
 for (int iy = 0; iy < f_p->getNY(); iy++)
  for (int iz = 0; iz < f_p->getNZ(); iz++)
   for (int ix = 0; ix < f_p->getNX(); ix++) {
    Cell *c_p = f_p->getCell(ix, iy, iz);
    Cell *c_t = f_t->getCell(ix, iy, iz);
    Cell *c_f = f_f->getCell(ix, iy, iz);
    c_p->getQ(Q_p);
    c_f->getQ(Q_f);
    c_t->getQ(Q_t);
    getEnergyMomentumTensor(Ttemp, Q_p, Q_f, Q_t);

    // calculation of the energy-momentum tensor
    TMatrixDSym T(4);
    for (int i=0; i<4; i++)
     for (int j=0; j<4; j++){
      T[i][j] = Ttemp[i][j]*gmunu[j][j];
    }
    // diagonalization of the energy-momentum tensor
    TMatrixDSymEigen Te(T);
    TVectorD eigenValues = Te.GetEigenValues();
    TMatrixD eigenVectors = Te.GetEigenVectors();

    double energyDensity;
    TVectorD v(4);
    for (int i=0; i<4; i++) {
     double vmuvmu = 0;
     energyDensity = eigenValues[i];
     v = TMatrixDColumn(eigenVectors,i);
     for (int j=0; j<4; j++) {
      vmuvmu += v[j]*v[j]*gmunu[j][j];
     }
     if (vmuvmu > 0 && energyDensity >= 0) {
      break;
     }
     else if (i == 3) {
      cout << "Multihydro: None of the eigenvectors is time-like, ";
      cout << "using largest eigenvalue for energy density." << endl;
      energyDensity = eigenValues[0];
      v = TMatrixDColumn(eigenVectors,0);
      break;
     }
    }

    // save computed energy density into private field
    MHeps[ix][iy][iz] = energyDensity;
   }
}

void MultiHydro::updateEnergyDensity()
{
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    MHepsPrev[ix][iy][iz] = MHeps[ix][iy][iz];
    MHeps[ix][iy][iz] = 0.0;
   }
}

void MultiHydro::resizeMHeps()
{
 double temp[nx][ny][nz];
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    if (2*(ix-(nx-1)/2)+(nx-1)/2 >= 0 && 2*(ix-(nx-1)/2)+(nx-1)/2 < nx && 2*(iy-(ny-1)/2)+(ny-1)/2 >=0 && 2*(iy-(ny-1)/2)+(ny-1)/2 < ny) {
     temp[ix][iy][iz] = MHeps[2*(ix-(nx-1)/2)+(nx-1)/2][2*(iy-(ny-1)/2)+(ny-1)/2][iz];
    } else {
     temp[ix][iy][iz] = 0.;
    }
   }
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    MHeps[ix][iy][iz] = temp[ix][iy][iz];
   }
}

void MultiHydro::findFreezeout()
{
 updateEnergyDensity();
 getEnergyDensity();

 int nelements = 0;
 int ne_pos = 0;
 double E=0., Efull = 0.;
 double Ttemp[4][4];

 // allocating corner points for Cornelius
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

 for (int ix = 2; ix < nx - 2; ix++)
  for (int iy = 2; iy < ny - 2; iy++)
   for (int iz = 2; iz < nz - 2; iz++) {
    double QCube_p[2][2][2][2][7], QCube_f[2][2][2][2][7], QCube_t[2][2][2][2][7];
    // array for storing full energy-momentum tensor of all three fluids at corners
    double TCube[2][2][2][2][4][4];
    double piSquare[2][2][2][10], PiSquare[2][2][2];

    // fill all corner cell with energy-momentum tensor
    for (int jx = 0; jx < 2; jx++)
     for (int jy = 0; jy < 2; jy++)
      for (int jz = 0; jz < 2; jz++) {
       ccube[0][jx][jy][jz] = MHepsPrev[ix + jx][iy + jy][iz + jz];
       ccube[1][jx][jy][jz] = MHeps[ix + jx][iy + jy][iz + jz];
       Cell *cc_p = f_p->getCell(ix + jx, iy + jy, iz + jz);
       Cell *cc_t = f_t->getCell(ix + jx, iy + jy, iz + jz);
       Cell *cc_f = f_f->getCell(ix + jx, iy + jy, iz + jz);
       double Qc_p[7], Qc_f[7], Qc_t[7];
       cc_p->getQ(QCube_p[1][jx][jy][jz]);
       cc_f->getQ(QCube_f[1][jx][jy][jz]);
       cc_t->getQ(QCube_t[1][jx][jy][jz]);
       getEnergyMomentumTensor(Ttemp, QCube_p[1][jx][jy][jz], QCube_f[1][jx][jy][jz], QCube_t[1][jx][jy][jz]);
       for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
         TCube[1][jx][jy][jz][i][j] = Ttemp[i][j];
       }
       cc_p->getQprev(QCube_p[0][jx][jy][jz]);
       cc_f->getQprev(QCube_f[0][jx][jy][jz]);
       cc_t->getQprev(QCube_t[0][jx][jy][jz]);
       getEnergyMomentumTensor(Ttemp, QCube_p[0][jx][jy][jz], QCube_f[0][jx][jy][jz], QCube_t[0][jx][jy][jz]);
       for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
         TCube[0][jx][jy][jz][i][j] = Ttemp[i][j];
       }
    }

    // cornelius
    cornelius->find_surface_4d(ccube);
    const int Nsegm = cornelius->get_Nelements();
    for (int isegm = 0; isegm < Nsegm; isegm++) {
     nelements++;
     fmhfreeze_p.precision(15);
     fmhfreeze_p << setw(24) << h_p->getTau() - h_p->getDtau() + cornelius->get_centroid_elem(isegm, 0)
               << setw(24) << f_p->getX(ix) + cornelius->get_centroid_elem(isegm, 1)
               << setw(24) << f_p->getY(iy) + cornelius->get_centroid_elem(isegm, 2)
               << setw(24) << f_p->getZ(iz) + cornelius->get_centroid_elem(isegm, 3);
     fmhfreeze_t.precision(15);
     fmhfreeze_t << setw(24) << h_p->getTau() - h_p->getDtau() + cornelius->get_centroid_elem(isegm, 0)
               << setw(24) << f_p->getX(ix) + cornelius->get_centroid_elem(isegm, 1)
               << setw(24) << f_p->getY(iy) + cornelius->get_centroid_elem(isegm, 2)
               << setw(24) << f_p->getZ(iz) + cornelius->get_centroid_elem(isegm, 3);
     fmhfreeze_f.precision(15);
     fmhfreeze_f << setw(24) << h_p->getTau() - h_p->getDtau() + cornelius->get_centroid_elem(isegm, 0)
               << setw(24) << f_p->getX(ix) + cornelius->get_centroid_elem(isegm, 1)
               << setw(24) << f_p->getY(iy) + cornelius->get_centroid_elem(isegm, 2)
               << setw(24) << f_p->getZ(iz) + cornelius->get_centroid_elem(isegm, 3);

     // interpolation procedure
     double vxC = 0., vyC = 0., vzC = 0., TC = 0., mubC = 0., muqC = 0.,
            musC = 0., piC[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            PiC = 0., nbC = 0., nqC = 0.;
     double QC_p[7] = {0., 0., 0., 0., 0., 0., 0.};
     double QC_t[7] = {0., 0., 0., 0., 0., 0., 0.};
     double QC_f[7] = {0., 0., 0., 0., 0., 0., 0.};
     double TmunuC[4][4];
     for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
       TmunuC[i][j] = 0;
     }
     double eC = 0., pC = 0.;
     for (int ii = 0; ii < 10; ii++) piC[ii] = 0.0;
     double wCenT[2] = {1. - cornelius->get_centroid_elem(isegm, 0) / h_p->getDtau(),
                        cornelius->get_centroid_elem(isegm, 0) / h_p->getDtau()};
     double wCenX[2] = {1. - cornelius->get_centroid_elem(isegm, 1) / dx,
                        cornelius->get_centroid_elem(isegm, 1) / dx};
     double wCenY[2] = {1. - cornelius->get_centroid_elem(isegm, 2) / dy,
                        cornelius->get_centroid_elem(isegm, 2) / dy};
     double wCenZ[2] = {1. - cornelius->get_centroid_elem(isegm, 3) / dz,
                        cornelius->get_centroid_elem(isegm, 3) / dz};

     for (int jt = 0; jt < 2; jt++)
      for (int jx = 0; jx < 2; jx++)
       for (int jy = 0; jy < 2; jy++)
        for (int jz = 0; jz < 2; jz++) {
         for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++) {
           TmunuC[i][j] += TCube[jt][jx][jy][jz][i][j] * wCenT[jt] * wCenX[jx] *
                   wCenY[jy] * wCenZ[jz];
         }
         for (int i = 0; i < 7; i++) {
          QC_p[i] += QCube_p[jt][jx][jy][jz][i] * wCenT[jt] * wCenX[jx] *
                   wCenY[jy] * wCenZ[jz];
          QC_t[i] += QCube_t[jt][jx][jy][jz][i] * wCenT[jt] * wCenX[jx] *
                   wCenY[jy] * wCenZ[jz];
          QC_f[i] += QCube_f[jt][jx][jy][jz][i] * wCenT[jt] * wCenX[jx] *
                   wCenY[jy] * wCenZ[jz];
         }
     }

     for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
       TmunuC[i][j] = TmunuC[i][j] / (h_p->getTau() - h_p->getDtau() + cornelius->get_centroid_elem(isegm, 0));
     for (int i = 0; i < 7; i++) {
      QC_p[i] = QC_p[i] / (h_p->getTau() - h_p->getDtau() + cornelius->get_centroid_elem(isegm, 0));
      QC_t[i] = QC_t[i] / (h_t->getTau() - h_t->getDtau() + cornelius->get_centroid_elem(isegm, 0));
      QC_f[i] = QC_f[i] / (h_f->getTau() - h_f->getDtau() + cornelius->get_centroid_elem(isegm, 0));
     }
     double _ns = 0.0;
     double ep, pp, nbp, nqp, nsp, vxp, vyp, vzp;
     double et, pt, nbt, nqt, nst, vxt, vyt, vzt;
     double ef, pf, nbf, nqf, nsf, vxf, vyf, vzf;
     transformPV(eos, QC_p, ep, pp, nbp, nqp, nsp, vxp, vyp, vzp);
     transformPV(eos, QC_t, et, pt, nbt, nqt, nst, vxt, vyt, vzt);
     transformPV(eos, QC_f, ef, pf, nbf, nqf, nsf, vxf, vyf, vzf);

     TMatrixDSym T(4);
     for (int i=0; i<4; i++)
      for (int j=0; j<4; j++){
       T[i][j] = TmunuC[i][j]*gmunu[j][j];
     }
     // diagonalization of the energy-momentum tensor
     TMatrixDSymEigen Te(T);
     TVectorD eigenValues = Te.GetEigenValues();
     TMatrixD eigenVectors = Te.GetEigenVectors();

     eC = eigenValues[0];
     TVectorD v(4);
     v = TMatrixDColumn(eigenVectors,0);
     vxC = v[1]/v[0];
     vyC = v[2]/v[0];
     vzC = v[3]/v[0];

     /*fmhfreeze << setw(24) << vxp << setw(24) << vyp << setw(24) << vzp
               << setw(24) << vxt << setw(24) << vyt << setw(24) << vzt
               << setw(24) << vxf << setw(24) << vyf << setw(24) << vzf
               << setw(24) << vxC << setw(24) << vyC << setw(24) << vzC
               << setw(24) << ep << setw(24) << et << setw(24) << ef << setw(24) << eC;*/

     //transformPV(eos, QC, eC, pC, nbC, nqC, _ns, vxC, vyC, vzC);
     nbC = nbp + nbt + nbf;
     nqC = nqp + nqt + nqf;
     _ns = nsp + nst + nsf;
     eos->eos(eC, nbC, nqC, _ns, TC, mubC, muqC, musC, pC);
     double TCp, mubCp, muqCp, musCp, pCp;
     double TCt, mubCt, muqCt, musCt, pCt;
     double TCf, mubCf, muqCf, musCf, pCf;
     eos->eos(ep, nbp, nqp, nsp, TCp, mubCp, muqCp, musCp, pCp);
     eos->eos(et, nbt, nqt, nst, TCt, mubCt, muqCt, musCt, pCt);
     eos->eos(ef, nbf, nqf, nsf, TCf, mubCf, muqCf, musCf, pCf);
     if (TC > 0.4 || fabs(mubC) > 0.85) {
      cout << "#### Error (multifluid surface): high T/mu_b (T=" << TC << "/mu_b=" << mubC << ") ####\n";
     }
     double v2C = vxC * vxC + vyC * vyC + vzC * vzC;
     if (v2C > 1.) {
      vxC *= sqrt(0.99 / v2C);
      vyC *= sqrt(0.99 / v2C);
      vzC *= sqrt(0.99 / v2C);
      v2C = 0.99;
     }
     double etaC = f_p->getZ(iz) + cornelius->get_centroid_elem(isegm, 3);
     transformToLab(etaC, vxC, vyC, vzC);  // viC is now in lab.frame!
     transformToLab(etaC, vxp, vyp, vzp);
     transformToLab(etaC, vxt, vyt, vzt);
     transformToLab(etaC, vxf, vyf, vzf);
     double gammaC = 1. / sqrt(1. - vxC * vxC - vyC * vyC - vzC * vzC);
     double gammaC_p = 1. / sqrt(1. - vxp * vxp - vyp * vyp - vzp * vzp);
     double gammaC_t = 1. / sqrt(1. - vxt * vxt - vyt * vyt - vzt * vzt);
     double gammaC_f = 1. / sqrt(1. - vxf * vxf - vyf * vyf - vzf * vzf);
     double uC[4] = {gammaC, gammaC * vxC, gammaC * vyC, gammaC * vzC};
     double uC_p[4] = {gammaC_p, gammaC_p * vxp, gammaC_p * vyp, gammaC_p * vzp};
     double uC_t[4] = {gammaC_t, gammaC_t * vxt, gammaC_t * vyt, gammaC_t * vzt};
     double uC_f[4] = {gammaC_f, gammaC_f * vxf, gammaC_f * vyf, gammaC_f * vzf};
     const double tauC = h_p->getTau() - h_p->getDtau() + cornelius->get_centroid_elem(isegm, 0);
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
     if (dVEff > 0) ne_pos++;
     vEff += dVEff;
     for (int ii = 0; ii < 4; ii++) {
      fmhfreeze_p << setw(24) << dsigma[ii];
      fmhfreeze_t << setw(24) << dsigma[ii];
      fmhfreeze_f << setw(24) << dsigma[ii];
     }
     for (int ii = 0; ii < 4; ii++) {
      fmhfreeze_p << setw(24) << uC_p[ii];
      fmhfreeze_t << setw(24) << uC_t[ii];
      fmhfreeze_f << setw(24) << uC_f[ii];
     }
     fmhfreeze_p << setw(24) << TCp << setw(24) << mubCp << setw(24) << muqCp
             << setw(24) << musCp;
     fmhfreeze_t << setw(24) << TCt << setw(24) << mubCt << setw(24) << muqCt
             << setw(24) << musCt;
     fmhfreeze_f << setw(24) << TCf << setw(24) << mubCf << setw(24) << muqCf
             << setw(24) << musCf;
#ifdef OUTPI
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
     for (int ii = 0; ii < 10; ii++) {
      fmhfreeze_p << setw(24) << picart[ii];
      fmhfreeze_t << setw(24) << picart[ii];
      fmhfreeze_f << setw(24) << picart[ii];
     }
     fmhfreeze_p << setw(24) << PiC;
     fmhfreeze_t << setw(24) << PiC;
     fmhfreeze_f << setw(24) << PiC;
#else
     fmhfreeze_p << setw(24) << dVEff;
     fmhfreeze_t << setw(24) << dVEff;
     fmhfreeze_f << setw(24) << dVEff;
#endif
     fmhfreeze_p << endl;
     fmhfreeze_t << endl;
     fmhfreeze_f << endl;
     double dEtotSurf[3] = {0., 0., 0.};
     dEtotSurf[0] = (ep + pCp) * uC_p[0] * dVEff - pCp * dsigma[0]; // projectile
     dEtotSurf[1] = (et + pCt) * uC_t[0] * dVEff - pCt * dsigma[0]; // target
     dEtotSurf[2] = (ef + pCf) * uC_f[0] * dVEff - pCf * dsigma[0]; // fireball
     EtotSurf[0] += dEtotSurf[0];
     EtotSurf[1] += dEtotSurf[1];
     EtotSurf[2] += dEtotSurf[2];
     if (dEtotSurf[0] > 0) EtotSurf_positive[0] += dEtotSurf[0];
     else EtotSurf_negative[0] += dEtotSurf[0];
     if (dEtotSurf[1] > 0) EtotSurf_positive[1] += dEtotSurf[1];
     else EtotSurf_negative[1] += dEtotSurf[1];
     if (dEtotSurf[2] > 0) EtotSurf_positive[2] += dEtotSurf[2];
     else EtotSurf_negative[2] += dEtotSurf[2];
    }
 }

 cout << setw(10) << h_p->getTau() << setw(10) << nelements << "\t" << ne_pos << "\t"
      << EtotSurf[0] << "\t" << EtotSurf_positive[0] << "\t" << EtotSurf_negative[0] << "\t"
      << EtotSurf[1] << "\t" << EtotSurf_positive[1] << "\t" << EtotSurf_negative[1] << "\t"
      << EtotSurf[2] << "\t" << EtotSurf_positive[2] << "\t" << EtotSurf_negative[2] << endl;

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
 if (nelements == 0) exit(0);
}
