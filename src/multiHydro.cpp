#include <cmath>
#include <fstream>
#include <iostream>
#include <TMatrixDSymEigen.h>
#include <TMatrixDSym.h>

#include "multiHydro.h"
#include "hdo.h"
#include "fld.h"
#include "eos.h"
#include "rmn.h"
#include "trancoeff.h"
#include "cll.h"
#include "xsect.h"
#include "cornelius.h"

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
 double arrayDx[4] = {dtau, f_p->getDx(), f_p->getDy(), f_p->getDz()};
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

void MultiHydro::initOutput(const char *dir) {
 string outfreeze = dir;
 outfreeze.append("/freezeout.dat");
 fmhfreeze.open(outfreeze.c_str());
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
    for(int i=0; i<4; i++){
     flux_p[i] += -nbp*nbt*(D_P*(up[i] - ut[i]) + D_E*(up[i] + ut[i]))*h_p->getDtau();
     flux_t[i] += -nbp*nbt*(D_P*(ut[i] - up[i]) + D_E*(up[i] + ut[i]))*h_p->getDtau();
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
     fmhfreeze.precision(15);
     fmhfreeze << setw(24) << h_p->getTau() + cornelius->get_centroid_elem(isegm, 0)
               << setw(24) << f_p->getX(ix) + cornelius->get_centroid_elem(isegm, 1)
               << setw(24) << f_p->getY(iy) + cornelius->get_centroid_elem(isegm, 2)
               << setw(24) << f_p->getZ(iz) + cornelius->get_centroid_elem(isegm, 3);

     // interpolation procedure
     double vxC = 0., vyC = 0., vzC = 0., TC = 0., mubC = 0., muqC = 0.,
            musC = 0., piC[10], PiC = 0., nbC = 0., nqC = 0.;
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
     double wCenT[2] = {1. - cornelius->get_centroid_elem(isegm, 0) / dtau,
                        cornelius->get_centroid_elem(isegm, 0) / dtau};
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
       TmunuC[i][j] = TmunuC[i][j] / (h_p->getTau() + cornelius->get_centroid_elem(isegm, 0));
     for (int i = 0; i < 7; i++) {
      QC_p[i] = QC_p[i] / (h_p->getTau() + cornelius->get_centroid_elem(isegm, 0));
      QC_t[i] = QC_t[i] / (h_p->getTau() + cornelius->get_centroid_elem(isegm, 0));
      QC_f[i] = QC_f[i] / (h_p->getTau() + cornelius->get_centroid_elem(isegm, 0));
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

     // condition for switching the sign of resulting velocity
     if ((vzp+vzt+vzf) * vzC + (vxp+vxt+vxf) * vxC + (vyp+vyt+vyf) * vyC < 0) {
      vxC *= -1;
      vyC *= -1;
      vzC *= -1;
     }

     /*fmhfreeze << setw(24) << vxp << setw(24) << vyp << setw(24) << vzp
               << setw(24) << vxt << setw(24) << vyt << setw(24) << vzt
               << setw(24) << vxf << setw(24) << vyf << setw(24) << vzf
               << setw(24) << vxC << setw(24) << vyC << setw(24) << vzC
               << setw(24) << ep << setw(24) << et << setw(24) << ef << setw(24) << eC;*/

     /*transformPV(eos, QC, eC, pC, nbC, nqC, _ns, vxC, vyC, vzC);
     eos->eos(eC, nbC, nqC, _ns, TC, mubC, muqC, musC, pC);
     if (TC > 0.4 || fabs(mubC) > 0.85) {
      cout << "#### Error (surface): high T/mu_b (T=" << TC << "/mu_b=" << mubC << ") ####\n";
     }*/
     double v2C = vxC * vxC + vyC * vyC + vzC * vzC;
     if (v2C > 1.) {
      vxC *= sqrt(0.99 / v2C);
      vyC *= sqrt(0.99 / v2C);
      vzC *= sqrt(0.99 / v2C);
      v2C = 0.99;
     }
     double etaC = f_p->getZ(iz) + cornelius->get_centroid_elem(isegm, 3);
     transformToLab(etaC, vxC, vyC, vzC);  // viC is now in lab.frame!
     double gammaC = 1. / sqrt(1. - vxC * vxC - vyC * vyC - vzC * vzC);

     double uC[4] = {gammaC, gammaC * vxC, gammaC * vyC, gammaC * vzC};
     const double tauC = h_p->getTau() + cornelius->get_centroid_elem(isegm, 0);
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
     for (int ii = 0; ii < 4; ii++) fmhfreeze << setw(24) << dsigma[ii];
     for (int ii = 0; ii < 4; ii++) fmhfreeze << setw(24) << uC[ii];
     //fmhfreeze << setw(24) << TC << setw(24) << mubC << setw(24) << muqC
     //        << setw(24) << musC;*/
     fmhfreeze << endl;
    }
 }

 cout << setw(10) << h_p->getTau() << setw(10) << nelements << " " << ne_pos << endl;

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
}
