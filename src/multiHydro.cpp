#include <cmath>
#include <fstream>
#include <TMatrixDSymEigen.h>

#include "multiHydro.h"
#include "hdo.h"
#include "fld.h"
#include "eos.h"
#include "trancoeff.h"
#include "cll.h"
#include "xsect.h"


MultiHydro::MultiHydro(Fluid *_f_p, Fluid *_f_t, Fluid *_f_f, Hydro *_h_p,
 Hydro *_h_t, Hydro *_h_f, EoS *_eos, TransportCoeff *_trcoeff)
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
   c_p->addFlux(flux_p[0], flux_p[1], flux_p[2], flux_p[3], 0., 0., 0.);
   c_t->addFlux(flux_t[0], flux_t[1], flux_t[2], flux_t[3], 0., 0., 0.);
   c_f->addFlux(-flux_p[0]-flux_t[0], -flux_p[1]-flux_t[1],
    -flux_p[2]-flux_t[2], -flux_p[3]-flux_t[3], 0., 0., 0.);
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

void MultiHydro::findFreezeout()
{
 const double gmunu[4][4] = {
     {1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, -1, 0}, {0, 0, 0, -1}};
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
    double T [4][4] = {0.};
    for (int i=0; i<4; i++)
     for (int j=0; j<0; j++){
      T[i][j] = (ep + pp) * up[i] * up[j] - pp * gmunu[i][j]
       + (et + pt) * ut[i] * ut[j] - pt * gmunu[i][j]
       + (ef + pf) * uf[i] * uf[j] - pf * gmunu[i][j];
     }
   }
}
