#include <fstream>
#include <iomanip>
#include <iostream>

#include "fld.h"
#include "eos.h"
#include "icBjorken.h"
#include "inc.h"

using namespace std;

ICBjorken::ICBjorken(double e) : epsilon(e) {
}

ICBjorken::~ICBjorken(void) {}

void ICBjorken::setIC(Fluid *f, EoS *eos, double tau) {
  double e, nb, nq, vx = 0., vy = 0., vz = 0.;
  Cell *c;

  double avv_num = 0., avv_den = 0.;
  double Etotal = 0.0;

  for (int ix = 0; ix < f->getNX(); ix++)
    for (int iy = 0; iy < f->getNY(); iy++)
      for (int iz = 0; iz < f->getNZ(); iz++) {
        c = f->getCell(ix, iy, iz);
        double x = f->getX(ix);
        double y = f->getY(iy);
        double eta = f->getZ(iz);

		e = epsilon;
		vx = vy = 0.;

        nb = nq = 0.0;
        vz = 0.0;

        avv_num += sqrt(vx * vx + vy * vy) * e;
        avv_den += e;

        c->setPrimVar(eos, tau, e, nb, nq, 0., vx, vy, vz);
        double _p = eos->p(e, nb, nq, 0., c->getTauP());
        const double gamma2 = 1.0 / (1.0 - vx * vx - vy * vy - vz * vz);
        Etotal +=
            ((e + _p) * gamma2 * (cosh(eta) + vz * sinh(eta)) - _p * cosh(eta));
        c->saveQprev();

        if (e > 0.) c->setAllM(1.);
      }
  cout << "average initial flow = " << avv_num / avv_den << endl;
  cout << "total energy = " << Etotal *f->getDx() * f->getDy() * f->getDz() *
                                   tau << endl;
}
