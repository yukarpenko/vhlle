#include <iostream>
#include <fstream>
#include <iomanip>

#include "fld.h"
#include "eos.h"
#include "ic.h"
#include "inc.h"

using namespace std;


IC::IC(char* icInputFile, double s0ScaleFactor)
{
}

IC::~IC(void) {}


void IC::setIC(Fluid *f, EoS *eos, double tau) {
  double e, nb, nq, vx = 0., vy = 0., vz = 0.;
  Cell *c;

  double avv_num = 0., avv_den = 0.;
  double Etotal = 0.0, Stotal = 0.0;

  for (int ix = 0; ix < f->getNX(); ix++)
    for (int iy = 0; iy < f->getNY(); iy++)
      for (int iz = 0; iz < f->getNZ(); iz++) {
        c = f->getCell(ix, iy, iz);
        double x = f->getX(ix);
        double y = f->getY(iy);
        double z = f->getZ(iz);




         e=30.*exp( - x*x -y*y - z*z ) ;
         nb = nq = 0.0;
         if(e<0.01) e=0. ;
         vx = vy = vz = 0. ;
        avv_num += sqrt(vx * vx + vy * vy) * e;
        avv_den += e;

        c->setPrimVar(eos, e, nb, nq, 0., vx, vy, vz);
        double _p = eos->p(e, nb, nq, 0.);
        double _s = eos->s(e, nb, nq, 0.);
        const double gamma2 = 1.0 / (1.0 - vx * vx - vy * vy - vz * vz);
        Etotal +=
            ((e + _p) * gamma2 - _p );
        c->saveQprev();
        Stotal += _s * sqrt(gamma2);

        if (e > 0.) c->setAllM(1.);
      }
  cout << "average initial flow = " << avv_num / avv_den << endl;
  cout << "total energy = " << Etotal *f->getDx() * f->getDy() * f->getDz() *
                                   tau << endl;
}
