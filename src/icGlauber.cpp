#include <fstream>
#include <iomanip>
#include <iomanip>
#include <TF1.h>
#include <TF2.h>
#include <TGraph.h>

#include "fld.h"
#include "eos.h"
#include "icGlauber.h"
#include "inc.h"

using namespace std;

// --------------------------------------------
//   Initial state from optical Glauber model
// --------------------------------------------

// Au nucleus parameters for optical Glauber
const double A = 197.0;    // mass number
const double Ra = 6.37;    // radius
const double dlt = 0.54;   // diffuseness
const double sigma = 4.0;  // NN cross section in fm^2

const double etam = 2.0;
const double etaflat = 1.0;
const double sigEta = 1.3;

const int nphi = 301;

ICGlauber::ICGlauber(double e, double impactPar, double _tau0) {
  epsilon = e;
  b = impactPar;
  tau0 = _tau0;
}

ICGlauber::~ICGlauber(void) {}

double fplu(double eta)
{
 if(eta<-etam)
  return 0.;
 else if(eta<=etam)
  return (eta+etam)/(2.0*etam);
 else
  return 1.;
}

double fmin(double eta)
{
 if(eta<-etam)
  return 1.;
 else if(eta<=etam)
  return (-eta+etam)/(2.0*etam);
 else
  return 0.;
}

double ICGlauber::eProfile(double x, double y, double eta) {
 iff->SetParameters(sqrt((x + b / 2.0) * (x + b / 2.0) + y * y), 0.0); // second arg for the sake of syntax
 const double tpp = iff->Integral(-3.0 * Ra, 3.0 * Ra, 1.0e-9);
 iff->SetParameters(sqrt((x - b / 2.0) * (x - b / 2.0) + y * y), 0.0); // second arg for the sake of syntax
 const double tmm = iff->Integral(-3.0 * Ra, 3.0 * Ra, 1.0e-9);
 double T1 = tpp * (1.0 - pow((1.0 - sigma * tmm / A), A));
 double T2 = tmm * (1.0 - pow((1.0 - sigma * tpp / A), A));
 return 2.0 * (T1 * fmin(eta) + T2 * fplu(eta));
}


void ICGlauber::setIC(Fluid *f, EoS *eos) {
  double e, nb, nq, vx = 0., vy = 0., vz = 0.;
  Cell *c;
  ofstream fvel("velocity_debug.txt");

  TF2 *ff = 0;
  double prms2[2], intgr2;
  cout << "finding normalization constant\n";
  ff = new TF2("ThicknessF", this, &ICGlauber::Thickness, -3.0 * Ra, 3.0 * Ra,
               -3.0 * Ra, 3.0 * Ra, 2, "IC", "Thickness");
  prms2[0] = Ra;
  prms2[1] = dlt;
  ff->SetParameters(prms2);
  intgr2 = ff->Integral(-3.0 * Ra, 3.0 * Ra, -3.0 * Ra, 3.0 * Ra, 1.0e-9);
  if (intgr2 == 0.0) {
    cerr << "IC::setICGlauber Error! ff->Integral == 0; Return -1\n";
    delete ff;
    exit(1);
  }
  delete ff;
  cout << "a = " << A / intgr2 << endl;
  iff = new TF1("WoodSaxonDF", this, &ICGlauber::WoodSaxon, -3.0 * Ra, 3.0 * Ra, 4, "IC", "WoodSaxon");
  rho0 = eProfile(0., 0., 0.);

  //--------------
  double avv_num = 0., avv_den = 0.;
  double Etotal = 0.0;

  for (int ix = 0; ix < f->getNX(); ix++)
    for (int iy = 0; iy < f->getNY(); iy++)
      for (int iz = 0; iz < f->getNZ(); iz++) {
        c = f->getCell(ix, iy, iz);
        double x = f->getX(ix);
        double y = f->getY(iy);
        double eta = f->getZ(iz);
        double etat = fabs(eta) - 0.5 * etaflat;
        double H = etat>0. ? exp(-etat*etat/(2.0*sigEta*sigEta)) : 1.0;
        e = epsilon * eProfile(x, y, eta) / rho0 * H;
        if (e < 0.05) e = 0.0;
        vx = vy = 0.0;
        nb = nq = 0.0;
        vz = 0.0;

      avv_num += sqrt(vx * vx + vy * vy) * e;
      avv_den += e;

        c->setPrimVar(eos, tau0, e, nb, nq, 0., vx, vy, vz);
        double _p = eos->p(e, nb, nq, 0.);
        const double gamma2 = 1.0 / (1.0 - vx * vx - vy * vy - vz * vz);
        Etotal +=
            ((e + _p) * gamma2 * (cosh(eta) + vz * sinh(eta)) - _p * cosh(eta));
        c->saveQprev();

        if (e > 0.) c->setAllM(1.);
      }
  fvel.close();
  cout << "average initial flow = " << avv_num / avv_den << endl;
  cout << "total energy = " << Etotal *f->getDx() * f->getDy() * f->getDz() *
                                   tau0 << endl;
}

double ICGlauber::Thickness(double *x, double *p) {
  // p[0]: Ra radius; p[1]: delta = 0.54fm
  double intgrl;
  TF1 *iff = 0;
  iff = new TF1("WoodSaxonDF", this, &ICGlauber::WoodSaxon, -3.0 * p[0], 3.0 * p[0], 4,
                "IC", "WoodSaxon");
  iff->SetParameters( sqrt(x[0] * x[0] + x[1] * x[1]), 1.0, p[0], p[1]);
  intgrl = iff->Integral(-3.0 * p[0], 3.0 * p[0], 1.0e-9);
  if (iff) delete iff;
  return intgrl;
}

double ICGlauber::WoodSaxon(double *x, double *p) {
  return p[1] / (exp((sqrt(x[0] * x[0] + p[0] * p[0]) - p[2]) / p[3]) + 1.0);
}
