#include <fstream>
#include <iostream>
#include <iomanip>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_miser.h>

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

const int nphi = 301;

// forward declarations
double WoodSaxon3D(double* x, size_t dim, void *params);
double WoodSaxon(double x, void *params);


ICGlauber::ICGlauber(double e, double impactPar, double _tau0) {
 epsilon = e;
 b = impactPar;
 tau0 = _tau0;
}

ICGlauber::~ICGlauber(void) {}

double ICGlauber::eProfile(double x, double y, WS_params* params) {
 gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
 gsl_function F;
 F.function = &WoodSaxon;
 double result, error;
 params->r = sqrt((x + b / 2.0) * (x + b / 2.0) + y * y);
 F.params = params;
 gsl_integration_qags(&F, -3.0 * Ra, 3.0 * Ra, 0, 1e-7, 1000, w, &result, &error);
 const double tpp = result;
 params->r = sqrt((x - b / 2.0) * (x - b / 2.0) + y * y);
 F.params = params;
 gsl_integration_qags(&F, -3.0 * Ra, 3.0 * Ra, 0, 1e-7, 1000, w, &result, &error);
 const double tmm = result;
 gsl_integration_workspace_free(w);
 return epsilon *
        pow(1. / rho0 * (tpp * (1.0 - pow((1.0 - sigma * tmm / A), A)) +
                         tmm * (1.0 - pow((1.0 - sigma * tpp / A), A))),
            1.0);
}

void ICGlauber::findRPhi(WS_params* params) {
 _rphi = new double[nphi];
 for (int iphi = 0; iphi < nphi; iphi++) {
  double phi = iphi * C_PI * 2. / (nphi - 1);
  double r = 0., r1 = 0., r2 = 2. * Ra;
  while (fabs((r2 - r1) / r2) > 0.001 && r2 > 0.001) {
   r = 0.5 * (r1 + r2);
   if (eProfile(r * cos(phi), r * sin(phi), params) > 0.5)
    r1 = r;
   else
    r2 = r;
  }
  _rphi[iphi] = r;
 }
}

double ICGlauber::rPhi(double phi) {
 const double cpi = C_PI;
 phi = phi - 2. * cpi * floor(phi / 2. / cpi);
 int iphi = (int)(phi / (2. * cpi) * (nphi - 1));
 int iphi1 = iphi + 1;
 if (iphi1 == nphi) iphi = nphi - 2;
 return _rphi[iphi] * (1. - (phi / (2. * cpi) * (nphi - 1) - iphi)) +
        _rphi[iphi1] * (phi / (2. * cpi) * (nphi - 1) - iphi);
}

void ICGlauber::setIC(Fluid *f, EoS *eos) {
 double e, nb, nq, vx = 0., vy = 0., vz = 0.;
 Cell *c;
 ofstream fvel("velocity_debug.txt");
 WS_params WSparams;
 WSparams.r = 0.;
 WSparams.norm = 1.0;
 WSparams.Ra = Ra;
 WSparams.delta = dlt;
 // === rewriting using Monte Carlo here
 double result, error;
 const gsl_rng_type *gsl_T;
 gsl_rng *gsl_r;
 gsl_monte_function gsl_integrand = { &WoodSaxon3D, 3, &WSparams};
 double xl[3] = {-3.0*Ra, -3.0*Ra, -3.0*Ra};
 double xu[3] = {3.0*Ra, 3.0*Ra, 3.0*Ra};
 size_t MC_calls = 500000;
 gsl_rng_env_setup ();
 gsl_T = gsl_rng_default;
 gsl_r = gsl_rng_alloc (gsl_T);
 gsl_monte_miser_state *gsl_miser_s = gsl_monte_miser_alloc(3);
 gsl_monte_miser_integrate(&gsl_integrand, xl, xu, 3, MC_calls, gsl_r, gsl_miser_s, &result, &error);
 double intgr2 = result;
 gsl_monte_miser_free(gsl_miser_s);
 if (intgr2 == 0.0) {
  std::cerr << "IC::setICGlauber Error! ff->Integral == 0; Return -1\n";
  exit(1);
 }
 std::cout << "a = " << A / intgr2 << endl;
 WSparams.r = 0.;
 WSparams.norm = A / intgr2;
 WSparams.Ra = Ra;
 WSparams.delta = dlt;
 // --- computing the norm === new code
 gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
 gsl_function F;
 F.function = &WoodSaxon;
 F.params = &WSparams;
 gsl_integration_qags(&F, -3.0 * Ra, 3.0 * Ra, 0, 1e-7, 1000, w, &result, &error);
 const double tpp = result;
 rho0 = 2.0 * tpp * (1.0 - pow((1.0 - sigma * tpp / A), A));
 gsl_integration_workspace_free(w);

 findRPhi(&WSparams);  // fill in R(phi) table
 cout << "R(phi) =  ";
 for (int jj = 0; jj < 5; jj++) cout << rPhi(jj * C_PI / 2.) << "  ";  // test
 cout << endl;
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
    double etaFactor;
    double eta1 = fabs(eta) < 1.3 ? 0.0 : fabs(eta) - 1.3;
    etaFactor = exp(-eta1 * eta1 / 2.1 / 2.1) * (fabs(eta) < 5.3 ? 1.0 : 0.0);
    e = eProfile(x, y, &WSparams) * etaFactor;
    if (e < 0.5) e = 0.0;
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
 cout << "total energy = "
      << Etotal * f->getDx() * f->getDy() * f->getDz() * tau0 << endl;
}

double WoodSaxon3D(double* x, size_t dim, void *params) {
 struct WS_params *p = (struct WS_params *)params;
 return p->norm / (exp((sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]) - p->Ra) / p->delta) + 1.0);
}

double WoodSaxon(double x, void *params) {
 struct WS_params *p = (struct WS_params *)params;
 return p->norm / (exp((sqrt(x*x + p->r*p->r) - p->Ra) / p->delta) + 1.0);
}
