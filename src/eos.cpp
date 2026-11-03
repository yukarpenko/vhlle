#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "eos.h"
#include "inc.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

const double bagp = pow(247.19 / 197.32, 4) / gevtofm;
const double bagt = pow(247.19 / 197.32, 4) / gevtofm;
const double deg = 16.0 + 3.0 * 12.0 * (7.0 / 8.0);

// EoS choise
#define TABLE  // Laine, etc
//#define SIMPLE  // p=e/3

double EoS::s(double e, double nb, double nq, double ns) {
 double T, mub, muq, mus, p;
 eos(e, nb, nq, ns, T, mub, muq, mus, p);
 if (T > 0.0)
  return (e + p - mub * nb - muq * nq - mus * ns) / T;
 else
  return 0.;
}

EoSs::EoSs(string fname, int ncols) {
#if defined TABLE || defined LAINE_CFO

 int edat = 10000;
 double* e = new double[edat];
 double* pGrid = new double[edat];
 double* tpGrid = new double[edat];
 double* muGrid = new double[edat];

 ifstream finput(fname.c_str(), ios::in);
 if (!finput) {
  cerr << "can't open input file \"" << fname.c_str() << "\"" << endl;
  exit(1);
 }
 edat = 0;
 while (!finput.eof()) {
  if (ncols == 3) {
   finput >> e[edat] >> pGrid[edat] >> tpGrid[edat];
   muGrid[edat] = 0.;
  } else {
   finput >> e[edat] >> pGrid[edat] >> tpGrid[edat] >> muGrid[edat];
  }
  if (pGrid[edat] < 0.) pGrid[edat] = 0.;
  // monotonicity check
  if(edat>0 and e[edat]<=e[edat-1]) {
   cout << "WARNING non-monotonic EoS table " << edat << "  " << e[edat-1] << "  " << e[edat] << endl;
  }
  edat++;
 }
 edat--;
 finput.close();

 // === initialising GSL interpolation machinery
 const gsl_interp_type *Tint = gsl_interp_cspline;
#ifdef _OPENMP
 nAcc = omp_get_max_threads();
#else
 nAcc = 1;
#endif
 acc_p.resize(nAcc);
 acc_T.resize(nAcc);
 acc_mu.resize(nAcc);
 for (int i = 0; i < nAcc; i++) {
  acc_p[i] = gsl_interp_accel_alloc();
  acc_T[i] = gsl_interp_accel_alloc();
  acc_mu[i] = gsl_interp_accel_alloc();
 }
 spline_p = gsl_spline_alloc(Tint, edat);
 gsl_spline_init(spline_p, e, pGrid, edat);
 spline_T = gsl_spline_alloc(Tint, edat);
 gsl_spline_init(spline_T, e, tpGrid, edat);
 spline_mu = gsl_spline_alloc(Tint, edat);
 gsl_spline_init(spline_mu, e, muGrid, edat);

#elif defined SIMPLE
// nothing
#endif
}

EoSs::~EoSs(void) {
#if defined TABLE
 gsl_spline_free(spline_p);
 gsl_spline_free(spline_T);
 gsl_spline_free(spline_mu);
 for (int i = 0; i < nAcc; i++) {
  gsl_interp_accel_free(acc_p[i]);
  gsl_interp_accel_free(acc_T[i]);
  gsl_interp_accel_free(acc_mu[i]);
 }
#endif
}

// index of the accelerator owned by the calling thread
static inline int accIdx() {
#ifdef _OPENMP
 return omp_get_thread_num();
#else
 return 0;
#endif
}

double EoSs::p(double e) {
#if defined TABLE
 return gsl_spline_eval(spline_p, e, acc_p[accIdx()]);
#elif defined SIMPLE
 return e / 3.;
#endif
}

double EoSs::dpe(double e) {
#if defined TABLE
 const int ia = accIdx();
 return (gsl_spline_eval(spline_p, e * 1.1, acc_p[ia]) -
         gsl_spline_eval(spline_p, e, acc_p[ia])) / (0.1 * e);
#elif defined SIMPLE
 return 1. / 3.;
#endif
}

double EoSs::t(double e) {
#if defined TABLE
 return gsl_spline_eval(spline_T, e, acc_T[accIdx()]);
#elif defined SIMPLE
 const double cnst =
     (16 + 0.5 * 21.0 * 2.5) * pow(C_PI, 2) / 30.0 / pow(0.197326968, 3);
 return e > 0. ? 1.0 * pow(e / cnst, 0.25) : 0.;
#endif
}

double EoSs::mu(double e) {
#if defined TABLE
 return gsl_spline_eval(spline_mu, e, acc_mu[accIdx()]);
#elif defined SIMPLE
 return 0.;
#endif
}
