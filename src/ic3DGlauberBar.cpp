// ===========================================================================
//  Ic3DGlauberBar : 3D initial state of arXiv:2211.16408 (Du, Shen, Jeon, Gale)
//                   built on GLISSANDO participant lists.
//
//  ------------------------------------------------------------------------
//  !!! Three typographic errors in the published paper are corrected here.
//
//  (1) yCM(x_perp).  The paper writes  yCM = arctan[r tanh(ybeam)].  The same
//      quantity is *derived* from local energy-momentum conservation in the
//      ancestor paper (Shen & Alzhrani, PRC 102, 014909, Eq. 5) as
//      yCM = arctanh[r tanh(ybeam)], which is the only form that returns a
//      genuine rapidity.  arctan and arctanh agree only for small arguments
//      and differ strongly as r tanh(ybeam) -> +-1 (peripheral regions).
//      IMPLEMENTED:  yCM = atanh[ r tanh(ybeam) ].
//
//  (2) Eq. (A.8).  As typeset, both theta-branches of f^B_pm carry the SAME
//      width sigma_{B,-+}, so f^B_pm would collapse to a plain Gaussian, and
//      then f^B_+ (width sigma_{B,-}) and f^B_- (width sigma_{B,+}) would NOT
//      be mirror images of each other.  For a symmetric system such as Au+Au
//      that produces a large spurious forward/backward net-baryon asymmetry
//      (at 19.6 GeV the projectile peak comes out 3x wider than the target
//      peak), and it contradicts Fig. 1(b) where the two curves are plainly
//      mirror images.  The structure is inherited from Denicol et al.
//      (arXiv:1804.10557) and Shen & Alzhrani Eqs. (17)-(18), where the two
//      branches carry DIFFERENT widths: narrow on the side facing the beam
//      ("out"), wide on the side facing midrapidity ("in").  Matching Table I
//      against Shen & Alzhrani's table identifies
//          sigma_{B,+} = sigma_out ,   sigma_{B,-} = sigma_in
//      (sigma_{B,-} is systematically smaller than Shen & Alzhrani's sigma_in,
//      exactly as the paper states, because the new plateau term now supplies
//      the midrapidity baryons).
//      IMPLEMENTED:  asymmetric Gaussian, f^B_-(eta) = f^B_+(-eta), with
//          f^B_+ :  width sigma_{B,+} for eta > +eta_0^B/2   (outward)
//                   width sigma_{B,-} for eta < +eta_0^B/2   (inward)
//          normalization N = 1/[ sqrt(pi/2) (sigma_{B,+} + sigma_{B,-}) ].
//
//  (3) Eq. (A.9).  Taken literally the plateau is flat for |eta| <= eta_0^B,
//      i.e. of full width 2*eta_0^B.  With the Table I values that plateau
//      extends BEYOND the beam rapidity at every one of the four energies
//      (2.30 vs ybeam 2.09 at 7.7 GeV; 3.20 vs 3.04 at 19.6; 5.40 vs 4.20 at
//      62.4; 7.00 vs 5.36 at 200), which is unphysical.  The main text says
//      the plateau width equals the distance between the peaks of f^B_+ and
//      f^B_-, which is eta_0^B since eta_0^{B,pm} = +-eta_0^B/2.  That is
//      self-consistent, stays inside ybeam at all energies, and reproduces
//      Fig. 1(b) where the plateau edges sit on the two peaks.
//      IMPLEMENTED:  eta_0^B is the FULL plateau width, i.e. f^B_c is flat for
//          |eta| <= eta_0^B/2, with half-Gaussian tails of width sigma_{eta,B}
//          and normalization N' = 1/[ eta_0^B + sqrt(2 pi) sigma_{eta,B} ].
//
//  As an independent check of the parameter reading, eta_0^B/2 reproduces the
//  baryon peak positions eta_{B,0} of Shen & Alzhrani almost exactly:
//      7.7 GeV 1.15 vs 1.05 | 19.6 GeV 1.60 vs 1.50
//     62.4 GeV 2.70 vs 2.70 |  200 GeV 3.50 vs 3.50
// ===========================================================================

#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "eos.h"
#include "fld.h"
#include "ic3DGlauberBar.h"
#include "inc.h"
#include "rmn.h"
#include "s95p.h"

using namespace std;

namespace {
const double mN = 0.939;  // nucleon mass [GeV]
}

// ---------------------------------------------------------------------------
//  construction
// ---------------------------------------------------------------------------

Ic3DGlauberBar::Ic3DGlauberBar(Fluid *f, const char *filename, double _tau0,
                               const char *setup) {
 cout << "=== Ic3DGlauberBar: 3D IC of arXiv:2211.16408 (Du-Shen-Jeon-Gale) "
         "===\n";
 nx = f->getNX();
 ny = f->getNY();
 nz = f->getNZ();
 dx = f->getDx();
 dy = f->getDy();
 dz = f->getDz();
 xmin = f->getX(0);
 xmax = f->getX(nx - 1);
 ymin = f->getY(0);
 ymax = f->getY(ny - 1);
 zmin = f->getZ(0);
 zmax = f->getZ(nz - 1);
 tau0 = _tau0;

 // ---- defaults ----
 sNN = 19.6;
 s0 = 5.85;
 eta0s = 1.5;
 sigEtaS = 0.25;
 NB = 0.6;
 sigBplus = 0.19;
 sigBminus = 0.62;
 yLB = 0.3;
 Nc = 0.65;
 etaB0 = 3.2;
 sigEtaB = 0.1;

 w = 0.5;
 ZoverA = 0.4;
 nCut = 4.0;

 recenter = 1;
 rotatePP = 0;
 autoNormB = 1;
 autoNormS = 0;
 eosInvert = 1;
 dumpFile = "";

 readParameters(setup);

 ybeam = acosh(sNN / (2. * mN));

 // sanity: the baryon plateau must not stick out beyond the beam rapidity
 if (0.5 * etaB0 > ybeam) {
  cout << "Ic3DGlauberBar: WARNING baryon plateau half width " << 0.5 * etaB0
       << " exceeds ybeam " << ybeam << " -- check etaB0 and sNN.\n";
 }

 cout << "---- collision system ----\n";
 cout << "  sqrt(s_NN) = " << sNN << " GeV,  ybeam = " << ybeam << endl;
 cout << "  tau0       = " << tau0 << " fm/c" << endl;
 cout << "---- Table I parameters (arXiv:2211.16408) ----\n";
 cout << "  s0    = " << s0 << "   eta0s      = " << eta0s
      << "   sigEtaS = " << sigEtaS << endl;
 cout << "  NB    = " << NB << "   sigma_{B,+} = " << sigBplus
      << "   sigma_{B,-} = " << sigBminus << endl;
 cout << "  y_L^B = " << yLB << "   Nc         = " << Nc
      << "   eta0B   = " << etaB0 << "   sigEtaB = " << sigEtaB << endl;
 cout << "---- technical ----\n";
 cout << "  nucleon width w = " << w << " fm,   nq/nb = " << ZoverA << endl;
 cout << "  recenter = " << recenter << ", rotatePP = " << rotatePP << endl;
 cout << "  autoNormB = " << autoNormB << ", autoNormS = " << autoNormS
      << ", eosInvert = " << eosInvert << endl;

 allocate();
 loadEvents(filename);
}

Ic3DGlauberBar::~Ic3DGlauberBar() {
 for (int ix = 0; ix < nx; ix++) {
  delete[] TA[ix];
  delete[] TB[ix];
 }
 delete[] TA;
 delete[] TB;
}

void Ic3DGlauberBar::allocate() {
 TA = new double *[nx];
 TB = new double *[nx];
 for (int ix = 0; ix < nx; ix++) {
  TA[ix] = new double[ny];
  TB[ix] = new double[ny];
  for (int iy = 0; iy < ny; iy++) {
   TA[ix][iy] = 0.0;
   TB[ix][iy] = 0.0;
  }
 }
}

// ---------------------------------------------------------------------------
//  parameter input.  'setup' is either the name of a built-in preset holding
//  the Table I row for one of the four energies of the paper, or the path to
//  a parameter file with "name value" lines.
// ---------------------------------------------------------------------------

void Ic3DGlauberBar::applyPreset(const string &name) {
 // Table I of arXiv:2211.16408, plus tau0 from Shen & Alzhrani (Ref. [32]).
 // NOTE: s0 and NB are quoted for the thickness function normalization used
 // in the paper.  The GLISSANDO based T_pm built here is normalized
 // differently, so these two numbers are only starting values -- use
 // autoNormB / autoNormS, or retune them against dN/deta and dN^{p-pbar}/dy.
 if (name == "AuAu7.7") {
  sNN = 7.7;  s0 = 2.45; eta0s = 0.75; sigEtaS = 0.17;
  NB = 0.8;   sigBplus = 0.15; sigBminus = 0.58; yLB = 0.2;
  Nc = 0.6;   etaB0 = 2.3;  sigEtaB = 0.07;
 } else if (name == "AuAu19.6") {
  sNN = 19.6; s0 = 5.85; eta0s = 1.5;  sigEtaS = 0.25;
  NB = 0.6;   sigBplus = 0.19; sigBminus = 0.62; yLB = 0.3;
  Nc = 0.65;  etaB0 = 3.2;  sigEtaB = 0.10;
 } else if (name == "AuAu62.4") {
  sNN = 62.4; s0 = 11.0; eta0s = 2.3;  sigEtaS = 0.28;
  NB = 0.65;  sigBplus = 0.25; sigBminus = 0.95; yLB = 0.25;
  Nc = 0.51;  etaB0 = 5.4;  sigEtaB = 0.22;
 } else if (name == "AuAu200") {
  sNN = 200.; s0 = 16.0; eta0s = 2.5;  sigEtaS = 0.58;
  NB = 0.55;  sigBplus = 0.3;  sigBminus = 1.15; yLB = 0.1;
  Nc = 0.41;  etaB0 = 7.0;  sigEtaB = 0.25;
 } else {
  cout << "Ic3DGlauberBar: unknown preset '" << name << "'\n";
  cout << "  known presets: AuAu7.7 AuAu19.6 AuAu62.4 AuAu200\n";
  cout << "  or give the path to a parameter file instead.\n";
  exit(1);
 }
 cout << "Ic3DGlauberBar: using built-in preset " << name << endl;
}

void Ic3DGlauberBar::readParameters(const string &setup) {
 if (setup.empty()) {
  cout << "Ic3DGlauberBar: no setup given, using built-in defaults\n";
  return;
 }
 ifstream fin(setup.c_str());
 if (!fin.is_open()) {
  // not a file -> treat as a preset name
  applyPreset(setup);
  return;
 }
 cout << "Ic3DGlauberBar: reading parameters from " << setup << endl;
 string line;
 while (getline(fin, line)) {
  // strip comments
  size_t hash = line.find_first_of("#!");
  if (hash != string::npos) line = line.substr(0, hash);
  istringstream sline(line);
  string parName, parValue;
  if (!(sline >> parName >> parValue)) continue;

  if (parName == "preset")            applyPreset(parValue);
  else if (parName == "sNN")          sNN = atof(parValue.c_str());
  else if (parName == "s0")           s0 = atof(parValue.c_str());
  else if (parName == "eta0s")        eta0s = atof(parValue.c_str());
  else if (parName == "sigEtaS")      sigEtaS = atof(parValue.c_str());
  else if (parName == "NB")           NB = atof(parValue.c_str());
  else if (parName == "sigBplus")     sigBplus = atof(parValue.c_str());
  else if (parName == "sigBminus")    sigBminus = atof(parValue.c_str());
  else if (parName == "yLB")          yLB = atof(parValue.c_str());
  else if (parName == "Nc")           Nc = atof(parValue.c_str());
  else if (parName == "etaB0")        etaB0 = atof(parValue.c_str());
  else if (parName == "sigEtaB")      sigEtaB = atof(parValue.c_str());
  else if (parName == "w")            w = atof(parValue.c_str());
  else if (parName == "ZoverA")       ZoverA = atof(parValue.c_str());
  else if (parName == "nCut")         nCut = atof(parValue.c_str());
  else if (parName == "recenter")     recenter = atoi(parValue.c_str());
  else if (parName == "rotatePP")     rotatePP = atoi(parValue.c_str());
  else if (parName == "autoNormB")    autoNormB = atoi(parValue.c_str());
  else if (parName == "autoNormS")    autoNormS = atoi(parValue.c_str());
  else if (parName == "eosInvert")    eosInvert = atoi(parValue.c_str());
  else if (parName == "dumpProfiles") dumpFile = parValue;
  else cout << "Ic3DGlauberBar: unknown parameter '" << parName << "'\n";
 }
}

// ---------------------------------------------------------------------------
//  GLISSANDO event reading.
//  Same format as IcGlissando:  x  y  C  w   with events separated by a line
//  that fails to parse.  C>0 -> projectile (T_+), C<0 -> target (T_-),
//  C==0 -> binary collision source, skipped.
// ---------------------------------------------------------------------------

void Ic3DGlauberBar::loadEvents(const string &filename) {
 ifstream fin(filename.c_str());
 if (!fin.good()) {
  cout << "Ic3DGlauberBar: I/O error with " << filename << endl;
  exit(1);
 }

 vector<double> xA, yA, xB, yB;
 nevents = 0;
 long npartTot = 0;
 string line;
 istringstream instream;

 while (!fin.eof()) {
  getline(fin, line);
  instream.str(line);
  instream.seekg(0);
  instream.clear();
  double xr, yr, wr;
  int cr;
  instream >> xr >> yr >> cr >> wr;
  if (!instream.fail()) {
   if (cr > 0) {
    xA.push_back(xr);
    yA.push_back(yr);
   } else if (cr < 0) {
    xB.push_back(xr);
    yB.push_back(yr);
   }
   // cr == 0 : binary collision entry, ignored
  } else if (xA.size() + xB.size() > 0) {
   npartTot += (long)(xA.size() + xB.size());
   addEvent(xA, yA, xB, yB);
   xA.clear(); yA.clear(); xB.clear(); yB.clear();
   nevents++;
  }
 }
 // flush a trailing event not followed by a blank line
 if (xA.size() + xB.size() > 0) {
  npartTot += (long)(xA.size() + xB.size());
  addEvent(xA, yA, xB, yB);
  nevents++;
 }
 fin.close();

 if (nevents == 0) {
  cout << "Ic3DGlauberBar: no events read from " << filename << endl;
  exit(1);
 }

 // ensemble average of the thickness functions -- done BEFORE the non-linear
 // longitudinal profiles are applied, exactly as in the paper.
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++) {
   TA[ix][iy] /= (double)nevents;
   TB[ix][iy] /= (double)nevents;
  }

 npartAvg = (double)npartTot / (double)nevents;
 cout << "Ic3DGlauberBar: read " << nevents << " GLISSANDO events, <Npart> = "
      << npartAvg << endl;

 // check of the transverse normalization: int T d^2x should give <Npart>/2
 double sumA = 0., sumB = 0.;
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++) {
   sumA += TA[ix][iy] * dx * dy;
   sumB += TB[ix][iy] * dx * dy;
  }
 cout << "  int T_+ d2x = " << sumA << ",  int T_- d2x = " << sumB
      << "   (sum = " << sumA + sumB << ", should match <Npart>)" << endl;
 if (sumA + sumB < 0.9 * npartAvg)
  cout << "  WARNING: transverse grid too small or nCut too small -- "
          "part of the thickness function falls outside the grid.\n";
}

void Ic3DGlauberBar::addEvent(vector<double> &xA, vector<double> &yA,
                              vector<double> &xB, vector<double> &yB) {
 const size_t nA = xA.size(), nB = xB.size();
 const size_t nP = nA + nB;
 if (nP == 0) return;

 // ---- recentre on the participant centre of mass ----
 if (recenter) {
  double xc = 0., yc = 0.;
  for (size_t i = 0; i < nA; i++) { xc += xA[i]; yc += yA[i]; }
  for (size_t i = 0; i < nB; i++) { xc += xB[i]; yc += yB[i]; }
  xc /= (double)nP;
  yc /= (double)nP;
  for (size_t i = 0; i < nA; i++) { xA[i] -= xc; yA[i] -= yc; }
  for (size_t i = 0; i < nB; i++) { xB[i] -= xc; yB[i] -= yc; }
 }

 // ---- rotate to the 2nd order participant plane ----
 if (rotatePP) {
  double s2 = 0., c2 = 0., r2 = 0.;
  for (size_t i = 0; i < nA; i++) {
   const double rr = xA[i] * xA[i] + yA[i] * yA[i];
   const double phi = atan2(yA[i], xA[i]);
   s2 += rr * sin(2. * phi); c2 += rr * cos(2. * phi); r2 += rr;
  }
  for (size_t i = 0; i < nB; i++) {
   const double rr = xB[i] * xB[i] + yB[i] * yB[i];
   const double phi = atan2(yB[i], xB[i]);
   s2 += rr * sin(2. * phi); c2 += rr * cos(2. * phi); r2 += rr;
  }
  if (r2 > 0.) {
   double psi2 = 0.5 * (atan2(s2, c2) + C_PI);
   // Psi2 is defined mod pi.  The remaining two-fold ambiguity would flip
   // x -> -x and average the tilt (hence v1) away, so it is fixed here by
   // demanding that the projectile sits at positive x relative to the target.
   double mxA = 0., mxB = 0.;
   for (size_t i = 0; i < nA; i++)
    mxA += xA[i] * cos(psi2) + yA[i] * sin(psi2);
   for (size_t i = 0; i < nB; i++)
    mxB += xB[i] * cos(psi2) + yB[i] * sin(psi2);
   if (nA > 0) mxA /= (double)nA;
   if (nB > 0) mxB /= (double)nB;
   if (mxA - mxB < 0.) psi2 += C_PI;

   const double cs = cos(psi2), sn = sin(psi2);
   for (size_t i = 0; i < nA; i++) {
    const double xr = xA[i] * cs + yA[i] * sn;
    const double yr = -xA[i] * sn + yA[i] * cs;
    xA[i] = xr; yA[i] = yr;
   }
   for (size_t i = 0; i < nB; i++) {
    const double xr = xB[i] * cs + yB[i] * sn;
    const double yr = -xB[i] * sn + yB[i] * cs;
    xB[i] = xr; yB[i] = yr;
   }
  }
 }

 // ---- Gaussian smearing onto the transverse grid ----
 const int nsx = (int)(nCut * w / dx) + 1;
 const int nsy = (int)(nCut * w / dy) + 1;
 const double norm = 1.0 / (2.0 * C_PI * w * w);
 const double inv2w2 = 1.0 / (2.0 * w * w);

 for (int side = 0; side < 2; side++) {
  vector<double> &xv = (side == 0) ? xA : xB;
  vector<double> &yv = (side == 0) ? yA : yB;
  double **T = (side == 0) ? TA : TB;
  for (size_t ip = 0; ip < xv.size(); ip++) {
   const int ixc = (int)((xv[ip] - xmin) / dx);
   const int iyc = (int)((yv[ip] - ymin) / dy);
   for (int ix = ixc - nsx; ix <= ixc + nsx; ix++) {
    if (ix < 0 || ix >= nx) continue;
    const double xd = xv[ip] - (xmin + ix * dx);
    for (int iy = iyc - nsy; iy <= iyc + nsy; iy++) {
     if (iy < 0 || iy >= ny) continue;
     const double yd = yv[ip] - (ymin + iy * dy);
     double g = norm * exp(-(xd * xd + yd * yd) * inv2w2);
     if (g != g || fabs(g) > DBL_MAX) g = 0.0;
     T[ix][iy] += g;
    }
   }
  }
 }
}

// ---------------------------------------------------------------------------
//  longitudinal profiles
// ---------------------------------------------------------------------------

// Eq. (A.7).  sign = +1 -> projectile f^s_+, sign = -1 -> target f^s_-
double Ic3DGlauberBar::fSpm(double eta, int sign) const {
 const double etaMax = ybeam;
 if (fabs(eta) >= etaMax) return 0.0;
 const double tilt = 1.0 + sign * eta / etaMax;  // linear tilt factor
 double plateau;
 const double ae = fabs(eta);
 if (ae <= eta0s)
  plateau = 1.0;
 else
  plateau = exp(-(ae - eta0s) * (ae - eta0s) / (2.0 * sigEtaS * sigEtaS));
 return tilt * plateau;
}

// Eq. (A.8), corrected [see (2) in the file header].  Asymmetric Gaussian:
// narrow width sigma_{B,+} on the outward side (facing the beam), wide width
// sigma_{B,-} on the inward side (facing midrapidity).  Mirror symmetric,
// f^B_-(eta) = f^B_+(-eta), and normalized to unit integral in eta.
// sign = +1 -> f^B_+ (peak at +eta0B/2), -1 -> f^B_- (peak at -eta0B/2).
double Ic3DGlauberBar::fBpm(double eta, int sign) const {
 const double peak = sign * 0.5 * etaB0;
 const double d = eta - peak;
 const bool outward = (sign > 0) ? (d > 0.0) : (d < 0.0);
 const double sig = outward ? sigBplus : sigBminus;
 const double nrm = 1.0 / (sqrt(0.5 * C_PI) * (sigBplus + sigBminus));
 return nrm * exp(-d * d / (2.0 * sig * sig));
}

// Eq. (A.9), corrected [see (3) in the file header].  eta0B is the FULL
// plateau width, so f^B_c is flat for |eta| <= eta0B/2 -- the plateau edges
// coincide with the two peaks of f^B_pm -- with half-Gaussian tails of width
// sigma_{eta,B}.  Normalized to unit integral in eta.
double Ic3DGlauberBar::fBc(double eta) const {
 const double half = 0.5 * etaB0;
 const double nrm = 1.0 / (etaB0 + sqrt(2.0 * C_PI) * sigEtaB);
 const double ae = fabs(eta);
 if (ae <= half) return nrm;
 return nrm * exp(-(ae - half) * (ae - half) / (2.0 * sigEtaB * sigEtaB));
}

double Ic3DGlauberBar::rAsym(double TAv, double TBv) const {
 const double sum = TAv + TBv;
 if (sum < 1e-12) return 0.0;
 return (TAv - TBv) / sum;
}

// Corrected form [see (1) in the file header]: arctanh, as derived from local
// energy-momentum conservation in Shen & Alzhrani, PRC 102, 014909, Eq. (5).
double Ic3DGlauberBar::yCM(double TAv, double TBv) const {
 const double r = rAsym(TAv, TBv);
 const double arg = r * tanh(ybeam);
 const double a = max(-0.999999, min(0.999999, arg));  // guard |arg| -> 1
 return atanh(a);
}

// ---------------------------------------------------------------------------
//  entropy density -> energy density.
//  s95p gives the nb = 0 estimate; when eosInvert is on this is refined by
//  bisecting the actual EoS relation s(e, nb, nq, ns) = sDens, which matters
//  at the low energies where mu_B is large.
// ---------------------------------------------------------------------------

double Ic3DGlauberBar::eFromS(EoS *eos, double sDens, double nb, double nq,
                              double ns) const {
 if (sDens <= 0.0) return 0.0;
 const double eGuess = s95p::s95p_e(sDens);
 if (!eosInvert || eGuess <= 0.0) return eGuess;

 // bracket around the zero-density guess (finite nb always raises e)
 double eLo = eGuess / 50.0, eHi = eGuess * 50.0;
 const double sLo = eos->s(eLo, nb, nq, ns);
 const double sHi = eos->s(eHi, nb, nq, ns);
 if (!(sLo < sDens && sHi > sDens)) return eGuess;  // fall back silently

 for (int it = 0; it < 200; it++) {
  const double eMid = 0.5 * (eLo + eHi);
  const double sMid = eos->s(eMid, nb, nq, ns);
  if (sMid < sDens)
   eLo = eMid;
  else
   eHi = eMid;
  if (eHi - eLo < 1e-12 * (eHi + eLo)) break;
 }
 return 0.5 * (eLo + eHi);
}

// ---------------------------------------------------------------------------
//  normalization helpers
// ---------------------------------------------------------------------------

double Ic3DGlauberBar::totalBaryon() const {
 // Nb = int n tau0 dx dy deta , with Bjorken flow (u^tau = 1)
 double Nb = 0.0;
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++) {
   const double TAv = TA[ix][iy], TBv = TB[ix][iy];
   if (TAv + TBv < 1e-12) continue;
   const double r = rAsym(TAv, TBv);
   const double ycm = yCM(TAv, TBv);
   const double supp = 1.0 / (cosh(r) * cosh(r));
   for (int iz = 0; iz < nz; iz++) {
    const double eta = zmin + iz * dz;
    const double etaL = eta - yLB * ycm;
    const double n = (NB / tau0) *
                     (fBpm(etaL, -1) * TBv + fBpm(etaL, +1) * TAv +
                      supp * Nc * fBc(eta) * (TAv + TBv));
    Nb += n * tau0 * dx * dy * dz;
   }
  }
 return Nb;
}

double Ic3DGlauberBar::totalEnergy(EoS *eos) const {
 double E = 0.0;
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++) {
   const double TAv = TA[ix][iy], TBv = TB[ix][iy];
   if (TAv + TBv < 1e-12) continue;
   for (int iz = 0; iz < nz; iz++) {
    const double eta = zmin + iz * dz;
    const double sDens = s0 * (fSpm(eta, -1) * TBv + fSpm(eta, +1) * TAv);
    if (sDens <= 0.0) continue;
    const double e = s95p::s95p_e(sDens);  // nb=0 estimate is enough here
    E += tau0 * e * dx * dy * dz * cosh(eta);
   }
  }
 return E;
}

// ---------------------------------------------------------------------------
//  main entry point
// ---------------------------------------------------------------------------

void Ic3DGlauberBar::setIC(Fluid *f, EoS *eos) {
 // ---- optional automatic normalizations ----
 if (autoNormB) {
  const double Nb0 = totalBaryon();
  if (Nb0 > 0.) {
   const double scale = npartAvg / Nb0;
   NB *= scale;
   cout << "Ic3DGlauberBar: autoNormB  NB rescaled by " << scale << " -> "
        << NB << endl;
   cout << "                (total baryon number matched to <Npart> = "
        << npartAvg << ")" << endl;
  }
 }
 if (autoNormS) {
  // target: all beam energy deposited into the fluid, as in IcGlissando
  const double Etarget = npartAvg * 0.5 * sNN;
  for (int it = 0; it < 100; it++) {
   const double E = totalEnergy(eos);
   if (E <= 0.) break;
   const double corr = pow(Etarget / E, 0.75);
   s0 *= corr;
   if (fabs(corr - 1.0) < 1e-5) break;
  }
  cout << "Ic3DGlauberBar: autoNormS  s0 -> " << s0
       << "  (target E = " << Etarget << " GeV)" << endl;
 }

 if (!dumpFile.empty()) dumpProfiles();

 // ---- fill the grid ----
 double E = 0.0, Px = 0.0, Py = 0.0, Pz = 0.0, Nb = 0.0, S = 0.0;
 double Jy0 = 0.0, Jint1 = 0.0, Jint3 = 0.0;
 double Xcm = 0.0, Ycm = 0.0, Zcm = 0.0, Tcm = 0.0;
 double E_midrap = 0.0, Jy0_midrap = 0.0;
 double eMax = 0.0, nbMax = 0.0;

 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++) {
   const double TAv = TA[ix][iy], TBv = TB[ix][iy];
   const double Tsum = TAv + TBv;
   // per-column quantities of the transverse-longitudinal coupling
   const double r = rAsym(TAv, TBv);
   const double ycm = yCM(TAv, TBv);
   const double supp = 1.0 / (cosh(r) * cosh(r));

   for (int iz = 0; iz < nz; iz++) {
    const double eta = zmin + iz * dz;

    // ---- entropy density, Eq. (A.6) ----
    const double sDens = s0 * (fSpm(eta, -1) * TBv + fSpm(eta, +1) * TAv);

    // ---- net baryon density, Eq. (1); note the shifted argument etaL ----
    const double etaL = eta - yLB * ycm;
    double nb = 0.0;
    if (Tsum > 1e-12)
     nb = (NB / tau0) * (fBpm(etaL, -1) * TBv + fBpm(etaL, +1) * TAv +
                         supp * Nc * fBc(eta) * Tsum);
    if (nb < 0.) nb = 0.;

    const double nq = ZoverA * nb;
    const double ns = 0.0;

    const double e = eFromS(eos, sDens, nb, nq, ns);
    if (e > eMax) eMax = e;
    if (nb > nbMax) nbMax = nb;

    // ---- Bjorken flow: u^mu = (1,0,0,0) in Milne components ----
    const double vx = 0.0, vy = 0.0, vz = 0.0;

    Cell *c = f->getCell(ix, iy, iz);
    c->setPrimVar(eos, tau0, e, nb, nq, ns, vx, vy, vz);
    if (e > 1e-10) c->setAllM(1.);

    // ---- bookkeeping ----
    if (e <= 0.) continue;
    const double p = eos->p(e, nb, nq, ns);
    const double coshEta = cosh(eta), sinhEta = sinh(eta);
    const double u0 = 1.0, u3 = 0.0;  // Milne components of Bjorken flow
    const double u0lab = u0 * coshEta + u3 * sinhEta;
    const double uzlab = u0 * sinhEta + u3 * coshEta;
    const double dE = tau0 * ((e + p) * u0 * u0lab - p * coshEta) * dx * dy * dz;
    E += dE;
    Pz += tau0 * ((e + p) * u0 * uzlab - p * sinhEta) * dx * dy * dz;
    Nb += nb * tau0 * dx * dy * dz;
    S += tau0 * eos->s(e, nb, nq, ns) * u0 * dx * dy * dz;

    const double t = tau0 * coshEta;
    const double z = tau0 * sinhEta;
    const double x = xmin + ix * dx;
    const double y = ymin + iy * dy;
    Xcm += x * dE; Ycm += y * dE; Zcm += z * dE; Tcm += t * dE;
    Jy0 += tau0 * (e + p) * u0 * (z * 0.0 - x * uzlab) * dx * dy * dz * gevtofm;
    Jint1 += 0.0;
    Jint3 += tau0 * (e + p) * u0 * uzlab * dx * dy * dz * gevtofm;
    if (iz > nz / 2 - 2 && iz < nz / 2 + 2) {
     E_midrap += dE;
     Jy0_midrap +=
         tau0 * (e + p) * u0 * (z * 0.0 - x * uzlab) * dx * dy * dz * gevtofm;
    }
   }
  }

 if (E > 0.) { Xcm /= E; Ycm /= E; Zcm /= E; Tcm /= E; }
 const double Jy = Jy0 - Zcm * Jint1 + Xcm * Jint3;

 cout << "---- Ic3DGlauberBar summary ----\n";
 cout << "  hydrodynamic E = " << E << " GeV,  Pz = " << Pz
      << ",  Nbar = " << Nb << endl;
 cout << "  Px = " << Px << "  Py = " << Py << endl;
 cout << "  initial_entropy S_ini = " << S << endl;
 cout << "  max energy density  = " << eMax << " GeV/fm^3" << endl;
 cout << "  max baryon density  = " << nbMax << " 1/fm^3" << endl;
 cout << "  Xcm: " << sqrt(fabs(Tcm * Tcm - Zcm * Zcm)) << "  " << Xcm << "  "
      << Ycm << endl;
 cout << "  initial/corrected J_y  " << Jy0 << " " << Jy << endl;
 cout << "  1/tau*dE/dy_ini: " << E_midrap / (3.0 * dz * tau0) << endl;
 cout << "  1/tau*dJ/dy_ini: " << Jy0_midrap / (3.0 * dz * tau0) << endl;
 cout << "--------------------------------\n";
}

// ---------------------------------------------------------------------------
//  diagnostic dump of the 1D longitudinal profiles (Fig. 1(a),(b) of the paper)
// ---------------------------------------------------------------------------

void Ic3DGlauberBar::dumpProfiles() const {
 ofstream fout(dumpFile.c_str());
 if (!fout.is_open()) {
  cout << "Ic3DGlauberBar: cannot open " << dumpFile << " for writing\n";
  return;
 }
 fout << "# eta   f^s_+   f^s_-   f^B_+   f^B_-   f^B_c\n";
 fout << "# arXiv:2211.16408 Fig.1(a),(b);  sqrt(s_NN) = " << sNN
      << " GeV, ybeam = " << ybeam << "\n";
 const int npts = 601;
 const double emin = -1.1 * ybeam, emax = 1.1 * ybeam;
 for (int i = 0; i < npts; i++) {
  const double eta = emin + (emax - emin) * i / (npts - 1);
  fout << setw(12) << eta << setw(14) << fSpm(eta, +1) << setw(14)
       << fSpm(eta, -1) << setw(14) << fBpm(eta, +1) << setw(14)
       << fBpm(eta, -1) << setw(14) << fBc(eta) << "\n";
 }
 fout.close();
 cout << "Ic3DGlauberBar: longitudinal profiles written to " << dumpFile
      << endl;
}
