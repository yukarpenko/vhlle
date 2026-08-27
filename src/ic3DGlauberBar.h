#pragma once
#include <string>
#include <vector>

class Fluid;
class EoS;

// ===========================================================================
//  Ic3DGlauberBar
//
//  3D initial state with tilted entropy source and tilted + plateau baryon
//  stopping, following
//
//    L. Du, C. Shen, S. Jeon, C. Gale,
//    "Probing initial baryon stopping and equation of state with
//     rapidity-dependent directed flow of identified particles",
//    Phys. Rev. C 108, L041901 (2023), arXiv:2211.16408
//
//  Transverse geometry (the thickness functions T_+ and T_-) is built from
//  participant nucleon lists produced by GLISSANDO, read in the same file
//  format as IcGlissando (x  y  C  w), where the sign of the integer flag C
//  selects projectile (C>0) or target (C<0); C==0 entries (binary collision
//  sources) are skipped.  Following the paper, the thickness functions are
//  averaged over the whole event ensemble FIRST, and the (non-linear)
//  longitudinal profiles are evaluated afterwards on the smooth averaged
//  T_+(x_perp), T_-(x_perp).
//
//  Entropy density                                       [paper Eq. (A.6-A.7)]
//    s(tau0,x_perp,eta) = s0 [ f^s_-(eta) T_- + f^s_+(eta) T_+ ]
//    f^s_pm(eta) = theta(etaMax-|eta|) (1 pm eta/etaMax)
//                  x [ theta(|eta|-eta0s) exp(-(|eta|-eta0s)^2/2 sigEtaS^2)
//                      + theta(eta0s-|eta|) ]
//    etaMax = ybeam.
//
//  Net baryon density                                        [paper Eq. (1)]
//    n(tau0,x_perp,eta) = (NB/tau0) { f^B_-[etaL] T_- + f^B_+[etaL] T_+
//                         + cosh^-2[r] Nc f^B_c(eta) (T_- + T_+) }
//    r(x_perp)    = (T_+ - T_-)/(T_+ + T_-)
//    yCM(x_perp)  = atanh[ r tanh(ybeam) ]
//    etaL(x_perp) = eta - yLB * yCM(x_perp)
//
//  The paper as typeset contains three typographic errors (in yCM, in Eq. A.8
//  and in Eq. A.9).  They are corrected here; see the header comment of
//  ic3DGlauberBar.cpp for the corrected forms and the reasoning.
// ===========================================================================

class Ic3DGlauberBar {
private:
 // ---- hydro grid (copied from Fluid) ----
 int nx, ny, nz;
 double xmin, xmax, ymin, ymax, zmin, zmax;
 double dx, dy, dz;
 double tau0;

 // ---- collision system ----
 double sNN;    // sqrt(s_NN) [GeV]
 double ybeam;  // beam rapidity, arccosh(sqrt(s_NN)/(2 m_N))

 // ---- Table I parameters of arXiv:2211.16408 ----
 double s0;         // s_0      : entropy normalization
 double eta0s;      // eta_0^s  : half width of entropy plateau
 double sigEtaS;    // sigma_{eta,s} : Gaussian fall-off of entropy plateau
 double NB;         // N_B      : baryon normalization
 double sigBplus;   // sigma_{B,+} : OUTWARD width of the baryon peaks
 double sigBminus;  // sigma_{B,-} : INWARD  width of the baryon peaks
 double yLB;        // y_L^B    : strength of the longitudinal shift for baryon density
 double Nc;         // N_c      : relative weight of the central plateau
 double etaB0;      // eta_0^B  : FULL width of plateau = peak separation
 double sigEtaB;    // sigma_{eta,B} : width of the plateau half-Gaussian tails

 // ---- auxiliary / technical parameters ----
 double w;       // Gaussian width of a single nucleon in transverse plane [fm]
 double ZoverA;  // n_q / n_b ratio used to initialize the electric charge
 double nCut;    // Gaussian smearing truncated at nCut*w

 // ---- event handling ----
 int recenter;  // shift every event to its participant centre of mass
 int rotatePP;  // rotate every event to its 2nd order participant plane

 // ---- normalization control ----
 int autoNormB;  // rescale NB so that the total baryon number = <Npart>
 int autoNormS;  // rescale s0 so that the total energy = <Npart>*sqrt(s)/2
 int eosInvert;  // 1: invert s(e,nb) with the actual EoS; 0: use s95p (nb=0)

 // ---- ensemble-averaged thickness functions, [nx][ny], units 1/fm^2 ----
 double **TA;  // T_+ , projectile (right moving)
 double **TB;  // T_- , target     (left moving)

 int nevents;
 double npartAvg;  // average number of participants per event

 std::string dumpFile;  // if non-empty, write the 1D profiles there

 // ---- internals ----
 void readParameters(const std::string &setup);
 void applyPreset(const std::string &name);
 void loadEvents(const std::string &filename);
 void addEvent(std::vector<double> &xA, std::vector<double> &yA,
               std::vector<double> &xB, std::vector<double> &yB);
 void allocate();

 // longitudinal profiles
 double fSpm(double eta, int sign) const;                 // Eq. (A.7)
 double fBpm(double eta, int sign) const;                 // Eq. (A.8)
 double fBc(double eta) const;                            // Eq. (A.9)
 double yCM(double TAv, double TBv) const;
 double rAsym(double TAv, double TBv) const;

 // entropy density -> energy density
 double eFromS(EoS *eos, double sDens, double nb, double nq, double ns) const;

 // normalization helpers
 double totalEnergy(EoS *eos) const;
 double totalBaryon() const;
 void dumpProfiles() const;

public:
 Ic3DGlauberBar(Fluid *f, const char *filename, double tau0,
                const char *setup);
 ~Ic3DGlauberBar();
 void setIC(Fluid *f, EoS *eos);
};
