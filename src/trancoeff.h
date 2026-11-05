class EoS;

// this class contains the information about the transport coefficients
// of the fluid: eta, zeta and the corresponding relaxation times,
// taupi (\tau_\pi) and tauPi (\tau_\Pi)
//
// etaSparam:  0 - constant eta/s = etaS0
//             1 - eta/s(T), linear in T with a kink at T0
//             3 - eta*T/(e+p) linear in T with a mu_B-dependent kink and scale
//             4 - eta*T/(e+p) = eta_tld(mu_B), arXiv:2408.00537, Eq. (9)
//            90 - [legacy, not intended for normal use] eta/s(e, n_B),
//                 linear in e with a kink at etaSEpsilonMin, linear in n_B
// zetaSparam: 0 - zeta/s = zetaS0 * (1/3 - c_s^2), cut off below T=0.16 GeV
//             1 - arXiv:1910.12930
//             2 - arXiv:1910.12930
//             3 - arXiv:2103.09848
//             4 - zeta*T/(e+p) = zeta_tld(T, mu_B), arXiv:2408.00537, Eq. (10),
//                 with zeta_max = zetaS0
//            90 - [legacy, not intended for normal use] zeta*T/(e+p) as an
//                 asymmetric Gaussian in e^{1/4}
class TransportCoeff {
public:
  const double etaS0, zetaS0;
  const int etaSparam, zetaSparam;
  EoS *eos;  // EoS instance is needed optionally for zeta/s or eta/s parametrization,
             // which depends on the speed of sound
private:
  // ---- parameters of the individual parametrizations -----------------------
  // Hardcoded defaults; they are not constructor arguments.
  // etaSparam==1, 3, 90: linear dependence on T (or e)
  double etaSMin = 0.08;   // value of eta/s (or eta_tld) at the kink
  double ah = 0.;          // slope above the kink
  double al = 0.;          // slope below the kink
  double T0 = 0.15;        // [GeV] kink temperature, etaSparam==1,3
  // etaSparam==3: mu_B dependence of the kink position and of the overall scale
  double etaSShiftMuB = 0.;
  double etaSScaleMuB = 0.;
  // etaSparam==90: e and n_B dependence
  double etaSEpsilonMin = 5.;  // [GeV/fm^3] position of the kink in e
  double aRho = 0.;            // [fm^3] slope in n_B
  // etaSparam==4: arXiv:2408.00537 Eq. (9), eta_tld = eta*T/(e+p) at
  // mu_B = 0, 0.2, 0.4 GeV; values: highest-likelihood set, Table V
  double etaTld0 = 0.045;  // eta_0
  double etaTld2 = 0.280;  // eta_2
  double etaTld4 = 0.287;  // eta_4
  // zetaSparam==4: arXiv:2408.00537 Eq. (10), zeta_tld = zeta*T/(e+p);
  // values: highest-likelihood set, Table V.  zeta_max is taken from zetaS0.
  double zetaH_max = 0.148;  // zeta/H, max
  double zetaTz0 = 0.214;    // [GeV] T_{zeta,0}
  double zetaTcurv = 0.15;   // [1/GeV] T_zeta(mu_B) = T_{zeta,0} - zetaTcurv*mu_B^2
  double zetaSigM = 0.040;   // [GeV] sigma_{zeta,-}
  double zetaSigP = 0.018;   // [GeV] sigma_{zeta,+}
  // zetaSparam==90: asymmetric Gaussian in e^{1/4}, peak height zetaS0
  double zetaSPeakEpsilon = 5.;  // [GeV/fm^3] position of the peak
  double zetaSScaleBeta = 0.103;
  double zetaSSigmaMinus = 0.1;
  double zetaSSigmaPlus = 0.1;
public:
  // zeta and eta [1/fm^3] (the same units as entropy density)
  double zeta(double e, double p, double nb, double s, double T, double mub);
  double eta(double e, double p, double nb, double s, double T, double mub);
  TransportCoeff(double _etaS0, double _zetaS0, int _etaSparam, int _zetaSparam, EoS *_eos)
   : etaS0(_etaS0), zetaS0(_zetaS0), etaSparam(_etaSparam), zetaSparam(_zetaSparam),
     eos(_eos) {};
  ~TransportCoeff(){};
  void printZetaT();
  // prints the parametrization choice and the parameters which are in use
  void printParameters();
  // returns (optionally temperature dependent) eta and zeta
 void getEta(double e, double p, double nb, double s, double T, double mub, double &_eta, double &_zeta);
  // returns shear and bulk relaxation times
 void getTau(double e, double p, double nb, double s, double T, double mub, double &_taupi, double &_tauPi);
  // deltapipi, taupipi, lambdapiPi * divided by tau_pi * !
  void getOther(double e, double nb, double nq, double ns,
    double &deltapipi, double &taupipi, double &lambdapiPi, double &phi7) {
    deltapipi = 4./3.;  taupipi = 10./7.;  lambdapiPi = 6./5.;
    phi7 = 9./70./eos->p(e, nb, nq, ns);
    if(std::isinf(phi7)) phi7=0.0;
  }
  void getOtherBulk(double e, double nb, double nq, double ns,
    double &delPiPi, double &lamPipi) {
    delPiPi = 2./3.;  lamPipi = 8./5.*(1. / 3. - eos->cs2(e));
  }
  // isViscous tells whether the fluid is viscous or inviscid
  inline bool isViscous() {
    if (etaS0 > 0. || zetaS0 > 0. || etaSparam>0 || zetaSparam>0)
    return true;
    else
    return false;
  }
};
