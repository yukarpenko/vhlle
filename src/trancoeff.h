class EoS;

// this class contains the information about the transport coefficients
// of the fluid: eta/s, zeta/s and the corresponding relaxation times,
// taupi (\tau_\pi) and tauPi (\tau_\Pi)
class TransportCoeff {
public:
  const double etaS0, zetaS0, ah, al, aRho, T0, etaSMin, etaSEpsilonMin, etaSShiftMuB,
    etaSScaleMuB, zetaSPeakEpsilon, zetaSScaleBeta, zetaSSigmaMinus, zetaSSigmaPlus;
  EoS *eos;  // EoS instance is needed optionally for zeta/s or eta/s parametrization,
  const int etaSparam, zetaSparam;
  double zetaS(double e, double T, double s, double P);
  double etaS(double e,double rho, double T, double muB, double s, double P);
  TransportCoeff(double _etaS0, double _zetaS0, double _ah, double _al, double _aRho, double _T0, 
    double _etaSMin, double _etaSEpsilonMin, double _etaSShiftMuB, double _etaSScaleMuB, double _zetaSPeakEpsilon,
    double _zetaSScaleBeta, double _zetaSSigmaMinus, double _zetaSSigmaPlus, EoS *_eos, 
    int _etaSparam, int _zetaSparam) : etaS0(_etaS0), zetaS0(_zetaS0),  ah(_ah), al(_al), aRho(_aRho), T0(_T0),
    etaSMin(_etaSMin), etaSEpsilonMin(_etaSEpsilonMin), etaSShiftMuB(_etaSShiftMuB), etaSScaleMuB(_etaSScaleMuB),
    zetaSPeakEpsilon(_zetaSPeakEpsilon), zetaSScaleBeta(_zetaSScaleBeta), zetaSSigmaMinus(_zetaSSigmaMinus), 
    zetaSSigmaPlus(_zetaSSigmaPlus), eos(_eos), etaSparam(_etaSparam), zetaSparam(_zetaSparam){};
  ~TransportCoeff(){};
  void printZetaT();
  // returns (optionally temperature dependent) eta/s and zeta/s
  void getEta(double e, double rho, double T, double muB, double s, double P, double &_etaS, double &_zetaS);
  // returns shear and bulk relaxation times
  void getTau(double e, double rho, double T, double muB, double s, double P, double &_taupi, double &_tauPi);
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
