#include <iostream>
#include <iomanip>
#include "eos.h"
#include "trancoeff.h"
#include "inc.h"

TransportCoeff::TransportCoeff(double _etaS, double _zetaS, EoS *_eos, int _etaSparam, double _ah, double _al, double _aRho, double _T0, double _etaSMin, double _eEtaSMin) {
 etaS0 = _etaS;
 zetaS0 = _zetaS;
 eos = _eos;
 etaSparam = _etaSparam;
 ah = _ah;
 al = _al;
 aRho =_aRho;
 etaSMin = _etaSMin;
 eEtaSMin = _eEtaSMin;
 T0=_T0;
}

void TransportCoeff::printZetaT()
{
 std::cout << "------zeta/s(T):\n";
 for(double e=0.1; e<3.0; e+=0.1){
  double T, mub, muq, mus, p;
  eos->eos(e, 0., 0., 0., T, mub, muq, mus, p);
  std::cout << std::setw(14) << T << std::setw(14) << zetaS(e, T) << std::endl;
 }
 std::cout << "---------------:\n";
}

double TransportCoeff::zetaS(double e, double T)
{
 return zetaS0 * (1. / 3. - eos->cs2(e)) / (exp((0.16 - T) / 0.001) + 1.);
}

void TransportCoeff::getEta(double e, double rho, double T, double &_etaS, double &_zetaS) {
 _etaS = etaS(e,rho, T);
 _zetaS = zetaS(e,T);
}

double TransportCoeff::etaS(double e,double rho, double T)
{
  if (etaSparam == 0){
      return etaS0;
  }
  else if (etaSparam == 1){
      return etaSMin +  ((T>T0) ? ah*(T-T0) :  al*(T-T0));
  }
  else if (etaSparam == 2){
      return std::max(0.0, etaSMin + ((e>eEtaSMin) ? ( (ah*(e-eEtaSMin)+aRho*rho) ): al*(e-eEtaSMin)+aRho*rho));
  }
}

void TransportCoeff::getTau(double e, double rho, double T, double &_taupi, double &_tauPi) {
 if (T > 0.) {
  _taupi = 5. / 5.068 * etaS(e,rho, T) / T;
  _tauPi = 1. / 5.068 * zetaS(e,T) / (15. * pow(0.33333-eos->cs2(e),2) * T);
 } else {
  _taupi = _tauPi = 0.;
 }
}
