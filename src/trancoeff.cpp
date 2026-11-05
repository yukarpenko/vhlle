#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "eos.h"
#include "trancoeff.h"
#include "inc.h"


void TransportCoeff::printZetaT()
{
 std::cout << "------zeta(T):\n";
 for(double e=0.1; e<3.0; e+=0.1){
  double T, mub, muq, mus, p;
  eos->eos(e, 0., 0., 0., T, mub, muq, mus, p);
  const double s = eos->s(e, 0, 0, 0);
  std::cout << std::setw(14) << T << std::setw(14) << zeta(e, p, 0., s, T, mub)/s << std::endl;
 }
 std::cout << "---------------:\n";
}

void TransportCoeff::printParameters()
{
 std::cout << "====== transport coefficients ======\n";
 std::cout << "etaSparam = " << etaSparam << std::endl;
 if (etaSparam == 0) {
  std::cout << "eta/s = " << etaS0 << std::endl;
 } else if (etaSparam == 1) {
  std::cout << "etaSMin = " << etaSMin << "  T0 = " << T0
            << "  ah = " << ah << "  al = " << al << std::endl;
 } else if (etaSparam == 3) {
  std::cout << "etaSMin = " << etaSMin << "  T0 = " << T0
            << "  ah = " << ah << "  al = " << al
            << "  etaSShiftMuB = " << etaSShiftMuB
            << "  etaSScaleMuB = " << etaSScaleMuB << std::endl;
 } else if (etaSparam == 4) {
  std::cout << "eta*T/(e+p) at mu_B=0, 0.2, 0.4 GeV: " << etaTld0 << "  "
            << etaTld2 << "  " << etaTld4 << "  [arXiv:2408.00537]" << std::endl;
 } else if (etaSparam == 90) {
  std::cout << "[legacy] etaSMin = " << etaSMin << "  etaSEpsilonMin = " << etaSEpsilonMin
            << "  ah = " << ah << "  al = " << al << "  aRho = " << aRho << std::endl;
 } else {
  std::cout << "WARNING: unknown etaSparam, eta = 0 will be used\n";
 }
 std::cout << "zetaSparam = " << zetaSparam << std::endl;
 if (zetaSparam == 4) {
  std::cout << "zetaH_max = " << zetaH_max << "  T_zeta0 = " << zetaTz0
            << "  curvature = " << zetaTcurv << "  sigma_- = " << zetaSigM
            << "  sigma_+ = " << zetaSigP << "  [arXiv:2408.00537]" << std::endl;
 } else if (zetaSparam == 90) {
  std::cout << "[legacy] zetaS0 = " << zetaS0 << "  zetaSPeakEpsilon = " << zetaSPeakEpsilon
            << "  zetaSScaleBeta = " << zetaSScaleBeta
            << "  zetaSSigmaMinus = " << zetaSSigmaMinus
            << "  zetaSSigmaPlus = " << zetaSSigmaPlus << std::endl;
 } else if (zetaSparam >= 0 && zetaSparam <= 3) {
  std::cout << "zetaS0 = " << zetaS0 << std::endl;
 } else {
  std::cout << "WARNING: unknown zetaSparam, zeta = 0 will be used\n";
 }
 std::cout << "====================================\n";
}

double TransportCoeff::zeta(double e, double p, double nb, double s, double T, double mub)
{
 double T_p=0.180;

 double T_peak=0.165;
 double T_width=0.010;
 double B_norm=0.24;
 double B_width=1.5;
 double T_peak2=0.160;
 double B_norm2=0.13;
 double B1=0.01;
 double B2=0.12;

 if(zetaSparam==0)
    return (zetaS0 * (1. / 3. - eos->cs2(e)) / (exp((0.16 - T) / 0.001) + 1.)) * s;
 else if(zetaSparam==1)
 {
    if(T<0.180)
       return (0.03+(0.08*exp(((T/T_p)-1.)/(0.0025)))+(0.22*exp(((T/T_p)-1)/(0.022)))) * s;
    else if(T>=0.180 && T<0.200)
       return (27.55*(T/T_p)-13.45-(13.77*(T/T_p)*(T/T_p))) * s;
    else if(T>=0.200)
       return (0.001+(0.9*exp(-((T/T_p)-1.)/(0.0025)))+(0.25*exp(-((T/T_p)-1.)/(0.13)))) * s;
 }
 else if(zetaSparam==2)
 {
    if(T>T_peak)
       return (B_norm*((B_width*B_width)/((((T/T_peak)-1.)*((T/T_peak)-1.))+(B_width*B_width)))) * s;
    else if(T<=T_peak)
       return (B_norm*(exp(-((T-T_peak)/T_width)*((T-T_peak)/T_width)))) * s;
 }
 else if(zetaSparam==3)
 {
    if(T<T_peak2)
       return (B_norm2*exp(-((T-T_peak2)*(T-T_peak2)/(B1*B1)))) * s;
    else if(T>=T_peak2)
       return (B_norm2*exp(-((T-T_peak2)*(T-T_peak2)/(B2*B2)))) * s;
 }
 else if(zetaSparam==4)
 {
    // arXiv:2408.00537 (Jahan, Roch, Shen), Eq. (10):
    //   zeta_tld = zeta*T/(e+p) = zeta_max * exp(-[T - T_zeta(mu_B)]^2 / (2 sigma_{zeta,-/+}^2))
    //   for T < T_zeta / T >= T_zeta,  T_zeta(mu_B) = T_{zeta,0} - 0.15/GeV * mu_B^2.
    // zeta_max = zetaS0; the other parameters are class members (Table V).
    if (T > 0.) {
       const double T_zeta = zetaTz0 - zetaTcurv * mub * mub;
       const double dT = T - T_zeta;
       const double sigma = (T < T_zeta) ? zetaSigM : zetaSigP;
       const double zeta_tld = zetaH_max * exp(-dT * dT / (2. * sigma * sigma));
       return std::max(0.0, zeta_tld) * (e + p) / T;
    } else
       return 0.;
 }
 else if(zetaSparam==90)
 {  // legacy parametrization (formerly zetaSparam==4), not intended for normal use
    if (abs(s) < 1e-10 || T < 1e-1) {
            return 0;
         }
    if ( e < zetaSPeakEpsilon)
      return ((e+p)/T)*zetaS0 * exp(-pow(zetaSScaleBeta*(pow(e,0.25)-pow(zetaSPeakEpsilon,0.25)),2) / (2.0*pow(zetaSSigmaMinus,2)));
    else
      return ((e+p)/T)*zetaS0 * exp(-pow(zetaSScaleBeta*(pow(e,0.25)-pow(zetaSPeakEpsilon,0.25)),2) / (2.0*pow(zetaSSigmaPlus,2)));
 }
 return 0.;
}

void TransportCoeff::getEta(double e, double p, double nb, double s, double T, double mub, double &_eta, double &_zeta) {
 _eta = eta(e, p, nb, s, T, mub);
 _zeta = zeta(e, p, nb, s, T, mub);
}

double TransportCoeff::eta(double e, double p, double nb, double s, double T, double mub)
{
   if (etaSparam == 0){
      return etaS0*s;
   }
   else if (etaSparam == 1){
      return (etaSMin +  ((T>T0) ? ah*(T-T0) :  al*(T-T0))) * s;
   }
   else if (etaSparam == 3){
         if (abs(s) < 1e-10 || T < 1e-10) {
            return 0;
         }
         double Tc=T0+etaSShiftMuB*(mub/0.938);
         return ((e+p)/T)*std::max(0.0, (etaSMin + ((T>Tc) ? ah*(T-Tc) :  al*(T-Tc)))*(1 + etaSScaleMuB*(mub/0.938)));
   }
   else if (etaSparam == 4){
      // arXiv:2408.00537 (Jahan, Roch, Shen), Eq. (9): eta_tld = eta*T/(e+p),
      // piecewise linear in mu_B with nodes at mu_B = 0, 0.2, 0.4 GeV and
      // constant above 0.4 GeV.  |mu_B| is used, i.e. eta_tld(-mu_B) = eta_tld(mu_B).
      if (T > 0.) {
         const double muB = std::fabs(mub);
         double eta_tld;
         if (muB < 0.2)
            eta_tld = etaTld0 + (etaTld2 - etaTld0) * muB / 0.2;
         else if (muB < 0.4)
            eta_tld = etaTld2 + (etaTld4 - etaTld2) * (muB - 0.2) / 0.2;
         else
            eta_tld = etaTld4;
         return std::max(0.0, eta_tld) * (e + p) / T;
      } else
         return 0.;
   }
   else if (etaSparam == 90){
      // legacy parametrization (formerly etaSparam==2), not intended for normal use
         return (std::max(0.0, etaSMin + ((e>etaSEpsilonMin) ? ( (ah*(e-etaSEpsilonMin)+aRho*nb) ): al*(e-etaSEpsilonMin)+aRho*nb))) * s;
   }
  return 0.;
}

void TransportCoeff::getTau(double e, double p, double nb, double s, double T, double mub, double &_taupi, double &_tauPi) {
 // tau = C*eta/(e+p) for all parametrizations; for etaSparam==4 and
 // zetaSparam==4 this equals C*eta_tld/T, Eq. (11) of arXiv:2408.00537.
 if (T > 1e-10 && abs(e + p) > 1.e-10) {
   //based on 1403.0962.  As no cs2 table is available, we approximate it by 0.228 (??)
  _taupi = (5./5.068)* eta(e, p, nb, s, T, mub) / (e + p);
  _tauPi =  (6./5.068)*zeta(e, p, nb, s, T, mub) / (e + p);
 } else {
  _taupi = _tauPi = 0.;
 }
}
