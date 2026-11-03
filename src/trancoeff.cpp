#include <iostream>
#include <iomanip>
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

 // --- zetaSparam==4: arXiv:2408.00537 (Jahan, Roch, Shen), Eq. (10).
 // Highest-likelihood parameter set, Table V.  Note that this parametrization
 // is written for zeta_tld = zeta*T/(e+p), NOT for zeta/s, and that the
 // Gaussian widths enter as 2*sigma^2 (unlike zetaSparam 2 and 3).
 // zeta_max is NOT hardcoded: it is taken from zetaS0.
 const double zetaTz0   = 0.214;  // T_{zeta,0} [GeV]
 const double zetaTcurv = 0.15;   // curvature of T_zeta(mu_B) [1/GeV]
 const double zetaSigM  = 0.040;  // sigma_{zeta,-} [GeV]
 const double zetaSigP  = 0.018;  // sigma_{zeta,+} [GeV]

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
    if (abs(s) < 1e-10 || T < 1e-1) {
            return 0;
         }
    if ( e < zetaSPeakEpsilon)
      return ((e+p)/T)*zetaS0 * exp(-pow(zetaSScaleBeta*(pow(e,0.25)-pow(zetaSPeakEpsilon,0.25)),2) / (2.0*pow(zetaSSigmaMinus,2)));
    else
      return ((e+p)/T)*zetaS0 * exp(-pow(zetaSScaleBeta*(pow(e,0.25)-pow(zetaSPeakEpsilon,0.25)),2) / (2.0*pow(zetaSSigmaPlus,2)));
 }
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
   else if (etaSparam == 2){
         return (std::max(0.0, etaSMin + ((e>etaSEpsilonMin) ? ( (ah*(e-etaSEpsilonMin)+aRho*nb) ): al*(e-etaSEpsilonMin)+aRho*nb))) * s;
   }else if (etaSparam == 3){
         if (abs(s) < 1e-10 || T < 1e-1) {
            return 0;
         }
         double Tc=T0+etaSShiftMuB*(mub/0.938);
         return ((e+p)/T)*std::max(0.0, (etaSMin + ((T>Tc) ? ah*(T-Tc) :  al*(T-Tc)))*(1 + etaSScaleMuB*(mub/0.938)));
   }
  return 0.;
}

void TransportCoeff::getTau(double e, double p, double nb, double s, double T, double mub, double &_taupi, double &_tauPi) {
 if (T > 1e-10 && abs(e + p) > 1.e-10) {
   //based on 1403.0962.  As no cs2 table is available, we approximate it by 0.228 (??)
  _taupi = (5./5.068)* eta(e, p, nb, s, T, mub) / (e + p);
  _tauPi =  (6./5.068)*zeta(e, p, nb, s, T, mub) / (e + p);
 } else {
  _taupi = _tauPi = 0.;
 }
}
