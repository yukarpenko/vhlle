#include <iostream>
#include <iomanip>
#include "eos.h"
#include "trancoeff.h"
#include "inc.h"

void TransportCoeff::printZetaT()
{
 std::cout << "------zeta/s(T):\n";
 for(double e=0.1; e<3.0; e+=0.1){
  double T, mub, muq, mus, p;
  eos->eos(e, 0., 0., 0., T, mub, muq, mus, p);
  double s=eos->s(e, 0., 0., 0.);
  std::cout << std::setw(14) << T << std::setw(14) << zetaS(e, T, s, p) << std::endl;
 }
 std::cout << "---------------:\n";
}

double TransportCoeff::zetaS(double e, double T, double s, double P)
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

 if(zetaSparam == 0)          // arXiv:0908.3656
    return zetaS0 * (1. / 3. - eos->cs2(e)) / (exp((0.16 - T) / 0.001) + 1.);
 else if(zetaSparam == 1)     // arXiv:1704.04216
 {
    if(T < 0.180)
       return  zetaS0 * (0.03 + ( 0.08 * exp(((T / T_p) - 1.) / 0.0025)) + ( 0.22 * exp(((T / T_p) - 1) / 0.022)));
    else if( T >= 0.180 && T < 0.200)
       return  zetaS0 * (27.55 * (T / T_p) - 13.45 - (13.77 * (T / T_p) * (T / T_p)));
    else if(T >= 0.200)
       return  zetaS0 * (0.001 + (0.9 * exp(-((T / T_p) - 1.) / 0.0025)) + ( 0.25 * exp(-((T / T_p) - 1.) / 0.13)));
 }
 else if(zetaSparam == 2)     // arXiv:1901.04378
 {
    if(T > T_peak)
       return zetaS0 * (B_norm * ((B_width * B_width)/((((T / T_peak) - 1.) * ((T / T_peak) - 1.)) + (B_width * B_width))));
    else if(T <= T_peak)
       return zetaS0 * (B_norm * (exp(-((T - T_peak) / T_width) * ((T - T_peak) / T_width))));
 }
 else if(zetaSparam == 3)     // arXiv:2005.14682
 {
    if(T < T_peak2)
       return zetaS0 * (B_norm2 * exp(-((T - T_peak2) * (T - T_peak2)/(B1 * B1))));
    else if(T >= T_peak2)
       return zetaS0 * (B_norm2 * exp(-((T - T_peak2) * (T - T_peak2)/(B2 * B2))));
 }else if(zetaSparam == 4)
 {  
    if (abs(s) < 1e-10 || T < 1e-1) {
            return 0;
         }
    if ( e < zetaSPeakEpsilon)
      return ((e+P)/(s*T))*zetaS0 * exp(pow(zetaSScaleBeta*(pow(e,0.25)-pow(zetaSPeakEpsilon,0.25)),2) / 2.0*pow(zetaSSigmaMinus,2));
    else
      return ((e+P)/(s*T))*zetaS0 * exp(pow(zetaSScaleBeta*(pow(e,0.25)-pow(zetaSPeakEpsilon,0.25)),2) / 2.0*pow(zetaSSigmaPlus,2));
 }
}

void TransportCoeff::getEta(double e, double rho, double T, double muB, double s, double P, double &_etaS, double &_zetaS) {
 _etaS = etaS(e,rho, T, muB, s, P);
 _zetaS = zetaS(e,T, s, P);
}

double TransportCoeff::etaS(double e,double rho, double T, double muB, double s, double P)
{
   if (etaSparam == 0){
         return etaS0;
   }
   else if (etaSparam == 1){
         return etaSMin +  ((T>T0) ? ah*(T-T0) :  al*(T-T0));
   }
   else if (etaSparam == 2){
         return std::max(0.0, etaSMin + ((e>etaSEpsilonMin) ? ( (ah*(e-etaSEpsilonMin)+aRho*rho) ): al*(e-etaSEpsilonMin)+aRho*rho));
   }else if (etaSparam == 3){
         if (abs(s) < 1e-10 || T < 1e-1) {
            return 0;
         }
         double Tc=T0+etaSShiftMuB*(muB/0.938);
         return ((e+P)/(s*T))*std::max(0.0, (etaSMin + ((T>Tc) ? ah*(T-Tc) :  al*(T-Tc)))*(1 + etaSScaleMuB*(muB/0.938)));
   }
}

void TransportCoeff::getTau(double e, double rho, double T, double muB, double s, double P, double &_taupi, double &_tauPi) {
 if (T > 1e-1 && abs(e + P) > 1.e-10) {
   //based on 1403.0962.  As no cs2 table is available, we approximate it by 0.228
  _taupi = (5./5.068)* s * etaS(e,rho, T, muB, s, P) / (e + P);
  _tauPi =  (6./5.068)*zetaS(e,T, s, P) * s / (e + P);
 } else {
  _taupi = _tauPi = 0.;
 }
}
