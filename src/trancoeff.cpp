#include <iostream>
#include <iomanip>
#include "eos.h"
#include "trancoeff.h"
#include "inc.h"

TransportCoeff::TransportCoeff(double _etaS, double _zetaS, EoS *_eos) {
 etaS = _etaS;
 zetaS0 = _zetaS;
 eos = _eos;
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
 int param=1; //0 - basic ,1 ,2 - arxiv:1910.12930, 3 arxiv:2103.09848

 double T_p=0.180;

 double T_peak=0.165;
 double T_width=0.010;
 double B_norm=0.24;
 double B_width=1.5;
 double T_peak2=0.160;
 double B_norm2=0.13;
 double B1=0.01;
 double B2=0.12;
 
 if(param==0)
    return zetaS0 * (1. / 3. - eos->cs2(e)) / (exp((0.16 - T) / 0.001) + 1.);
 else if(param==1)
 {
    if(T<0.180)
       return 0.03+(0.08*exp(((T/T_p)-1.)/(0.0025)))+(0.22*exp(((T/T_p)-1)/(0.0022)));
    else if(T>=0.180 && T<0.200)
       return 27.55*(T/T_p)-13.45-(13.77*(T/T_p)*(T/T_p));
    else if(T>=0.200)
       return 0.001+(0.9*exp(-((T/T_p)-1.)/(0.0025)))+(0.25*exp(-((T/T_p)-1.)/(0.13)));
 }
 else if(param==2)
 {
    if(T>T_peak)
       return B_norm*((B_width*B_width)/((((T/T_peak)-1.)*((T/T_peak)-1.))+(B_width*B_width)));
    else if(T<=T_peak)
       return B_norm*(exp(-((T-T_peak)/T_width)*((T-T_peak)/T_width)));
 }
 else if(param==3)
 {
    if(T<T_peak2)
       return B_norm2*exp(-((T-T_peak2)*(T-T_peak2)/(B1*B1)));
    else if(T>=T_peak2)
       return B_norm2*exp(-((T-T_peak2)*(T-T_peak2)/(B2*B2)));
 }
}


void TransportCoeff::getEta(double e, double T, double &_etaS, double &_zetaS) {
 _etaS = etaS;
 _zetaS = zetaS(e,T);
}

void TransportCoeff::getTau(double e, double T, double &_taupi, double &_tauPi) {
 if (T > 0.) {
  _taupi = 5. / 5.068 * etaS / T;
  _tauPi = 6.0 / 5.068 * zetaS(e,T) / T;;
 } else {
  _taupi = _tauPi = 0.;
 }
}
