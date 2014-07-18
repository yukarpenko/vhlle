class EoS ;

class TransportCoeff
{
 double etaS, zetaS, taupi, tauPi ;
 EoS *eos ;
 public:
 TransportCoeff(double _etaS, double _zetaS, EoS *_eos) ;
 ~TransportCoeff() {} ;
 void getEta(double e, double T, double &_etaS, double &_zetaS) ;
 void getTau(double T, double &_taupi, double &_tauPi) ;
 inline bool isViscous() { if(etaS>0. || zetaS>0.) return true ; else return false ; }
} ;
