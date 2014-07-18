
class EoS ;
class Fluid ;
class TF1 ;

class IC
{
 TF1 *iff ;
 double rho0, prms[4] ;
 double *_rphi ;
 void findRPhi(void) ;
 double rPhi(double phi) ;
 
  double WoodSaxon(double *x, double *p);
  double Thickness(double *x, double *p);  
double epsilon, alpha, b ;
public:
	IC(double e, double impactPar, double a);
	~IC(void);
	double eProfile(double x, double y) ;
	void setIC(Fluid *f, EoS* eos, double tau) ;
};
