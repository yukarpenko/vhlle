#include <cmath>
//#include "MersenneTwister.h"
#include "TRandom3.h"
#include "inc.h"

class EoS;
class Fluid;

class ICGlauberMC
{
	double cR;
	double ca;
	double rmax;
	double fmax;
	//MTRand ran;
	TRandom3 ran;

	// epsilon: normalization of the initial energy density
    // alpha: parameter relevant to the initial transverse flow
    // b: impact parameter (for optical Glauber)
    double epsilon, b, tau0;

	// Nucleus parameters
	int A;
	double Ra, dlt, sigma;
public:
	ICGlauberMC(double e, double impactPar, double _tau0, int A_ = 208, double R_ = 6.5, double dlt_ = 0.54, double sigma_ = 4.);
	~ICGlauberMC(void);
	double WS_r(double r, double R, double a)				// Woods-Saxon r density distribution function (unnormalized)
	{
		return r * r / ( 1 + exp( (r - R) / a) );
	}
	double WS_root_function(double r, double R, double a)   // This function root corresponds to rmax in WS distribution
	{
		return 2. + exp( (r - R) / a ) * ( 2. - r / a);
	}
	double Gauss(double x, double y, double xa, double ya, double w) const {
		return 1. / 2. / C_PI / w / w * exp(-(x-xa)*(x-xa)/2./w/w - (y-ya)*(y-ya)/2./w/w);
	}
	void WS_initparameters(double R, double a);						// Computates neccesary parameters for WS random number generation
	double Random_WS_r(double R, double a);							// Generates WS r
	// setIC: initializes entire hydro grid at a given initial proper time tau
    void setIC(Fluid *f, EoS *eos);
};

