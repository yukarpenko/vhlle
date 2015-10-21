#include "eos.h"
#include "eosuQGP.h"
#include "inc.h"


EoSuQGP::EoSuQGP(double taus_, double tau0_):taus(taus_), tau0(tau0_)
{
	sigq  = 7.0 * C_PI * C_PI * 3.0 / 20.0;
	siggl = 8.0 * C_PI * C_PI / 15.0; 
}


EoSuQGP::~EoSuQGP(void)
{
}

double EoSuQGP::t(double e, double tau) {
	double la = 1. - exp((tau0-tau)/taus);
	if (taus==0.) la = 1.;
	return e > 0. ? 1.0 * pow(e * pow(0.197326968, 3) / (sigq*la + siggl), 0.25) : 0.;
}

double EoSuQGP::s(double e, double nb, double nq, double ns, double tau) {
	double la   = 1. - exp((tau0-tau)/taus);
	if (taus==0.) la = 1.;
	double Temp = t(e, tau);
	if (e <= 0.) return 0.;
	double nqqbar = 32.0 * C_PI * C_PI / 45.0 * Temp * Temp * Temp * 0.156 * 3.0 * la / pow(0.197326968, 3);
	double ret = 4.0 * e / 3.0 / Temp;
	if (la > 0.) ret += - nqqbar * log(la);
	return ret;
}