#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "eoSimpleSpline.h"
#include "inc.h"

EoSimpleSpline::EoSimpleSpline(std::string fname, double taus_, double tau0_) : taus(taus_), tau0(tau0_)
{
	ifstream finput(fname.c_str(), ios::in);

	std::vector<double> evec(0), pvec(0), Tvec(0);

	if (!finput) {
		cerr << "can't open input file \"" << fname.c_str() << "\"" << endl;
		exit(1); 
	}

	double ein, pin, Tin;
	while (!finput.eof()) {
		finput >> ein >> pin >> Tin;
		evec.push_back(ein);
		pvec.push_back(pin);
		Tvec.push_back(Tin);
	}
	finput.close();

	splPE.fill(evec, pvec);
	splTE.fill(evec, Tvec);

	splPT.fill(Tvec, pvec);
	splET.fill(Tvec, evec);

	sigq  = 7.0 * C_PI * C_PI * 3.0 / 20.0;
	siggl = 8.0 * C_PI * C_PI / 15.0; 
}


EoSimpleSpline::~EoSimpleSpline(void)
{
}

double EoSimpleSpline::p(double e) {
	return splPE.f(e);
}

double EoSimpleSpline::dpe(double e) {
	return splPE.df(e);
}

double EoSimpleSpline::t(double e, double tau) {
	double ret = splTE.f(e);
	if (taus>0.0) {
		double lambda = 1. - exp((tau0-tau)/taus);
		ret *= pow((sigq + siggl) / (sigq*lambda + siggl), 0.25);
	}
	return ret;
}

double EoSimpleSpline::mu(double e) {
	return 0.;
}
