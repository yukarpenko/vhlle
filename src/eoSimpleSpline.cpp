#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "eos.h"
#include "eoSimpleSpline.h"

EoSimpleSpline::EoSimpleSpline(std::string fname)
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

double EoSimpleSpline::t(double e) {
	return splTE.f(e);
}

double EoSimpleSpline::mu(double e) {
	return 0.;
}