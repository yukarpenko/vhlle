#include "eosLinearCombination.h"
#include "eoSimpleSpline.h"
//#include "eos.h"


EoSLinearCombination::EoSLinearCombination(EoSimpleSpline *eos1_, EoSimpleSpline *eos2_, double taus_, double tau0_) : eos1(eos1_), eos2(eos2_), taus(taus_), tau0(tau0_)
{
}


EoSLinearCombination::~EoSLinearCombination(void)
{
}

double EoSLinearCombination::p(double la, double T) {
	//return eos1->fPT(T);
	return la * eos1->fPT(T) + (1. - la) * eos2->fPT(T);
}

double EoSLinearCombination::en(double la, double T) {
	return la * eos1->fET(T) + (1. - la) * eos2->fET(T);
}

double EoSLinearCombination::dpe(double e, double tau) {
	double la = lambda(tau);
	double de = 0.001;
	return (p(la, t(e+de,tau)) - p(la, t(e,tau))) / de;
}

double EoSLinearCombination::lambda(double tau) {
	if (taus==0.) return 1.;
	if (tau==tau0) return 0.;
	return 1. - exp((tau0 - tau) / taus);
}

double EoSLinearCombination::t(double e, double la) {
	//return eos1->t(e, tau);
	double left = 1e-9, right = 10., center = 0.;
	//double la = lambda(tau);
	double valleft = en(la, left) - e, valright = en(la, right) - e, valcenter = 0.;
	if (valleft * valright > 0.) return left;
	while ((right-left)/right > 1.e-5) {
		center    = (left + right) / 2.;
		valcenter = en(la, center) - e;
		if (valcenter * valleft > 0) {
			left    = center;
			valleft = valcenter;
		}
		else {
			right    = center;
			valright = valcenter;
		}
	}
	return center;
}

double EoSLinearCombination::mu(double e) {
	return 0.;
}


double EoSLinearCombination::s(double e, double nb, double nq, double ns, double tau) {
	//return eos1->s(e, nb, nq, ns, tau);
	double la = lambda(tau);
	double T  = t(e, la);
	double ret = (e + p(la, T)) / T;
	if (T>1.e-5 && la>0.) ret += la * log(la) * (eos2->fPT(T) - eos1->fPT(T)) / T;
	return ret;
}
