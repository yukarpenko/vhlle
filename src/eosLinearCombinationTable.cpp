#include "eosLinearCombinationTable.h"
#include "eoSimpleSpline.h"
//#include "eos.h"

#include <iostream>


EoSLinearCombinationTable::EoSLinearCombinationTable(EoSimpleSpline *eos1_, EoSimpleSpline *eos2_, double taus_, double tau0_) : taus(taus_), tau0(tau0_)
{
	eoslin = new EoSLinearCombination(eos1_, eos2_, taus_, tau0_);
	dlam = 0.01;
	eosla.resize(0);
	for(double la=0.;la<1.+1e-5;la+=dlam) {
		EoSimpleSpline teos = *eos2_;
		for(int i=0;i<teos.splPE.vals.size();++i) {
			
			double Temp = teos.splPT.vals[i].first;
			double tP = eoslin->p(la, Temp);
			double tE = eoslin->en(la, Temp);
			teos.splPT.vals[i].second = tP;
			teos.splPE.vals[i].second = tP;
			teos.splPE.vals[i].first  = tE;
			//teos.splTE.vals[i].second = tP;
			teos.splTE.vals[i].first  = tE;
			teos.splET.vals[i].second = tE;
		}
		eosla.push_back(teos);
		//eosla.push_back(EoSimpleSpline
	}
}


EoSLinearCombinationTable::~EoSLinearCombinationTable(void)
{
	delete eoslin;
}

double EoSLinearCombinationTable::p(double la, double T) {
	//return eos1->fPT(T);
	//return la * eos1->fPT(T) + (1. - la) * eos2->fPT(T);
	if (la<=0.) return eosla[0].fPT(T);
	if (la>=1.) return eosla[eosla.size()-1].fPT(T);
	int ind = index(la);
	if (ind==eosla.size()-1) return eosla[ind].fPT(T);
	else {
		double dist = la - ind*dlam;
		return (1. - dist / dlam) * eosla[ind].fPT(T) + (dist / dlam) * eosla[ind+1].fPT(T);
	}
}

double EoSLinearCombinationTable::en(double la, double T) {
	int ind = index(la);
	if (la<=0.) return eosla[0].fET(T);
	if (la>=1.) return eosla[eosla.size()-1].fET(T);
	if (ind==eosla.size()-1) return eosla[ind].fET(T);
	else {
		double dist = la - ind*dlam;
		return (1. - dist / dlam) * eosla[ind].fET(T) + (dist / dlam) * eosla[ind+1].fET(T);
	}
}

double EoSLinearCombinationTable::dpe(double e, double tau) {
	double la = lambda(tau);
	double de = 0.001;
	return (p(la, t(e+de,tau)) - p(la, t(e,tau))) / de;
}

double EoSLinearCombinationTable::lambda(double tau) {
	if (taus==0.) return 1.;
	if (tau==tau0) return 0.;
	return 1. - exp((tau0 - tau) / taus);
}

double EoSLinearCombinationTable::t(double e, double la) {
	int ind = index(la);
	if (la<=0.) return eosla[0].fTE(e);
	if (la>=1.) return eosla[eosla.size()-1].fTE(e);
	if (ind==eosla.size()-1) return eosla[ind].fTE(e);
	else {
		double dist = la - ind*dlam;
		return (1. - dist / dlam) * eosla[ind].fTE(e) + (dist / dlam) * eosla[ind+1].fTE(e);
	}
}

double EoSLinearCombinationTable::mu(double e) {
	return 0.;
}

//#include <iostream>

double EoSLinearCombinationTable::s(double e, double nb, double nq, double ns, double tau) {
	//return eos1->s(e, nb, nq, ns, tau);
	double la = lambda(tau);
	double T  = t(e, la);
	double ret = (e + p(la, T)) / T;
	if (T>1.e-5 && la>0.) ret += la * log(la) * (eosla[0].fPT(T) - eosla[eosla.size()-1].fPT(T)) / T;
	if (T<1.e-9) ret = 0.;
	if (ret!=ret) std::cout << "AAAAAA  " << t(e, la) << " " << e << " " << p(la, T) << " " << (eosla[0].fPT(T) - eosla[eosla.size()-1].fPT(T)) / T << "\n";
	//if (ret!=ret) ret = 0.;
	return ret;
}
