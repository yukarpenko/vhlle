#ifndef EOSLINEARCOMBINATION_H
#define EOSLINEARCOMBINATION_H
#include "SplineFunction.h"
#include "eoSimpleSpline.h"
#include "eos.h"

//class EoS;
//class EoSimpleSpline;


class EoSLinearCombination :
	public EoS
{
	EoSimpleSpline *eos1;
	EoSimpleSpline *eos2;
	double ctau;
	double taus, tau0;
public:
	EoSLinearCombination(EoSimpleSpline *eos1_, EoSimpleSpline *eos2_, double taus_ = 0., double tau0_ = 0.1);
	virtual ~EoSLinearCombination(void);

	virtual inline void eos(double e, double nb, double nq, double ns, double &T,
                          double &mub, double &muq, double &mus, double &_p, double tau = 1.) {
	 ctau = tau;
	 double la = lambda(tau);
	 T = t(e, la);
     _p = p(la, T);
     mub = muq = mus = 0.;
    }
    virtual inline double p(double e, double nb, double ns, double nq, double tau = 1.) {
	  double la = lambda(tau);
      return p(la, t(e, la));
    }

	double lambda(double tau);

    double p(double la, double T);
	double en(double la, double T);
    double dpe(double e, double tau);
    double t(double e, double la);
    double mu(double e);

	virtual double s(double e, double nb, double nq, double ns, double tau = 1.);

    virtual double cs2(double e) { return dpe(e, ctau); }
};

#endif