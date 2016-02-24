#ifndef eossimplespline_h
#define eossimplespline_h

#include "SplineFunction.h"
#include "eos.h"

//class EoS;


class EoSimpleSpline :
	public EoS
{
	double taus, tau0;
	double sigq, siggl;
public:
	SplineFunction splPE;
	SplineFunction splTE;
	SplineFunction splET;
	SplineFunction splPT;
	EoSimpleSpline() { }
	EoSimpleSpline(std::string fname, double taus_ = 0., double tau0_ = 0.1);
	virtual ~EoSimpleSpline(void);

	virtual inline void eos(double e, double nb, double nq, double ns, double &T,
                          double &mub, double &muq, double &mus, double &_p, double tau = 1.) {
    _p = p(e);
    T = t(e, tau);
    mub = muq = mus = 0.;
  }
  virtual inline double p(double e, double nb, double ns, double nq, double tau = 1.) {
    return p(e);
  }

  double p(double e);
  double dpe(double e);
  double t(double e, double tau);
  double mu(double e);

  virtual double cs2(double e) { return dpe(e); }

  double fPT(double T) { return splPT.f(T); }
  double fET(double T) { return splET.f(T); }
  double fTE(double E) { return splTE.f(E); }

};


#endif
