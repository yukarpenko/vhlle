#include "SplineFunction.h"

class EoS;


class EoSimpleSpline :
	public EoS
{
	SplineFunction splPE;
	SplineFunction splTE;
public:
	EoSimpleSpline(std::string fname);
	~EoSimpleSpline(void);

	virtual inline void eos(double e, double nb, double nq, double ns, double &T,
                          double &mub, double &muq, double &mus, double &_p, double tau = 1.) {
    _p = p(e);
    T = t(e);
    mub = muq = mus = 0.;
  }
  virtual inline double p(double e, double nb, double ns, double nq) {
    return p(e);
  }

  double p(double e);
  double dpe(double e);
  double t(double e);
  double mu(double e);

  virtual double cs2(double e) { return dpe(e); }
};

