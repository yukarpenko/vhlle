
class EoS;

class EoSuQGP :
	public EoS
{
	double taus, tau0;
	double sigq, siggl;
public:
	EoSuQGP(double taus_, double tau0_);
	~EoSuQGP(void);

	virtual inline void eos(double e, double nb, double nq, double ns, double &T,
                        double &mub, double &muq, double &mus, double &_p, double tau = 1.e9) {
		_p = e / 3.;
		T = t(e, tau);
		mub = muq = mus = 0.;
	}

  virtual inline double p(double e, double nb, double ns, double nq) {
    return e / 3.;
  }

  double t(double e, double tau = 1.e9);
  double mu(double e) { return 0.; }

  virtual double s(double e, double nb, double nq, double ns, double tau = 1.e9);

  virtual double cs2(double e) { return 1./3.; }
};

