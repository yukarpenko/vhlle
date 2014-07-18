#pragma once
#include <cmath>
#include <string>

class TGraph ;

class EoS
{
public :
	virtual ~EoS() { return ; } // why did I do this?
	virtual void eos(double e, double nb, double nq, double ns,
		double &T, double &mub, double &muq, double &mus, double &p) = 0 ;
	virtual double p(double e, double nb, double ns, double nq) = 0;
	double s(double e, double nb, double nq, double ns) ;
	inline double cs2(void) { return 1./3. ; }
	inline double cs(void) { return sqrt(1./3.) ; }
	virtual inline double cs2(double e) { return 1./3. ; };
};

class EoSs : public EoS
{
private:
 TGraph *gp, *gT, *gmu ;
public:
	EoSs(std::string fname, int ncols);
	~EoSs();

	virtual inline void eos(double e, double nb, double nq, double ns,
		double &T, double &mub, double &muq, double &mus, double &_p){
			_p = p(e) ;
			T = t(e) ;
			mub = muq = mus = 0. ;
	}
	virtual inline double p(double e, double nb, double ns, double nq)
	{ return p(e) ; }

	double p(double e) ;
	double dpe(double e) ;
	double t(double e) ;
	double mu(double e) ;

	virtual double cs2(double e) { return dpe(e) ; }
	//virtual double cs(double e) { return sqrt(dpe(e)) ; }

};
