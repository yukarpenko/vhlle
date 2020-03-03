class EoS;
class EoSSmash : public EoS {
private:
 double emax, nmax, emin, nmin;
 int ne, nn;
 double **ptab, **Ttab, **mubtab, **mustab;//, **stab;

public:
 EoSSmash(char* filename, int Ne, int Nn);
 ~EoSSmash(void);

 virtual void eos(double e, double nb, double nq, double ns, double &_T,
                  double &_mub, double &_muq, double &_mus, double &_p);
 virtual double p(double e, double nb, double nq, double ns);
};
