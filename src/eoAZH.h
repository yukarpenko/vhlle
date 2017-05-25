class EoS;
class EoSAZHaux;
class EoSAZH : public EoS {
private:
 EoSAZHaux *p1, *p2, *T1, *T2, *mu1, *mu2;

public:
 EoSAZH(void);
 ~EoSAZH(void);

 virtual void eos(double e, double nb, double nq, double ns, double &_T,
                  double &_mub, double &_muq, double &_mus, double &_p);
 virtual double p(double e, double nb, double nq, double ns);
};
