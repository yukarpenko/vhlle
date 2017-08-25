class EoS;
class EoSaux;
class EoSChiral : public EoS {
private:
 EoSaux *eossmall, *eosbig;

public:
 EoSChiral(void);
 ~EoSChiral(void);

 virtual void eos(double e, double nb, double nq, double ns, double &_T,
                  double &_mub, double &_muq, double &_mus, double &_p);
 virtual double p(double e, double nb, double nq, double ns);
};
