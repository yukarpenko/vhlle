class EoS;
class EoSCMFeaux;
class EoSCMFe : public EoS {
private:
 EoSCMFeaux *eossmall, *eosbig;

public:
 EoSCMFe(void);
 ~EoSCMFe(void);
 const double EPS_0 = 146.51751415742*0.87326281;
  const double N_0 = 0.15891*0.87272727;

 virtual void eos(double e, double nb, double nq, double ns, double &_T,
                  double &_mub, double &_muq, double &_mus, double &_p);
 virtual double p(double e, double nb, double nq, double ns);
};
