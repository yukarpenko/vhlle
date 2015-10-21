class EoS;

class EoS1f : public EoS {
 private:
  double emax, nmax, e0, n0;
  int ne, nn;
  double *egrid, *ngrid, *lgegrid, *lgngrid;
  double *T, *pre, *mub;
  inline int index(int ie, int ib) { return ie + ne * ib; }

 public:
  EoS1f(char *filename);
  ~EoS1f(void);

  void eosranges(double &_emax, double &_e0, double &_nmax, double &_n0,
                 int &_ne, int &_nn);
  void getue(double e, int &ixe, double &ue);
  void getun(double n, int &ixn, double &un);
  virtual void eos(double e, double nb, double nq, double ns, double &_T,
                   double &_mub, double &_muq, double &_mus, double &_p, double tau = 1.);
  virtual inline double p(double e, double nb, double nq, double ns) {
    double T, mub, muq, mus, pp;
    eos(e, nb, nq, ns, T, mub, muq, mus, pp);
    return pp;
  }
};
