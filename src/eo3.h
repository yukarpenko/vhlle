
class EoS;

class EoS3f : public EoS {
 private:
  double emax, nmax, e0, n0;
  int ne, nn;
  double B, volex0, delta0, aaa, bbb;
  double *egrid, *ngrid, *lgegrid, *lgngrid;
  double *T, *pre, *mub, *muq, *mus;
  inline int index(int ie, int ib, int iq, int is) {
    return ie + ne * ib + ne * nn * iq + ne * nn * nn * is;
  }

 public:
  EoS3f(char *filename, double B, double volex0, double delta0, double aaa,
        double bbb);
  ~EoS3f(void);

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
