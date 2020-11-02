class EoS;
class EoSSmash : public EoS {
  // Bounds (upper and lower) of energy density, baryon density and charge density
  double emax, emin, nbmax, nbmin, nqmax, nqmin;
  // Number of points in energy density, baryon density, charge density
  int ne, nnb, nnq;
  // Tabularised values for pressure, Temperature, muB, muQ, muS
  // Dimension: ne * nnb * nnq
  double *ptab, *Ttab, *mubtab, *muqtab, *mustab;

  // Get absolute index from given indices in energy density, baryon density
  // and charge density
  inline int index3(int ie, int inb, int inq) {
   return inq + nnq * inb + nnq * nnb * ie;
  }

public:
 EoSSmash(char* filename, int Ne, int Nnb, int Nq);
 ~EoSSmash(void);

 // Set (T, muB, muQ, muS, p) from (e, nb, nq, ns) according to EoS
 virtual void eos(double e, double nb, double nq, double ns, double &_T,
                  double &_mub, double &_muq, double &_mus, double &_p);
 // Calculate pressure from (e, nb, nq, ns)
 virtual double p(double e, double nb, double nq, double ns);
};
