class EoS ;

class EoSHadron : public EoS
{
  double e0, n0, logemin, logemax, lognmax ;
  int ne, nnb, nnq ;
  double nb_abs_min, nq_abs_min, e_min ;
  double *ptab, *Ttab, *mubtab, *muqtab, *mustab ;
  int *statustab ;
  inline int index3(int ie, int inb, int inq)
  { return inq+nnq*inb+nnq*nnb*ie ; }
public:
	EoSHadron(char *filename);
	~EoSHadron(void);

	virtual void eos(double e, double nb, double nq, double ns,
		double &_T, double &_mub, double &_muq, double &_mus, double &_p) ;
	virtual double p(double e, double nb, double nq, double ns) ;
};
