class EoS;
class EoSaux;
class EoSCMF : public EoS {
public:
 EoSCMF(void);
 EoSCMF(std::string filename, int N_T, int N_nB);
 ~EoSCMF(void);
    double Tmax, nmax, Tmin, nmin;
    double dN, dT;
    int N_T, N_nB;
    double **ptab, **Ttab, **etab, **mubtab, **mustab, **stab, **qfractab;
    double EPS_0, N_0;

 virtual void eos(double e, double nb, double nq, double ns, double &_T,
                  double &_mub, double &_muq, double &_mus, double &_p);
 virtual double p(double e, double nb, double nq, double ns);

private:
//    EoSaux *eos_aux;
    int binary_search(double *a, int first, int last, double search_num);
    int get_ind(double eps, double n);
    void get(double e, double nb, double& p, double& T, double& mub, double& mus);
protected:
    double get_temp(double eps, double n);
};
