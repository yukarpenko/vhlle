namespace icEpos {

extern double tau0, dt;

void loadIC(const char* fnParams, const char* fnGrid);
void getIC(double x, double y, double eta, double &e, double &nb, double &nq,
           double &vx, double &vy, double &vz);
void setIC(Fluid *f, EoS *eos);
void deleteIC(void);
}
