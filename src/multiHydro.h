
class Fluid;
class Hydro;
class EoS;
class TransportCoeff;
class CrossSections;
class Cornelius;

class MultiHydro {
 Fluid *f_p, *f_t, *f_f; // projectile and target fluids
 Hydro *h_p, *h_t, *h_f; // solutions of hydro eqs. for projectile and target
 EoS *eos;
 Cornelius *cornelius;
 TransportCoeff *trcoeff;
 CrossSections *xsect;
 std::ofstream fmhfreeze;
 double ***MHeps, ***MHepsPrev;
 double ecrit, vEff;
 int nx, ny, nz;
 double dx, dy, dz, dtau;
 const double gmunu[4][4] = {
     {1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, -1, 0}, {0, 0, 0, -1}};

public:
 MultiHydro(Fluid *f_p, Fluid *f_t, Fluid *f_f, Hydro *h_p, Hydro *h_t,
  Hydro* h_f, EoS *eos, TransportCoeff *trcoeff, double dtau, double eCrit);
 ~MultiHydro(void);
 void initOutput(const char *dir);
 void performStep();
 void frictionSubstep();
 double** getEnergyMomentumTensor(double Q_p[7], double Q_f[7], double Q_t[7]);
 void getEnergyDensity();
 void updateEnergyDensity();
 void findFreezeout();
};
