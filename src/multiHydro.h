
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
 double ecrit;
 int nx, ny, nz;

public:
 MultiHydro(Fluid *f_p, Fluid *f_t, Fluid *f_f, Hydro *h_p, Hydro *h_t,
  Hydro* h_f, EoS *eos, TransportCoeff *trcoeff, double dtau, double eCrit);
 ~MultiHydro(void);
 void initOutput(const char *dir);
 void performStep();
 void frictionSubstep();
 void getEnergyDensity();
 void updateEnergyDensity();
 void findFreezeout();
};
