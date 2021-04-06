
class EoS;
class Fluid;

// this class takes care of the initial conditions for hydrodynamic evolution
class IC3F {
private:
 int nx, ny, nz;
 int nevents;
 double xmin, xmax, ymin, ymax, zmin, zmax;
 double dx, dy, dz;
 double ***T00_p, ***T0z_p, ***QB_p, ***QE_p;
 double ***T00_t, ***T0z_t, ***QB_t, ***QE_t;
 double snn;
 int projA;
 int targA;
 int projZ;
 int targZ;
 double WSdelta = 0.5; // diffuseness in Wood-Saxon
 const double nucleon_mass = 0.938;
 double tau0;

 double Rgx, Rgy, Rgz;
 int nsmoothx;  // smoothly distribute to +- this many cells
 int nsmoothy;
 int nsmoothz;
 void makeSmoothPart(double x, double y, double eta, int Charge, double rap, bool isProjectile);

public:
 IC3F(Fluid *f_p, Fluid *f_t, double tau, int _nevents, double _snn, int _projA, int _targA, int _projZ, int _targZ, double _Rg);
 ~IC3F(void);
 // setIC: initializes entire hydro grid at a given initial proper time tau
 void setIC(Fluid* f_p, Fluid *f_t, EoS* eos);
};
