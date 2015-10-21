
class EoS;
class Fluid;
class TF1;

// this class takes care of the initial conditions for hydrodynamic evolution
class ICGlauber {
  TF1 *iff;  // instance of TF1 class used for numerical integration relevant to
             // the initial conditions from optical Glauber approach
  double rho0, prms[4];  // relevant to the initial conditions from optical
                         // Glauber approach
  double *_rphi;
  void findRPhi(void);
  double rPhi(double phi);

  double WoodSaxon(double *x, double *p);  // Wood-Saxon profile
  double Thickness(double *x, double *p);  // nuclear thickness profile function
                                           // for optical Glauber approach
  // epsilon: normalization of the initial energy density
  // alpha: parameter relevant to the initial transverse flow
  // b: impact parameter (for optical Glauber)
  double epsilon, b, tau0;

  // Nucleus parameters
  double A, Ra, dlt, sigma;
 public:
  ICGlauber(double e, double impactPar, double _tau0, double A_ = 208, double R_ = 6.5, double dlt_ = 0.54, double sigma_ = 4.);
  ~ICGlauber(void);
  // energy density profile at given point in transverse plane
  double eProfile(double x, double y);
  // Init params
  void init();
  // setIC: initializes entire hydro grid at a given initial proper time tau
  void setIC(Fluid *f, EoS *eos);
  double getrho0() const { return rho0; }
};
