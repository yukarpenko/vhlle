
class EoS;
class Fluid;

struct WS_params {   // parameters of the Woods-Saxon profile
    double r, norm, Ra, delta;
};

// this class takes care of the initial conditions for hydrodynamic evolution
class ICGlauber {
 double rho0;  // relevant to the initial conditions from optical
                        // Glauber approach
 double *_rphi;
 void findRPhi(WS_params* params);
 double rPhi(double phi);

 //double WoodSaxon(double x, void *params);  // Wood-Saxon profile
 //double Thickness(double x, double y, void *params);  // nuclear thickness profile function
                                          // for optical Glauber approach
 // epsilon: normalization of the initial energy density
 // alpha: parameter relevant to the initial transverse flow
 // b: impact parameter (for optical Glauber)
 double epsilon, b, tau0;

public:
 ICGlauber(double e, double impactPar, double _tau0);
 ~ICGlauber(void);
 // energy density profile at given point in transverse plane
 double eProfile(double x, double y, WS_params* params);
 // setIC: initializes entire hydro grid at a given initial proper time tau
 void setIC(Fluid *f, EoS *eos);
};
