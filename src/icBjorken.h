
class EoS;
class Fluid;

// initial state from Gubser analytical solution (ideal hydro)
class ICBjorken {
  double epsilon;
 public:
  ICBjorken(double e);
  ~ICBjorken(void);
  // setIC: initializes entire hydro grid at a given initial proper time tau
  void setIC(Fluid *f, EoS *eos, double tau);
};
