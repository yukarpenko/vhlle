
class Cell;
class Fluid;
class EoS;
class TransportCoeff;

// this class implements the hydrodynamic evolution and
// contains the hydrodynamic algorithm
class Hydro {
 private:
  Fluid *f;
  EoS *eos;
  TransportCoeff *trcoeff;
  double dt, time;

 public:
  Hydro(Fluid *_f, EoS *_eos, TransportCoeff *_trcoeff, double _t0, double _dt);
  ~Hydro();
  void setDt(double _dt);  // change the timestep

  // HLLE (ideal)flux between two neighbouring cells in a given direction
  // mode: PREDICT = used in predictor step; calculates fluxes for dt/2
  // CORRECT = used in corrector step, calculates fluxes based on predicted
  // half-step quantities
  void hlle_flux(Cell *left, Cell *right, int direction, int mode);
  // viscous flux \delta F
  void visc_flux(Cell *left, Cell *right, int direction);
  // shear stress tensor and bulk pressure in Navier-Stokes (NS) limit
  // plus \partial_\mu u^\nu matrix (dmu) and
  // expansion scalar \partial_mu u^\mu (du)
  // for a given cell (ix,iy,iz)
  void NSquant(int ix, int iy, int iz, double pi[][4], double &Pi,
               double dmu[4][4], double &du);
  // sets the values of shear stress/bulk pressure in NS limit in all hydro grid
  void setNSvalues();
  // advances numerical solution for shear/bulk in a whole grid over one
  // timestep
  void ISformal();
  void setQfull();
  // advances numerical solution for Q (including ideal and viscous fluxes and
  // source terms) over one timestep
  void performStep(void);
  inline double getTime() { return time; }
};
