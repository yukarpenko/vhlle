
class Cell;
class Fluid;
class EoS;
class TransportCoeff;

class Hydro {
 private:
  Fluid *f;
  EoS *eos;
  TransportCoeff *trcoeff;
  double dt, time;

 public:
  Hydro(Fluid *_f, EoS *_eos, TransportCoeff *_trcoeff, double _t0, double _dt);
  ~Hydro();
  void setDt(double _dt);

  void hlle_flux(Cell *left, Cell *right, int direction, int mode);
  void visc_flux(Cell *left, Cell *right, int direction);
  void NSquant(int ix, int iy, int iz, double pi[][4], double &Pi,
               double dmu[4][4], double &du);
  void setNSvalues();
  void ISformal();
  void setQfull();
  void performStep(void);
  inline double getTime() { return time; }
};
