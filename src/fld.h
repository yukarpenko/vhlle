#include <iosfwd>
#include "cll.h"

class EoS;
class TransportCoeff;
class Cornelius;

#ifdef _MSC_VER

bool isinf(const double &op);

bool isnan(const double &op);

#endif

#include <vector>

// This structure contains output of physical quantities at a given tau timestep
// Quantities: e, s, p, T, nb, nq, ns, taup, vx, vy
struct TauOutput {
	std::vector< std::vector<double> > XOut;	// At Y=Eta=0
	std::vector< std::vector<double> > YOut;	// At X=Eta=0
	std::vector< std::vector<double> > XYOut;	// At Eta=0
	std::vector< std::vector<double> > EtaOut;	// At X=Y=0
	double E, Efull, S, dSdEta;
};

// this class contains the information and methods related to the hydro grid
// 'z' direction actually denotes eta direction, as well as
// dz == d(eta), getZ returns eta coordinate etc.
class Fluid {
 private:
  EoS *eos, *eosH;          // equation(s) of state
  TransportCoeff *trcoeff;  // transport coefficients for visc fluid
  Cornelius *cornelius;  // instance of Cornelius to calculate the hypersurface
  Cell *cell;            // 3D hydro grid, packed in 1D array
  Cell *cell0;           // reference to cell containing all zero quantities
  int nx, ny, nz;        // dimensions of the grid
                         // physical dimensions of the grid
  double minx, maxx, miny, maxy, minz, maxz;
  double dx, dy, dz, dt; // physical sizes of the hydro cell and timestep
  double ecrit;
  double vEff, EtotSurf;  // cumulative effective volume and energy flux
  double SSurf;		// cumulative entropy flux
  double SSurfdeta;	// cumulative dS/dEta
  bool SurfCriterion; // 0 - eCrit, 1 - TCrit
  std::ofstream foutkw, foutkw_dim, foutxvisc, foutyvisc, foutdiagvisc, foutx,
      fouty, foutdiag, foutz, fout_aniz, fout2d, ffreeze;
  int compress2dOut;

  std::fstream foutXTau;
  std::fstream foutYTau;
  std::fstream foutTau;

 public:
  Fluid(EoS *_eos, EoS *_eosH, TransportCoeff *_trcoeff, int _nx, int _ny,
        int _nz, double _minx, double _maxx, double _miny, double _maxy,
        double _minz, double _maxz, double dt, double eCrit, bool SurfCrit = 0);
  ~Fluid();
  void initOutput(char *dir, int maxstep, double tau0, int cmpr2dOut);
  int output_xy_spacing;
  int output_eta_points;
  int output_tau_spacing;
  int output_eta_spacing;

  int output_nt, output_nx, output_ny;

  inline int getNX() { return nx; }
  inline int getNY() { return ny; }
  inline int getNZ() { return nz; }
  inline double getDx() { return dx; }
  inline double getDy() { return dy; }
  inline double getDz() { return dz; }
  inline double getX(int ix) { return minx + ix * dx; }
  inline double getY(int iy) { return miny + iy * dy; }
  inline double getZ(int iz) { return minz + iz * dz; }

  void getCMFvariables(Cell *c, double tau, double &e, double &nb, double &nq,
                       double &ns, double &vx, double &vy, double &Y);

  inline Cell *getCell(int ix, int iy, int iz) {
    ix = ix > 0 ? ix : 0;
    ix = ix < nx ? ix : nx - 1;
    iy = iy > 0 ? iy : 0;
    iy = iy < ny ? iy : ny - 1;
    iz = iz > 0 ? iz : 0;
    iz = iz < nz ? iz : nz - 1;
    return &cell[ix + nx * iy + nx * ny * iz];
  }

  void correctImagCells(void);      // only ideal hydro part
  void correctImagCellsFull(void);  // correct ideal+visc
  void updateM(double tau, double dt);

  void outputGnuplot(double tau);
  void outputSurface(double tau);
  void outputCorona(double tau);

  void calcTotals(double tau);

  void outputXTau(double tau);
  void outputYTau(double tau);
  void outputTau(double tau);

  //TauOutput output(double tau);
  void output(double tau, TauOutput & outp);
};
