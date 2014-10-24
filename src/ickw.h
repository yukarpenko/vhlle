#include <iosfwd>

class Fluid;
class EoS;

class IC_KW {
 private:
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double dx, dy, dz;
  int nheadlines;  // number of header lines in input file
  char **header;   // header lines in input file

 public:
  int nx, ny, nz;
  double ***e;
  double ***vx;
  double ***vy;
  double ***vz;
  double ***nu;
  double ***nd;
  double ***ns;
  IC_KW(const char *filename);
  ~IC_KW(void);

  void writeHeader(std::ofstream &fout);

  void getICs(double x, double y, double z, double &e, double &vx, double &vy,
              double &vz, double &_nu, double &_nd, double &_ns);
  void setIC(Fluid *f, EoS *eos, double tau);
};
