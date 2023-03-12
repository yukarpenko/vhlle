
class Fluid;
class EoS;

class IcTrento3d {
private:
 double A; // initial shear flow constant
 int nx, ny, nz;
 int n_grid,n_grid_eta;
 double xminG, xmaxG, yminG, ymaxG, etaminG, etamaxG; // grid sizes
 double xmin, xmax, ymin, ymax, zmin, zmax;
 double dx, dy, dz;
 double ***rho;
 double ***source;
 static const int NP = 10000;  // dimension for particle arrays
 double tau0;
 void makeSmoothTable(int npart);

public:
 IcTrento3d(Fluid *f, const char *filename, double tau0, const char* setup);
 ~IcTrento3d();
 double interpolateGrid(double x, double y,double eta);
 void setIC(Fluid *f, EoS *eos);
};
