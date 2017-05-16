
class Fluid;
class EoS;

class IcGlissando {
private:
 int nx, ny, nz, nevents;
 double xmin, xmax, ymin, ymax, zmin, zmax;
 double dx, dy, dz;
 double ***rho;
 static const int NP = 10000;  // dimension for particle arrays
 // auxiliary particle arrays
 double X[NP], Y[NP], W[NP];
 int C[NP];

 double tau0;
 int nsmoothx;  // smoothly distribute to +- this many cells
 int nsmoothy;
 int nsmoothz;
 void makeSmoothTable(int npart);

public:
 IcGlissando(Fluid *f, char *filename, double tau0);
 ~IcGlissando();
 void setIC(Fluid *f, EoS *eos);
};
