
class Fluid;
class EoS;

class IcPartUrqmd {
private:
 int nx, ny, nz, nevents;
 double xmin, xmax, ymin, ymax, zmin, zmax;
 double dx, dy, dz;
 double ***T00, ***T0x, ***T0y, ***T0z, ***QB, ***QE;
 static const int NP = 10000;  // dimension for particle arrays
 // auxiliary particle arrays
 double Tau[NP], X[NP], Y[NP], Eta[NP], Mt[NP], Px[NP], Py[NP], Rap[NP];
 int Id[NP], Charge[NP];

 double tau0;
 double Rgx, Rgy, Rgz;
 int nsmoothx;  // smoothly distribute to +- this many cells
 int nsmoothy;
 int nsmoothz;
 void makeSmoothTable(int npart);

public:
 IcPartUrqmd(Fluid *f, const char *filename, double _Rgt, double _Rgz, double tau0);
 ~IcPartUrqmd();
 void setIC(Fluid *f, EoS *eos);
};
