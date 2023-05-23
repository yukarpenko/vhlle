
class Fluid;
class EoS;
#include <vector>

class IcTable3d {
private:
 const int nx_grid=261, ny_grid=261, neta_grid=64;
 const double xminG=-13., xmaxG=13., yminG=-13., ymaxG=13., etaminG=-7., etamaxG=7.; // grid sizes
 //double ed_grid [nx_grid][ny_grid][neta_grid];
 std::vector<std::vector<std::vector<double>>> ed_grid;
 double tau0;

public:
 IcTable3d(Fluid *f, const char *filename, double tau0);
 ~IcTable3d();
 double interpolateGrid(double x, double y, double eta);
 void setIC(Fluid *f, EoS *eos);
};
