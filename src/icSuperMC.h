#include <vector>

class Fluid;
class EoS;

class ICSuperMC{
    private:
        std::string name_file;
        double tau0;
        // std::string setup;
        double sNN, eta0, sigmaeta, w, eff, etaB, sigmaIN, sigmaOUT;
        double const ZoverA = 0.4; //ratio between atomic and mass number (it's the same for Pb and Au more or less...)
        double const mN = 0.939; //GeV nucleon mass


        double T(double x, double y, std::vector<double> &xa, std::vector<double> &ya, double xcdm, double ycdm, double maxR);
        double yCM(double TA, double TB, double ybeam);
        double Mass(double TA, double TB, double ybeam, double mN);
        double normalization(double M, double C, double sigmaeta, double eta0);
        double C(double eta0, double sigmaeta);
        double energy(double N, double etas, double ycm, double eta0, double sigmaeta, double ybeam);
        double baryon_density_profile(double etas, double etaB, double sigmaIN, double sigmaOUT, int id);
        double max_distance_from_center(std::vector<double> x, std::vector<double> y, double xcdm, double ycdm);
        void initialize_superMC_parameters(std::string setup_data);


    public:
        ICSuperMC(std::string filename, double tau0, std::string setup_data);
        ~ICSuperMC();
        void setIC(Fluid *f, EoS *eos);
};
