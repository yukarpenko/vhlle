#include <vector>

class Fluid;
class EoS;
class Cell;

struct DileptonSpectrumEntry {
	double M, y;
	double yield;//, v1, v2, v3;
	double v2;
};

struct DileptonSpectrum {
	std::vector<DileptonSpectrumEntry> Entries;
};

// Base class for calculating the yield of thermal dileptons
class Dileptons
{
protected:	
	std::vector<double> Ms;     // Set of pt values to calculate spectrum at
	std::vector<double> ys;      // Set of y  values to calculate spectrum at
	std::vector<double> yield;   // Dilepton yields at corresponding y and M
	//std::vector<double> v1;      // Dilepton directed   flow at corresponding y and pt (unnormalized)
	std::vector<double> v2;      // Dilepton elliptic   flow at corresponding y and pt (unnormalized)
	//std::vector<double> v3;      // Dilepton triangular flow at corresponding y and pt (unnormalized)
	double Tcut;				 // Cutoff temperature
public:
	Dileptons(double Tcut_ = 0.155);
	Dileptons(const std::vector<double> & Min, const std::vector<double> & yin, double Tcut_ = 0.155);
	virtual ~Dileptons();
	void addPtY(double pt, double y);
	virtual void addCell(double tau, double dtau, Fluid *f, EoS *eos, Cell *c) = 0;   // (3+1)D
	virtual void addCellBI(double tau, double dtau, Fluid *f, EoS *eos, Cell *c) = 0; // (2+1)D: boost-invariant
	virtual void addCellBISymm(double tau, double dtau, Fluid *f, EoS *eos, Cell *c) = 0; // (1+1)D: boost-invariant + azimuthally symmetric
	DileptonSpectrum GetSpectrum() const;
};

class DileptonsQGP : public Dileptons
{
	double coef;
	double coefBI;
	//double a, b, cc;
	double Fq;
	//double alphas;
	double tau0, taus;                // Undersaturated QGP
	std::vector<double> xphi, wphi;   // Legendre quadrature to integrate over phi
	std::vector<double> xeta, weta;   // Laguerre quadrature to integrate over eta (for boost-invariant case)
public:
	DileptonsQGP(double tau0_=0.1, double taus_=0., double Tcut_ = 0.155) : Dileptons(Tcut_), tau0(tau0_), taus(taus_) { init(); }
	DileptonsQGP(const std::vector<double> & Min, const std::vector<double> & yin, double tau0_=0.1, double taus_=0., double Tcut_ = 0.155);
	void init();
	virtual ~DileptonsQGP() { }
	virtual void addCell(double tau, double dtau, Fluid *f, EoS *eos, Cell *c);
	virtual void addCellBI(double tau, double dtau, Fluid *f, EoS *eos, Cell *c);
	virtual void addCellBISymm(double tau, double dtau, Fluid *f, EoS *eos, Cell *c);
};
