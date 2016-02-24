#include <vector>

class Fluid;
class EoS;
class Cell;

struct PhotonSpectrumEntry {
	double pt, y;
	double yield, v1, v2, v3;
};

struct PhotonSpectrum {
	std::vector<PhotonSpectrumEntry> Entries;
};

// Base class for calculating the yield of thermal photons
class Photons
{
protected:	
	std::vector<double> pts;     // Set of pt values to calculate spectrum at
	std::vector<double> ys;      // Set of y  values to calculate spectrum at
	std::vector<double> yield;   // Photon yields at corresponding y and pt
	std::vector<double> v1;      // Photon directed   flow at corresponding y and pt (unnormalized)
	std::vector<double> v2;      // Photon elliptic   flow at corresponding y and pt (unnormalized)
	std::vector<double> v3;      // Photon triangular flow at corresponding y and pt (unnormalized)
	double Tcut;				 // Cutoff temperature
public:
	Photons(double Tcut_ = 0.155);
	Photons(const std::vector<double> & ptin, const std::vector<double> & yin, double Tcut_ = 0.155);
	virtual ~Photons();
	void addPtY(double pt, double y);
	virtual void addCell(double tau, double dtau, Fluid *f, EoS *eos, Cell *c) = 0;   // (3+1)D
	virtual void addCellBI(double tau, double dtau, Fluid *f, EoS *eos, Cell *c) = 0; // (2+1)D: boost-invariant
	virtual void addCellBISymm(double tau, double dtau, Fluid *f, EoS *eos, Cell *c) = 0; // (1+1)D: boost-invariant + azimuthally symmetric
	PhotonSpectrum GetSpectrum() const;
	virtual void setMode(int mode_) = 0;
};

class PhotonsQGP : public Photons
{
	double coef;
	double coefBI;
	double a, b, cc;
	double alphas;
	double tau0, taus;                // Undersaturated QGP
	std::vector<double> xphi, wphi;   // Legendre quadrature to integrate over phi
	std::vector<double> xeta, weta;   // Laguerre quadrature to integrate over eta (for boost-invariant case)
public:
	PhotonsQGP(double tau0_=0.1, double taus_=0., double Tcut_ = 0.155) : Photons(Tcut_), tau0(tau0_), taus(taus_) { init(); }
	PhotonsQGP(const std::vector<double> & ptin, const std::vector<double> & yin, double tau0_=0.1, double taus_=0., double Tcut_ = 0.155);
	void init();
	virtual ~PhotonsQGP() { }
	virtual void addCell(double tau, double dtau, Fluid *f, EoS *eos, Cell *c);
	virtual void addCellBI(double tau, double dtau, Fluid *f, EoS *eos, Cell *c);
	virtual void addCellBISymm(double tau, double dtau, Fluid *f, EoS *eos, Cell *c);
};

class PhotonsAMY : public Photons
{
	double coef;
	double coefBI;
	double A1, A2;
	double B1, B2;
	double Fq;
	double Nf;
	double Tc;
	//double alphas;
	double aEM;
	double tau0, taus;                // Undersaturated QGP
	std::vector<double> xphi, wphi;   // Legendre quadrature to integrate over phi
	std::vector<double> xeta, weta;   // Laguerre quadrature to integrate over eta (for boost-invariant case)
	int fMode;
public:
	PhotonsAMY(double tau0_=0.1, double taus_=0., double Tcut_ = 0.155) : Photons(Tcut_), tau0(tau0_), taus(taus_), fMode(0) { init(); }
	PhotonsAMY(const std::vector<double> & ptin, const std::vector<double> & yin, double tau0_=0.1, double taus_=0., double Tcut_ = 0.155);
	void init();
	void setMode(int mode_ = 0) { fMode = mode_; }
	virtual ~PhotonsAMY() { }
	virtual void addCell(double tau, double dtau, Fluid *f, EoS *eos, Cell *c);
	virtual void addCellBI(double tau, double dtau, Fluid *f, EoS *eos, Cell *c);
	virtual void addCellBISymm(double tau, double dtau, Fluid *f, EoS *eos, Cell *c);
	double alphas(double T) const;
	double C12(double E) const;
	double C34(double E) const;
	double G1(double E, double T) const;
	double G2(double E, double T) const;
	double G(double E, double T) const;
	double GL(double E, double T, double la, int mode = 0) const;
};
