#include "dileptons.h"
#include "inc.h"
#include "NumericalIntegration.h"
#include "fld.h"
#include "eos.h"
#include "cll.h"
#include "TMath.h"

Dileptons::Dileptons(double Tcut_) : Tcut(Tcut_)
{
}

Dileptons::Dileptons(const std::vector<double> & Min, const std::vector<double> & yin, double Tcut_) : Ms(Min), ys(yin), Tcut(Tcut_)
{
	yield.resize(Min.size());
	//v1.resize(ptin.size());
	v2.resize(Min.size());
	//v3.resize(ptin.size());
}


Dileptons::~Dileptons()
{
}

void Dileptons::addPtY(double M, double y) { 
	Ms.push_back(M); 
	ys.push_back(y); 
	yield.push_back(0.); 
	//v1.push_back(0.); 
	v2.push_back(0.); 
	//v3.push_back(0.); 
}

DileptonSpectrum Dileptons::GetSpectrum() const {
	DileptonSpectrum ret;
	for(int ic=0;ic<Ms.size();++ic) {
		DileptonSpectrumEntry ent;
		ent.M     = Ms[ic];
		ent.y     = ys[ic];
		ent.yield = yield[ic];
		//ent.v1    = v1[ic] / yield[ic];
		ent.v2    = v2[ic] / yield[ic];
		//ent.v3    = v3[ic] / yield[ic];
		ret.Entries.push_back(ent);
	}
	return ret;
}

DileptonsQGP::DileptonsQGP(const std::vector<double> & Min, const std::vector<double> & yin, double tau0_, double taus_, double Tcut_) : Dileptons(Min, yin, Tcut_), tau0(tau0_), taus(taus_)
{
}

void DileptonsQGP::init() {
	Fq = 2./3.;
	coef = 7.297352e-3 * 7.297352e-3 * Fq / 4. / pow(C_PI, 3.);//7.297352e-3 * 0.3 / pow(C_PI, 5.) * (2. / 3.);
	coefBI = coef * 2.;
	//a = 0.197; b = 0.987; cc = 0.884;
	//alphas = 0.3;
	

	GetCoefsIntegrateLegendre32(0., 2*C_PI, xphi, wphi);
	GetCoefsIntegrateLaguerre32(xeta, weta);
}

// TODO v2 for (3+1)D
void DileptonsQGP::addCell(double tau, double dtau, Fluid *f, EoS *eos, Cell *c) {
	double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz, s, Q[7], tauP;
	tauP = c->getTauP();
	f->getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
	eos->eos(e, nb, nq, ns, t, mub, muq, mus, p, tauP);
	double lambda = 1.;
	if (taus>0.0) lambda = 1. - exp((tau0-tauP)/taus);
	//s = eos->s(e, nb, nq, ns, tauP);

	if (t<Tcut) return;

	// Sum over all M and Y values
	for(int ic = 0; ic < Ms.size(); ++ic) {
		double M = Ms[ic], Y = ys[ic];
		double Mtil = M / t;
		double sum = 0.;//, sumv1 = 0., sumv2 = 0., sumv3 = 0.;
		double sumv2 = 0.;
		double eta = c->getZ();
		double gamma = 1. / (1. - vx*vx - vy*vy - tanh(vz)*tanh(vz));
		double vT = sqrt(vx*vx + vy*vy);

		for(int iqt = 0; iqt < xeta.size(); ++iqt) {
			double tqt  = xeta[iqt];
			double tsum = lambda * lambda * t * t * tqt * TMath::BesselI0(gamma * vT * tqt) * exp(-gamma * sqrt(Mtil*Mtil+tqt*tqt) * (cosh(Y) - tanh(vz) * sinh(Y)));
			sum += weta[iqt] * tsum;
		}

		double tcoef = coef * f->getDx() * f->getDy() * f->getDz() * tau * dtau * pow(gevtofm, 4);
		sum  *= tcoef;
		yield[ic] += sum;
	}
}

void DileptonsQGP::addCellBI(double tau, double dtau, Fluid *f, EoS *eos, Cell *c) {
	double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz, s, Q[7], tauP;
	tauP = c->getTauP();
	f->getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
	eos->eos(e, nb, nq, ns, t, mub, muq, mus, p, tauP);
	double lambda = 1.;
	if (taus>0.0) lambda = 1. - exp((tau0-tauP)/taus);
	//s = eos->s(e, nb, nq, ns, tauP);

	if (t<Tcut) return;

	// Sum over all M and Y values
	for(int ic = 0; ic < Ms.size(); ++ic) {
		double M = Ms[ic], Y = ys[ic];
		double Mtil = M / t;
		double sum = 0.;//, sumv1 = 0., sumv2 = 0., sumv3 = 0.;
		double sumv2 = 0.;
		double eta = c->getZ();
		double gamma = 1. / sqrt(1. - vx*vx - vy*vy);
		double vT = sqrt(vx*vx + vy*vy);

		/*if (0) {
			for(int iqt = 0; iqt < xeta.size(); ++iqt) {
				double tqt  = xeta[iqt];
				double tsum = lambda * lambda * t * t * tqt * TMath::BesselI0(gamma * vT * tqt) * TMath::BesselK0(gamma * sqrt(Mtil*Mtil+tqt*tqt));
				sum += weta[iqt] * tsum;
			}
		}
		else*/ 
		{
			sum = lambda * lambda * t * M * TMath::BesselK1(M / t);
			if (vT>1.e-6) sumv2 = lambda * lambda * (vx*vx - vy*vy) / (vx*vx + vy*vy) * (t * M * TMath::BesselK1(M / t) - 2. * t * t / (gamma*gamma-1.) * (TMath::BesselK0(M / t)-TMath::BesselK0(gamma * M / t) ) ) ;
			else sumv2 = lambda * lambda * (vx*vx - vy*vy) / 4. * M * M * TMath::BesselK(2, M / t);
		}

		double tcoef = coefBI * f->getDx() * f->getDy() * tau * dtau * pow(gevtofm, 4);
		sum   *= tcoef;
		sumv2 *= tcoef;
		yield[ic] += sum;
		v2[ic]    += sumv2;
	}
}

void DileptonsQGP::addCellBISymm(double tau, double dtau, Fluid *f, EoS *eos, Cell *c) {
	addCellBI(tau, dtau, f, eos, c);
}
