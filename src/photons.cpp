#include "photons.h"
#include "inc.h"
#include "NumericalIntegration.h"
#include "fld.h"
#include "eos.h"
#include "cll.h"

Photons::Photons(double Tcut_) : Tcut(Tcut_)
{
}

Photons::Photons(const std::vector<double> & ptin, const std::vector<double> & yin, double Tcut_) : pts(ptin), ys(yin), Tcut(Tcut_)
{
	yield.resize(ptin.size());
	v1.resize(ptin.size());
	v2.resize(ptin.size());
	v3.resize(ptin.size());
}


Photons::~Photons()
{
}

void Photons::addPtY(double pt, double y) { 
	pts.push_back(pt); 
	ys.push_back(y); 
	yield.push_back(0.); 
	v1.push_back(0.); 
	v2.push_back(0.); 
	v3.push_back(0.); 
}

PhotonSpectrum Photons::GetSpectrum() const {
	PhotonSpectrum ret;
	for(int ic=0;ic<pts.size();++ic) {
		PhotonSpectrumEntry ent;
		ent.pt    = pts[ic];
		ent.y     = ys[ic];
		ent.yield = yield[ic];
		ent.v1    = v1[ic] / yield[ic];
		ent.v2    = v2[ic] / yield[ic];
		ent.v3    = v3[ic] / yield[ic];
		ret.Entries.push_back(ent);
	}
	return ret;
}

PhotonsQGP::PhotonsQGP(const std::vector<double> & ptin, const std::vector<double> & yin, double tau0_, double taus_, double Tcut_) : Photons(ptin, yin, Tcut_), tau0(tau0_), taus(taus_)
{
}

void PhotonsQGP::init() {
	coef = 7.297352e-3 * 0.3 / pow(C_PI, 3.) / 4. * (2. / 3.);
	coefBI = coef * 2.;
	//a = 0.197; b = 0.987; cc = 0.884;
	a = 0.232; b = 0.987; cc = 0.884;
	alphas = 0.3;

	GetCoefsIntegrateLegendre32(0., 2*C_PI, xphi, wphi);
	GetCoefsIntegrateLaguerre32(xeta, weta);
}

void PhotonsQGP::addCell(double tau, double dtau, Fluid *f, EoS *eos, Cell *c) {
	double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz, s, Q[7], tauP;
	tauP = c->getTauP();
	f->getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
	eos->eos(e, nb, nq, ns, t, mub, muq, mus, p, tauP);
	double lambda = 1.;
	if (taus>0.0) lambda = 1. - exp((tau0-tauP)/taus);
	//s = eos->s(e, nb, nq, ns, tauP);

	if (t<Tcut) return;

	// Sum over all pT and Y values
	for(int ic = 0; ic < pts.size(); ++ic) {
		double pt = pts[ic], Y = ys[ic];
		double sum = 0., sumv1 = 0., sumv2 = 0., sumv3 = 0.;
		double eta = c->getZ();
		double gamma = 1. / sqrt(1. - vx*vx - vy*vy - tanh(vz)*tanh(vz));
		double Etil1 = gamma * pt * (cosh(Y) - tanh(vz) * sinh(Y)); 
		// Integral over phi
		for(int iphi = 0; iphi < xphi.size(); ++iphi) {
			double tphi = xphi[iphi];
			double cosphi = cos(tphi), sinphi = sin(tphi);
			double sumphi = 0.;
			double Etil2 = -gamma * pt * (vx*cosphi + vy*sinphi);
			double etil = Etil1 + Etil2;
			sumphi += exp(-etil / t) * ( lambda*lambda*(log(a*etil/alphas/t) + b*etil/t) + lambda*log(cc*etil/alphas/t) );
			sum    += wphi[iphi]                * sumphi;
			sumv1  += wphi[iphi] * cos(tphi)    * sumphi;
			sumv2  += wphi[iphi] * cos(2.*tphi) * sumphi;
			sumv3  += wphi[iphi] * cos(3.*tphi) * sumphi;
		}

		double tcoef = coef * f->getDx() * f->getDy() * f->getDz() * tau * dtau * t * t * pow(gevtofm, 4);
		sum       *= tcoef;
		sumv1     *= tcoef;
		sumv2     *= tcoef;
		sumv3     *= tcoef;
		yield[ic] += sum;
		v1[ic]    += sumv1;
		v2[ic]    += sumv2;
		v3[ic]    += sumv3;
	}
}

void PhotonsQGP::addCellBI(double tau, double dtau, Fluid *f, EoS *eos, Cell *c) {
	double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz, s, Q[7], tauP;
	tauP = c->getTauP();
	f->getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
	eos->eos(e, nb, nq, ns, t, mub, muq, mus, p, tauP);
	double lambda = 1.;
	if (taus>0.0) lambda = 1. - exp((tau0-tauP)/taus);
	//s = eos->s(e, nb, nq, ns, tauP);

	if (t<Tcut) return;

	// Sum over all pT values
	for(int ic = 0; ic < pts.size(); ++ic) {
		double pt = pts[ic];
		double sum = 0., sumv1 = 0., sumv2 = 0., sumv3 = 0.;
		// Integral over phi
		for(int iphi = 0; iphi < xphi.size(); ++iphi) {
			double tphi = xphi[iphi];
			double cosphi = cos(tphi), sinphi = sin(tphi);
			double sumphi = 0.;
			// Integral over eta
			for(int ieta = 0; ieta < xeta.size(); ++ieta) {
				double etil = (pt * cosh(xeta[ieta]) - vx * pt * cosphi - vy * pt * sinphi) / sqrt(1. - vx*vx - vy*vy);
				sumphi += weta[ieta] * exp(-etil / t) * ( lambda*lambda*(log(a*etil/alphas/t) + b*etil/t) + lambda*log(cc*etil/alphas/t) );
			}
			sum   += wphi[iphi]                * sumphi;
			sumv1 += wphi[iphi] * cos(tphi)    * sumphi;
			sumv2 += wphi[iphi] * cos(2.*tphi) * sumphi;
			sumv3 += wphi[iphi] * cos(3.*tphi) * sumphi;
		}

		double tcoef = coefBI * f->getDx() * f->getDy() * tau * dtau * t * t * pow(gevtofm, 4);
		sum   *= tcoef;
		sumv1 *= tcoef;
		sumv2 *= tcoef;
		sumv3 *= tcoef;
		yield[ic] += sum;
		v1[ic]    += sumv1;
		v2[ic]    += sumv2;
		v3[ic]    += sumv3;
	}
}

void PhotonsQGP::addCellBISymm(double tau, double dtau, Fluid *f, EoS *eos, Cell *c) {
	double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz, s, Q[7], tauP;
	tauP = c->getTauP();
	f->getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
	eos->eos(e, nb, nq, ns, t, mub, muq, mus, p, tauP);
	double lambda = 1.;
	if (taus>0.0) lambda = 1. - exp((tau0-tauP)/taus);
	//s = eos->s(e, nb, nq, ns, tauP);

	if (t<Tcut) return;

	// Sum over all pT values
	for(int ic = 0; ic < pts.size(); ++ic) {
		double pt = pts[ic];
		double sum = 0., sumv1 = 0., sumv2 = 0., sumv3 = 0.;
		
		double sumphi = 0.;
		// Integral over eta
		for(int ieta = 0; ieta < xeta.size(); ++ieta) {
			double etil = (pt * cosh(xeta[ieta]) - vx * pt) / sqrt(1. - vx*vx - vy*vy);
			sumphi += weta[ieta] * exp(-etil / t) * ( lambda*lambda*(log(a*etil/alphas/t) + b*etil/t) + lambda*log(cc*etil/alphas/t) );
		}
		sum   += 2. * C_PI * sumphi;
		sumv1 += 0.;
		sumv2 += 0.;
		sumv3 += 0.;

		double tcoef = coefBI * f->getDx() * f->getDy() * tau * dtau * t * t * pow(gevtofm, 4);
		//double tcoef = coefBI * tau * dtau * t * t * pow(gevtofm, 4);
		sum   *= tcoef;
		sumv1 *= tcoef;
		sumv2 *= tcoef;
		sumv3 *= tcoef;
		yield[ic] += sum;
		v1[ic]    += sumv1;
		v2[ic]    += sumv2;
		v3[ic]    += sumv3;
	}
}

PhotonsAMY::PhotonsAMY(const std::vector<double> & ptin, const std::vector<double> & yin, double tau0_, double taus_, double Tcut_) : Photons(ptin, yin, Tcut_), tau0(tau0_), taus(taus_)
{
}

void PhotonsAMY::init() {
	A1  = 1. / 6. / C_PI / C_PI;
	A2  = 1. / 3. / C_PI / C_PI;
	B1  = 1.000;
	B2  = 0.112;
	Nf  = 3.;
	Tc  = 0.170;
	Fq  = 2. / 3.;
	aEM = 7.297352e-3;
	
	//coef = 7.297352e-3 * 0.3 / pow(C_PI, 3.) / 4. * (2. / 3.);
	//coefBI = coef * 2.;
	//a = 0.197; b = 0.987; cc = 0.884;
	//a = 0.232; b = 0.987; cc = 0.884;
	//alphas = 0.3;

	GetCoefsIntegrateLegendre32(0., 2.*C_PI, xphi, wphi);
	GetCoefsIntegrateLaguerre32(xeta, weta);
}

void PhotonsAMY::addCell(double tau, double dtau, Fluid *f, EoS *eos, Cell *c) {
	double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz, s, Q[7], tauP;
	tauP = c->getTauP();
	f->getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
	eos->eos(e, nb, nq, ns, t, mub, muq, mus, p, tauP);
	double lambda = 1.;
	if (taus>0.0) lambda = 1. - exp((tau0-tauP)/taus);
	//s = eos->s(e, nb, nq, ns, tauP);

	if (t<Tcut) return;

	// Sum over all pT and Y values
	for(int ic = 0; ic < pts.size(); ++ic) {
		double pt = pts[ic], Y = ys[ic];
		double sum = 0., sumv1 = 0., sumv2 = 0., sumv3 = 0.;
		double eta = c->getZ();
		double gamma = 1. / sqrt(1. - vx*vx - vy*vy - tanh(vz)*tanh(vz));
		double Etil1 = gamma * pt * (cosh(Y) - tanh(vz) * sinh(Y)); 
		// Integral over phi
		for(int iphi = 0; iphi < xphi.size(); ++iphi) {
			double tphi = xphi[iphi];
			double cosphi = cos(tphi), sinphi = sin(tphi);
			double sumphi = 0.;
			double Etil2 = -gamma * pt * (vx*cosphi + vy*sinphi);
			double etil = Etil1 + Etil2;
			//sumphi += exp(-etil / t) * ( lambda*lambda*(log(a*etil/alphas/t) + b*etil/t) + lambda*log(cc*etil/alphas/t) );
			sumphi += GL(etil/t, t, lambda, fMode);
			sum    += wphi[iphi]                * sumphi;
			sumv1  += wphi[iphi] * cos(tphi)    * sumphi;
			sumv2  += wphi[iphi] * cos(2.*tphi) * sumphi;
			sumv3  += wphi[iphi] * cos(3.*tphi) * sumphi;
		}

		//double tcoef = coef * f->getDx() * f->getDy() * f->getDz() * tau * dtau * t * t * pow(gevtofm, 4);
		double tcoef = f->getDx() * f->getDy() * f->getDz() * tau * dtau * pow(gevtofm, 4) / 2. / C_PI;
		sum       *= tcoef;
		sumv1     *= tcoef;
		sumv2     *= tcoef;
		sumv3     *= tcoef;
		yield[ic] += sum;
		v1[ic]    += sumv1;
		v2[ic]    += sumv2;
		v3[ic]    += sumv3;
	}
}

void PhotonsAMY::addCellBI(double tau, double dtau, Fluid *f, EoS *eos, Cell *c) {
	double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz, s, Q[7], tauP;
	tauP = c->getTauP();
	f->getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
	eos->eos(e, nb, nq, ns, t, mub, muq, mus, p, tauP);
	double lambda = 1.;
	if (taus>0.0) lambda = 1. - exp((tau0-tauP)/taus);
	//s = eos->s(e, nb, nq, ns, tauP);

	if (t<Tcut) return;

	// Sum over all pT values
	for(int ic = 0; ic < pts.size(); ++ic) {
		double pt = pts[ic];
		double sum = 0., sumv1 = 0., sumv2 = 0., sumv3 = 0.;
		// Integral over phi
		for(int iphi = 0; iphi < xphi.size(); ++iphi) {
			double tphi = xphi[iphi];
			double cosphi = cos(tphi), sinphi = sin(tphi);
			double sumphi = 0.;
			// Integral over eta
			for(int ieta = 0; ieta < xeta.size(); ++ieta) {
				double etil = (pt * cosh(xeta[ieta]) - vx * pt * cosphi - vy * pt * sinphi) / sqrt(1. - vx*vx - vy*vy);
				//sumphi += weta[ieta] * exp(-etil / t) * ( lambda*lambda*(log(a*etil/alphas/t) + b*etil/t) + lambda*log(cc*etil/alphas/t) );
				sumphi += weta[ieta] * GL(etil/t, t, lambda, fMode);
			}
			sum   += wphi[iphi]                * sumphi;
			sumv1 += wphi[iphi] * cos(tphi)    * sumphi;
			sumv2 += wphi[iphi] * cos(2.*tphi) * sumphi;
			sumv3 += wphi[iphi] * cos(3.*tphi) * sumphi;
		}

		//double tcoef = coefBI * f->getDx() * f->getDy() * tau * dtau * t * t * pow(gevtofm, 4);
		double tcoef = 2. * f->getDx() * f->getDy() * tau * dtau * pow(gevtofm, 4) / 2. / C_PI;
		sum   *= tcoef;
		sumv1 *= tcoef;
		sumv2 *= tcoef;
		sumv3 *= tcoef;
		yield[ic] += sum;
		v1[ic]    += sumv1;
		v2[ic]    += sumv2;
		v3[ic]    += sumv3;
	}
}

void PhotonsAMY::addCellBISymm(double tau, double dtau, Fluid *f, EoS *eos, Cell *c) {
	double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz, s, Q[7], tauP;
	tauP = c->getTauP();
	f->getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz);
	eos->eos(e, nb, nq, ns, t, mub, muq, mus, p, tauP);
	double lambda = 1.;
	if (taus>0.0) lambda = 1. - exp((tau0-tauP)/taus);
	//s = eos->s(e, nb, nq, ns, tauP);

	if (t<Tcut) return;

	// Sum over all pT values
	for(int ic = 0; ic < pts.size(); ++ic) {
		double pt = pts[ic];
		double sum = 0., sumv1 = 0., sumv2 = 0., sumv3 = 0.;
		
		double sumphi = 0.;
		// Integral over eta
		for(int ieta = 0; ieta < xeta.size(); ++ieta) {
			double etil = (pt * cosh(xeta[ieta]) - vx * pt) / sqrt(1. - vx*vx - vy*vy);
			//sumphi += weta[ieta] * exp(-etil / t) * ( lambda*lambda*(log(a*etil/alphas/t) + b*etil/t) + lambda*log(cc*etil/alphas/t) );
			sumphi += weta[ieta] * GL(etil/t, t, lambda, fMode);
		}
		sum   += 2. * C_PI * sumphi;
		sumv1 += 0.;
		sumv2 += 0.;
		sumv3 += 0.;

		//double tcoef = coefBI * f->getDx() * f->getDy() * tau * dtau * t * t * pow(gevtofm, 4);
		//double tcoef = coefBI * tau * dtau * t * t * pow(gevtofm, 4);
		double tcoef = 2. * f->getDx() * f->getDy() * tau * dtau * pow(gevtofm, 4) / 2. / C_PI;
		sum   *= tcoef;
		sumv1 *= tcoef;
		sumv2 *= tcoef;
		sumv3 *= tcoef;
		yield[ic] += sum;
		v1[ic]    += sumv1;
		v2[ic]    += sumv2;
		v3[ic]    += sumv3;
	}
}

double PhotonsAMY::alphas(double T) const {
	return 6. * C_PI / (33. - 2. * Nf) / log(8. * T / Tc);
}

double PhotonsAMY::G1(double E, double T) const {
	return A1 * Fq * aEM * alphas(T) * T * T * exp(-E) * log(B1*E/alphas(T));
}

double PhotonsAMY::G2(double E, double T) const {
	return A2 * Fq * aEM * alphas(T) * T * T * exp(-E) * log(B2*E/alphas(T));
}

double PhotonsAMY::G(double E, double T) const {
	return 1. / C_PI / C_PI * Fq * aEM * alphas(T) * T * T / (exp(E) + 1) * (0.5 * log(3.*E/2./C_PI/alphas(T)) + C12(E) + C34(E));
}

double PhotonsAMY::C12(double E) const {
	return (0.041 / E) - 0.3615 + 1.01 * exp(-1.35 * E);
}

double PhotonsAMY::C34(double E) const {
	return sqrt(1. + Nf/6.) * (0.548 / pow(E, 3./2.) * log(12.28 + 1./E) + 0.133 * E / sqrt(1. + E/16.27));
}

double PhotonsAMY::GL(double E, double T,  double lambda, int mode) const {
	//return 4. / 3. / pow(C_PI, 4.) * aEM * 0.3 * T * T * exp(-E) * (2. * log(0.417*E/0.3) + 0.987 * E);
	if (mode==0) return lambda * G1(E,T) + lambda * lambda * (G(E,T) - G1(E,T));
	else return lambda * lambda * G2(E,T) + lambda * (G(E,T) - G2(E,T));
}