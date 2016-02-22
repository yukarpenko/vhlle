#include <fstream>
#include "icGlauberMC.h"
#include "fld.h"
#include "eos.h"
#include <iostream>
#include <vector>
#include <ctime>

using namespace std;

ICGlauberMC::ICGlauberMC(double e, double impactPar, double _tau0, int A_, double R_, double dlt_, double sigma_) :
	epsilon(e), b(impactPar), tau0(_tau0),
	A(A_), Ra(R_), dlt(dlt_), sigma(sigma_)
{
	cR = -1.;
	ca = -1.;
	rmax = fmax = -1.;
	ran = TRandom3(time(0));
}


ICGlauberMC::~ICGlauberMC(void)
{
}

void ICGlauberMC::WS_initparameters(double R, double a)
{
	cR = R;
	ca = a;
	
	double left = 0., right = R, center = 0.;
	while ( WS_root_function(right, R, a) > 0 ) right *= 2;

	while ( (right - left) / right > 1e-9)
	{
		center = (left + right) / 2.;
		if ( WS_root_function(center, R, a) > 0 ) left = center;
		else right = center;
	}

	rmax = (left + right) / 2.;
	fmax = WS_r(rmax, R, a);
}

double ICGlauberMC::Random_WS_r(double R, double a)
{
	if (abs(R-cR)>1e-9 || abs(a-ca)>1e-9) WS_initparameters(R, a);	// Check if parameters are not initialized
	
	double x = 0., y = 0.;

	while (1)
	{
		//x = ran.rand( 3 * R );
		//y = ran.rand( fmax );
		x = ran.Uniform(3 * R);
		y = ran.Uniform(fmax);
		if ( WS_r(x, R, a) > y ) return x;
	}
	return x;
}

void ICGlauberMC::setIC(Fluid *f, EoS *eos) {
  double e, nb, nq, vx = 0., vy = 0., vz = 0.;
  Cell *c;
  
  double rhad = 1.0;
  vector<double> xA(A), yA(A);
  vector<double> xB(A), yB(A);

  double xav = 0., yav = 0.;

  for(int i=0;i<A;++i) {
	double tr   = Random_WS_r(Ra, dlt);
	//double tcth = 2. * ran.rand() - 1.;
	//double tphi = 2. * ran.rand() * C_PI;
	double tcth = 2. * ran.Uniform() - 1.;
	double tphi = 2. * ran.Uniform() * C_PI;
	xA[i] = tr * sqrt(1. - tcth*tcth) * cos(tphi) + b / 2.;
	yA[i] = tr * sqrt(1. - tcth*tcth) * sin(tphi);

	xav += xA[i]; yav += yA[i];

	tr   = Random_WS_r(Ra, dlt);
	//tcth = 2. * ran.rand() - 1.;
	//tphi = 2. * ran.rand() * C_PI;
	tcth = 2. * ran.Uniform() - 1.;
	tphi = 2. * ran.Uniform() * C_PI;
	xB[i] = tr * sqrt(1. - tcth*tcth) * cos(tphi) - b / 2.;
	yB[i] = tr * sqrt(1. - tcth*tcth) * sin(tphi);
  }


  cout << "Average Xa and Ya: " <<  xav / A << " " << yav / A << "\n";

  vector<int> NWA(A, 0), NWB(A, 0);
  for(int ia=0;ia<A;++ia)
	  for(int ib=0;ib<A;++ib) {
		  if (((xA[ia]-xB[ib])*(xA[ia]-xB[ib]) + (yA[ia]-yB[ib])*(yA[ia]-yB[ib])) < sigma / C_PI) {
			NWA[ia]++;
			NWB[ib]++;
		  }
	  }

  double avv_num = 0., avv_den = 0.;
  double Etotal = 0.0, Etot2 = 0.0;

  for (int ix = 0; ix < f->getNX(); ix++)
    for (int iy = 0; iy < f->getNY(); iy++)
      for (int iz = 0; iz < f->getNZ(); iz++) {
        c = f->getCell(ix, iy, iz);
        double x = f->getX(ix);
        double y = f->getY(iy);
        double eta = f->getZ(iz);
        double etaFactor;
        //double eta1 = fabs(eta) < 1.3 ? 0.0 : fabs(eta) - 1.3;
        //etaFactor = exp(-eta1 * eta1 / 2.1 / 2.1) * (fabs(eta) < 5.3 ? 1.0 : 0.0);

		e = 0.;
		for(int ia=0;ia<A;++ia)
			if (NWA[ia]>0) e += Gauss(x, y, xA[ia], yA[ia], rhad);
		for(int ib=0;ib<A;++ib)
			if (NWB[ib]>0) e += Gauss(x, y, xB[ib], yB[ib], rhad);

		if (fabs(eta)<3.0) etaFactor = 1.;
		else etaFactor = exp(-(fabs(eta)-3.0)*(fabs(eta)-3.0)/2./0.4/0.4);

        e *= etaFactor * epsilon;
        //if (e < 0.5) e = 0.0;
        vx = vy = 0.0;
		nb = nq = 0.0;
        //nb = nq = eProfile(x, y) / 0.5;
        vz = 0.0;

      avv_num += sqrt(vx * vx + vy * vy) * e;
      avv_den += e;

        c->setPrimVar(eos, tau0, e, nb, nq, 0., vx, vy, vz);
        double _p = eos->p(e, nb, nq, 0., c->getTauP());
        const double gamma2 = 1.0 / (1.0 - vx * vx - vy * vy - vz * vz);
        Etotal +=
            ((e + _p) * gamma2 * (cosh(eta) + vz * sinh(eta)) - _p * cosh(eta));
        c->saveQprev();

        if (e > 0.) c->setAllM(1.);
      }

  cout << "average initial flow = " << avv_num / avv_den << endl;
  cout << "total energy = " << Etotal *f->getDx() * f->getDy() * f->getDz() *
                                   tau0 << endl;

  for(int ia=0;ia<A;++ia)
		if (NWA[ia]>0)  Etot2 += 1.;
  for(int ib=0;ib<A;++ib)
		if (NWB[ib]>0)  Etot2 += 1.;

  Etot2 *= epsilon;

  cout << "total energy2 = " << Etot2 * f->getDz() * 
                                   tau0 << endl;
}
