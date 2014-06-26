#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "hdo.h"
#include "inc.h"
#include "rmn.h"
#include "fld.h"
#include "eos.h"
#include "cll.h"
#include "trancoeff.h"

using namespace std ;

extern bool debugRiemann ;

extern double sign(double x) ;

// this version contains does contain pre-advection for dt/2 for the IS solution,
// which cause some oscillations at the edges in diagonal direction at late times

// use formal solution for the relaxation part of Israel-Stewart equations (not recommended)
//#define FORMAL_SOLUTION
// else: use 1st order finite difference update

//ofstream foutNS ;


Hydro::Hydro(Fluid *_f, EoS *_eos, TransportCoeff *_trcoeff, double _t0, double _dt)
{
	eos = _eos ;
	trcoeff = _trcoeff ;
	f = _f ;
	dt = _dt ;
	tau = _t0 ;
	//foutNS.open("outputNS.dat") ;
}

Hydro::~Hydro()
{
}


void Hydro::setDtau(double deltaTau)
{
  dt = deltaTau ;
  if(dt>f->getDx()/2. || dt>f->getDy()/2. || dt>tau*f->getDz()/2.){
    cout << "too big delta_tau " << dt << "  " << f->getDx() << "  " << f->getDy() << "  " << tau*f->getDz() << endl ;
    exit(1) ;
  }
}


void Hydro::hlle_flux(Cell *left, Cell *right, int direction, int mode)
{
	double el, er, pl, pr, nbl, nql, nsl, nbr, nqr, nsr, vxl, vxr, vyl, vyr, vzl, vzr, bl, br, csb, vb, El, Er, dx ;
	double Ftl, Fxl, Fyl, Fzl, Fbl, Fql, Fsl, Ftr, Fxr, Fyr, Fzr, Fbr, Fqr, Fsr ;
	double U1l, U2l, U3l, U4l, Ubl, Uql, Usl, U1r, U2r, U3r, U4r, Ubr, Uqr, Usr ;
	double flux [7] ;
	const double dta = mode==0 ? dt/2. : dt ;
	double tauFactor ; // fluxes are also multiplied by tau
	if(mode==PREDICT){
	left->getPrimVarRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl, direction) ;
	right->getPrimVarLeft(eos, tau, er, pr, nbr, nqr, nsr, vxr, vyr, vzr, direction) ;
	El = (el+pl)/(1-vxl*vxl-vyl*vyl-vzl*vzl) ;
	Er = (er+pr)/(1-vxr*vxr-vyr*vyr-vzr*vzr) ;
	tauFactor = tau ;
	}else{
	left->getPrimVarHRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl, direction) ;
	right->getPrimVarHLeft(eos, tau, er, pr, nbr, nqr, nsr, vxr, vyr, vzr, direction) ;
	El = (el+pl)/(1-vxl*vxl-vyl*vyl-vzl*vzl) ;
	Er = (er+pr)/(1-vxr*vxr-vyr*vyr-vzr*vzr) ;
	tauFactor = tau + dt/2. ;
	}

	if(el<0.) { el=0. ; pl=0. ; }
	if(er<0.) { er=0. ; pr=0. ; }
	
	if(el>1e10){
		cout << "e>1e10; debug info below:\n" ;
		left->Dump(tau) ;
		debugRiemann = true ;
	if(mode==PREDICT)
	left->getPrimVarRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl, direction) ;
	else
	left->getPrimVarHRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl, direction) ;
	debugRiemann = false ;
	exit(0) ;
	}

	//vacuum treatment
	if(el==0. && er==0.) return ;
	if(pr<0.){
		cout << "Negative pressure" << endl ;
		left->getPrimVarRight(eos, tau, el, pl, nbl, nql, nsl, vxl, vyl, vzl, direction) ;
		right->getPrimVarLeft(eos, tau, er, pr, nbr, nqr, nsr, vxr, vyr, vzr, direction) ;
	}

	//
	if(left->getM(direction)<1. && right->getM(direction)<1.) return ;

	double gammal = 1./sqrt(1-vxl*vxl-vyl*vyl-vzl*vzl) ;
	double gammar = 1./sqrt(1-vxr*vxr-vyr*vyr-vzr*vzr) ;
	U1l = gammal*gammal*(el+pl)*vxl ;
	U2l = gammal*gammal*(el+pl)*vyl ;
	U3l = gammal*gammal*(el+pl)*vzl ;
	U4l = gammal*gammal*(el+pl)-pl ;
	Ubl = gammal*nbl ;
	Uql = gammal*nql ;
	Usl = gammal*nsl ;

	U1r = gammar*gammar*(er+pr)*vxr ;
	U2r = gammar*gammar*(er+pr)*vyr ;
	U3r = gammar*gammar*(er+pr)*vzr ;
	U4r = gammar*gammar*(er+pr)-pr ;
	Ubr = gammar*nbr ;
	Uqr = gammar*nqr ;
	Usr = gammar*nsr ;

	if(direction==X_){
		Ftl = U4l*vxl+pl*vxl ;
		Fxl = U1l*vxl+pl ;
		Fyl = U2l*vxl ;
		Fzl = U3l*vxl ;
		Fbl = Ubl*vxl ;
		Fql = Uql*vxl ;
		Fsl = Usl*vxl ;

		Ftr = U4r*vxr+pr*vxr ;
		Fxr = U1r*vxr+pr ;
		Fyr = U2r*vxr ;
		Fzr = U3r*vxr ;
		Fbr = Ubr*vxr ;
		Fqr = Uqr*vxr ;
		Fsr = Usr*vxr ;

			// for the case of constant c_s only
		csb = sqrt(eos->cs2() + 0.5*sqrt(El*Er)/pow(sqrt(El)+sqrt(Er),2)*pow(vxl-vxr,2)) ;
		vb = (sqrt(El)*vxl+sqrt(Er)*vxr)/(sqrt(El)+sqrt(Er)) ;
		bl = min(0., min((vb-csb)/(1-vb*csb), (vxl-eos->cs())/(1-vxl*eos->cs()) )) ;
		br = max(0., max((vb+csb)/(1+vb*csb), (vxr+eos->cs())/(1+vxr*eos->cs()) )) ;

		dx = f->getDx() ;

			// boundary with vacuum
		if(el==0.) bl = -1. ;
		if(er==0.) br = 1. ;
	}
	if(direction==Y_){
		Ftl = U4l*vyl+pl*vyl ;
		Fxl = U1l*vyl ;
		Fyl = U2l*vyl+pl ;
		Fzl = U3l*vyl ;
		Fbl = Ubl*vyl ;
		Fql = Uql*vyl ;
		Fsl = Usl*vyl ;

		Ftr = U4r*vyr+pr*vyr ;
		Fxr = U1r*vyr ;
		Fyr = U2r*vyr+pr ;
		Fzr = U3r*vyr ;
		Fbr = Ubr*vyr ;
		Fqr = Uqr*vyr ;
		Fsr = Usr*vyr ;

			// for the case of constant c_s only
		csb = sqrt(eos->cs2() + 0.5*sqrt(El*Er)/pow(sqrt(El)+sqrt(Er),2)*pow(vyl-vyr,2)) ;
		vb = (sqrt(El)*vyl+sqrt(Er)*vyr)/(sqrt(El)+sqrt(Er)) ;
		bl = min(0., min((vb-csb)/(1-vb*csb), (vyl-eos->cs())/(1-vyl*eos->cs()) )) ;
		br = max(0., max((vb+csb)/(1+vb*csb), (vyr+eos->cs())/(1+vyr*eos->cs()) )) ;

		dx = f->getDy() ;

			// boundary with vacuum
		if(el==0.) bl = -1. ;
		if(er==0.) br = 1. ;
	}
	if(direction==Z_){
		double tau1 = tau_z ;
		Ftl = U4l*vzl/tau1+pl*vzl/tau1 ;
		Fxl = U1l*vzl/tau1 ;
		Fyl = U2l*vzl/tau1 ;
		Fzl = U3l*vzl/tau1+pl/tau1 ;
		Fbl = Ubl*vzl/tau1 ;
		Fql = Uql*vzl/tau1 ;
		Fsl = Usl*vzl/tau1 ;

		Ftr = U4r*vzr/tau1+pr*vzr/tau1 ;
		Fxr = U1r*vzr/tau1 ;
		Fyr = U2r*vzr/tau1 ;
		Fzr = U3r*vzr/tau1+pr/tau1 ;
		Fbr = Ubr*vzr/tau1 ;
		Fqr = Uqr*vzr/tau1 ;
		Fsr = Usr*vzr/tau1 ;

		// for the case of constant c_s only
		// factor 1/tau accounts for eta-coordinate

		//different estimate
		csb = sqrt(eos->cs2() + 0.5*sqrt(El*Er)/pow(sqrt(El)+sqrt(Er),2)*pow(vzl-vzr,2)) ;
		vb = (sqrt(El)*vzl+sqrt(Er)*vzr)/(sqrt(El)+sqrt(Er)) ;
		bl = 1./tau*min(0., min((vb-csb)/(1-vb*csb), (vzl-eos->cs())/(1-vzl*eos->cs()) )) ;
		br = 1./tau*max(0., max((vb+csb)/(1+vb*csb), (vzr+eos->cs())/(1+vzr*eos->cs()) )) ;

		dx = f->getDz() ;

			// boundary with vacuum
		if(el==0.) bl = -1./tau ;
		if(er==0.) br = 1./tau ;
	}
//	bl = min(0., (vxl-EoS::cs(el))/(1-vxl*EoS::cs(el)) ) ;
//	br = max(0., (vxr+EoS::cs(er))/(1+vxr*EoS::cs(er)) ) ;
	
	if(bl==0. && br==0.) return ;


//	if(br>bl){
//	double Qr [7] ; right->getQ(Qr) ; // incorrect !
//	double Ql [7] ; left->getQ(Ql) ;
	flux[T_] = tauFactor*dta/dx*( -bl*br*(U4l-U4r) + br*Ftl - bl*Ftr ) / (-bl + br) ;
	flux[X_] = tauFactor*dta/dx*( -bl*br*(U1l-U1r) + br*Fxl - bl*Fxr ) / (-bl + br) ;
	flux[Y_] = tauFactor*dta/dx*( -bl*br*(U2l-U2r) + br*Fyl - bl*Fyr ) / (-bl + br) ;
	flux[Z_] = tauFactor*dta/dx*( -bl*br*(U3l-U3r) + br*Fzl - bl*Fzr ) / (-bl + br) ;
	flux[NB_] = tauFactor*dta/dx*( -bl*br*(Ubl-Ubr) + br*Fbl - bl*Fbr ) / (-bl + br) ;
	flux[NQ_] = tauFactor*dta/dx*( -bl*br*(Uql-Uqr) + br*Fql - bl*Fqr ) / (-bl + br) ;
	flux[NS_] = tauFactor*dta/dx*( -bl*br*(Usl-Usr) + br*Fsl - bl*Fsr ) / (-bl + br) ;
//	}
//	else{// br<bl
//		cout << "HLLE : failed\n" ; exit(1) ;
//		if(fabs(bl)>fabs(br)) { flux[T_] = dta/f->getDx()*Ftl ; flux[X_] = dta/f->getDx()*Fxl ; flux[Y_] = flux[Z_] = 0. ; }
//		else { flux[T_] = dta/f->getDx()*Ftr ; flux[X_] = dta/f->getDx()*Fxr ; flux[Y_] = flux[Z_] = 0. ; }
//	}

	if(!(flux[NB_]>=0. || flux[NB_]<0.)){
		cout<<"---- error in hlle_flux: f_nb undefined!\n" ;
		cout << setw(12) << U4l << setw(12) << U1l << setw(12) << U2l << setw(12) << U3l << endl ;
		cout << setw(12) << U4r << setw(12) << U1r << setw(12) << U2r << setw(12) << U3r << endl ;
		cout << setw(12) << Ubl << setw(12) << Uql << setw(12) << Usl << endl ;
		cout << setw(12) << Ubr << setw(12) << Uqr << setw(12) << Usr << endl ;
		cout << setw(12) << Ftl << setw(12) << Fxl << setw(12) << Fyl << setw(12) << Fzl << endl ;
		cout << setw(12) << Ftr << setw(12) << Fxr << setw(12) << Fyr << setw(12) << Fzr << endl ;
		exit(1) ;
	}
	
	left->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], -flux[NB_], -flux[NQ_], -flux[NS_]) ;
	right->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], flux[NB_], flux[NQ_], flux[NS_]) ;

}


void Hydro::source(double tau1, double x, double y, double z, double e, double p, double nb, double nq, double ns,
				   double vx, double vy, double vz, double S[7])
{
	double Q [7] ;
	transformCV(e, p, nb, nq, ns, vx, vy, vz, Q) ;
	S[T_] = -Q[T_]*vz*vz - p*(1.+vz*vz) ;
	S[X_] = 0. ;
	S[Y_] = 0. ;
	S[Z_] = -Q[Z_] ;
	S[NB_] = 0. ;
	S[NQ_] = 0. ;
	S[NS_] = 0. ;
}


void Hydro::source_step(int ix, int iy, int iz, int mode)
{
	double _dt ;
	if(mode==PREDICT) 
		_dt = dt/2. ;
	else 
		_dt = dt ;

	double e, p, nb, nq, ns, vx, vy, vz ;
	double k[7] ;

	double x = f->getX(ix), y = f->getY(iy), z = f->getZ(iz) ;
	Cell *c = f->getCell(ix, iy, iz) ;

	if(mode==PREDICT){
	c->getPrimVar(eos, tau, e, p, nb, nq, ns, vx, vy, vz) ;
//	c->getQ(Q) ;
	}else{
	c->getPrimVarHCenter(eos, tau+dt/2., e, p, nb, nq, ns, vx, vy, vz) ;
//	c->getQh(Q) ;
	}

	source(tau, x, y, z, e, p, nb, nq, ns, vx, vy, vz, k) ;  // setting k1
	for(int i=0; i<7; i++) k[i] *= _dt ;

	if(!(k[NB_]>=0. || k[NB_]<0.)){
		cout<<"---- error in source_step: k_nb undefined!\n" ;
		cout << setw(12) << k[0] << setw(12) << k[1] << setw(12) << k[2] << setw(12) << k[3] << endl ;
		cout << setw(12) << k[4] << setw(12) << k[5] << setw(12) << k[6] << endl ;
		exit(1) ;
	}
		c->addFlux(k[T_], k[X_], k[Y_], k[Z_], k[NB_], k[NQ_], k[NS_]) ;
}


void Hydro::visc_source_step(int ix, int iy, int iz)
{
	double e, p, nb, nq, ns, vx, vy, vz ;
	double uuu[4] ;
	double k[7] ;

	Cell *c = f->getCell(ix, iy, iz) ;

	c->getPrimVarHCenter(eos, tau-dt/2., e, p, nb, nq, ns, vx, vy, vz) ;
	if(e<=0.) return ;
	uuu[0]=1./sqrt(1.-vx*vx-vy*vy-vz*vz) ;
	uuu[1]=uuu[0]*vx ;
	uuu[2]=uuu[0]*vy ;
	uuu[3]=uuu[0]*vz ;

	k[T_] = -(c->getpiH(0,0)+c->getpiH(3,3)-c->getPiH()*(-uuu[0]*uuu[0]-uuu[3]*uuu[3]))/(tau-0.5*dt) ;
	k[X_] = -(c->getpiH(0,1)+c->getPiH()*uuu[0]*uuu[1])/(tau-0.5*dt) ;
	k[Y_] = -(c->getpiH(0,2)+c->getPiH()*uuu[0]*uuu[2])/(tau-0.5*dt) ;
	k[Z_] = -2.0*(c->getpiH(0,3)+c->getPiH()*uuu[0]*uuu[3])/(tau-0.5*dt) ;
	for(int i=0; i<4; i++) k[i] *= dt ;

//	if(!(k[NB_]>=0. || k[NB_]<0.)){
//		cout<<"---- error in source_step: k_nb undefined!\n" ;
//		cout << setw(12) << k[0] << setw(12) << k[1] << setw(12) << k[2] << setw(12) << k[3] << endl ;
//		cout << setw(12) << k[4] << setw(12) << k[5] << setw(12) << k[6] << endl ;
//		exit(1) ;
//	}
//----- debug
	//if(ix==36&&iy==50&&iz==2){
	 //cout<<"visc_source_step:" ;
	 //for(int i=0; i<4; i++)cout<<setw(14)<<k[i]; cout<<endl ;
	//}
		c->addFlux(k[T_], k[X_], k[Y_], k[Z_], 0., 0., 0.) ;
}

// dv/d(tau) = v^{t+dt}_ideal - v^{t}
// dv/dx_i ~ v^{x+dx}-v{x-dx}
// makes sense after non-viscous step
void Hydro::NSquant(int ix, int iy, int iz, double pi[4][4], double& Pi, double dmu[4][4], double &du)
{
	const double VMIN=1e-2 ;
	const double UDIFF=3.0 ;
	double e0, e1, p, nb, nq, ns, vx1, vy1, vz1, vx0, vy0, vz0, vxH, vyH, vzH ;
	double ut0, ux0, uy0, uz0, ut1, ux1, uy1, uz1 ;
//	double dmu [4][4] ; // \partial_\mu u^\nu matrix
	// coordinates: 0=tau, 1=x, 2=y, 3=eta
	double Z [4][4][4][4] ; // Z[mu][nu][lambda][rho]
	double uuu [4] ; // the 4-velocity
	double gmunu [4][4] = {{1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1}} ; // omit 1/tau^2 in g^{eta,eta}
	Cell *c = f->getCell(ix, iy, iz) ;
	double dx=f->getDx(), dy=f->getDy(), dz=f->getDz() ;
	// check if the cell is next to vacuum from +-x, +-y side:
	if(c->getNext(X_)->getMaxM()<=0.9 || c->getNext(Y_)->getMaxM()<=0.9 ||
	c->getPrev(X_)->getMaxM()<=0.9 || c->getPrev(Y_)->getMaxM()<=0.9 ||
	f->getCell(ix+1,iy+1,iz)->getMaxM()<=0.9 || f->getCell(ix+1,iy-1,iz)->getMaxM()<=0.9 ||
	f->getCell(ix-1,iy+1,iz)->getMaxM()<=0.9 || f->getCell(ix-1,iy-1,iz)->getMaxM()<=0.9){
	 for(int i=0; i<4; i++)
	 for(int j=0; j<4; j++){ pi[i][j]=0. ; dmu[i][j]=0. ; }
	 Pi=du=0. ;
	 return ;
	}
	// calculation of \partial_\mu u^\nu matrix
	// mu=first index, nu=second index
  // centered differences with respect to the values at (it+1/2, ix, iy, iz)
		// d_tau u^\mu
		c->getPrimVarPrev(eos, tau-dt, e0, p, nb, nq, ns, vx0, vy0, vz0) ;
		c->getPrimVar(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1) ;
    c->getPrimVarHCenter(eos, tau-0.5*dt, e1, p, nb, nq, ns, vxH, vyH, vzH) ;
	//############## get transport coefficients
		double T, mub, muq, mus ;
		double etaS, zetaS ;
    double s = eos->s(e1,nb,nq,ns) ; // entropy density in the current cell
		eos->eos(e1,nb,nq,ns,T,mub,muq,mus,p) ;
		trcoeff->getEta(e1, T, etaS, zetaS) ;
	//##############
		//if(e1<0.00004) s=0. ; // negative pressure due to pi^zz for small e
		ut0=1.0/sqrt(1.0-vx0*vx0-vy0*vy0-vz0*vz0) ;
		ux0=ut0*vx0 ; uy0=ut0*vy0 ; uz0=ut0*vz0 ;
		ut1=1.0/sqrt(1.0-vx1*vx1-vy1*vy1-vz1*vz1) ;
		ux1=ut1*vx1 ; uy1=ut1*vy1 ; uz1=ut1*vz1 ;
    uuu[0]=1.0/sqrt(1.0-vxH*vxH-vyH*vyH-vzH*vzH);
    uuu[1]=uuu[0]*vxH; uuu[2]=uuu[0]*vyH; uuu[3]=uuu[0]*vzH ;
    
		dmu[0][0]=(ut1*ut1-ut0*ut0)/2./uuu[0]/dt ; dmu[0][1]=(ux1*ux1-ux0*ux0)/2./uuu[1]/dt ;
		dmu[0][2]=(uy1*uy1-uy0*uy0)/2./uuu[2]/dt ; dmu[0][3]=(uz1*uz1-uz0*uz0)/2./uuu[3]/dt ;
		if(fabs(0.5*(ut1+ut0)/ut1)>UDIFF) dmu[0][0]=(ut1-ut0)/dt ;
		if(fabs(ux1)<VMIN || fabs(0.5*(ux1+ux0)/ux1)>UDIFF) dmu[0][1]=(ux1-ux0)/dt ;
		if(fabs(uy1)<VMIN || fabs(0.5*(uy1+uy0)/uy1)>UDIFF) dmu[0][2]=(uy1-uy0)/dt ;
		if(fabs(uz1)<VMIN || fabs(0.5*(uz1+uz0)/uz1)>UDIFF) dmu[0][3]=(uz1-uz0)/dt ;
		if(e1<=0. || e0<=0.){ // matter-vacuum
		dmu[0][0]=dmu[0][1]=dmu[0][2]=dmu[0][3]=0. ;
		}
		//if(fabs(dmu[0][0])>200){
			//cout<<"|dmu[0][0]|>2000  triggered"<<ut1<<"  "<<ut0<<endl ;
			//exit(1) ;
		//}
		// d_x u^\mu
		f->getCell(ix+1, iy, iz)->getPrimVarHCenter(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1) ;
		f->getCell(ix-1, iy, iz)->getPrimVarHCenter(eos, tau, e0, p, nb, nq, ns, vx0, vy0, vz0) ;
		if(e1>0. && e0>0.){
		ut0=1.0/sqrt(1.0-vx0*vx0-vy0*vy0-vz0*vz0) ;
		ux0=ut0*vx0 ; uy0=ut0*vy0 ; uz0=ut0*vz0 ;
		ut1=1.0/sqrt(1.0-vx1*vx1-vy1*vy1-vz1*vz1) ;
		ux1=ut1*vx1 ; uy1=ut1*vy1 ; uz1=ut1*vz1 ;
		dmu[1][0]=0.25*(ut1*ut1-ut0*ut0)/uuu[0]/dx ; dmu[1][1]=0.25*(ux1*ux1-ux0*ux0)/uuu[1]/dx ;
		dmu[1][2]=0.25*(uy1*uy1-uy0*uy0)/uuu[2]/dx ; dmu[1][3]=0.25*(uz1*uz1-uz0*uz0)/uuu[3]/dx ;
		if(fabs(0.5*(ut1+ut0)/uuu[0])>UDIFF) dmu[1][0]=0.5*(ut1-ut0)/dx ;
		if(fabs(uuu[1])<VMIN || fabs(0.5*(ux1+ux0)/uuu[1])>UDIFF) dmu[1][1]=0.5*(ux1-ux0)/dx ;
		if(fabs(uuu[2])<VMIN || fabs(0.5*(uy1+uy0)/uuu[2])>UDIFF) dmu[1][2]=0.5*(uy1-uy0)/dx ;
		if(fabs(uuu[3])<VMIN || fabs(0.5*(uz1+uz0)/uuu[3])>UDIFF) dmu[1][3]=0.5*(uz1-uz0)/dx ;
		}else{ // matter-vacuum
		dmu[1][0]=dmu[1][1]=dmu[1][2]=dmu[1][3]=0. ;
		}
		if(fabs(dmu[1][3])>1e+10)cout<<"dmu[1][3]:  "<<uz1<<"  "<<uz0<<"  "<<uuu[3]<<endl ;
		// d_y u^\mu
		f->getCell(ix, iy+1, iz)->getPrimVarHCenter(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1) ;
		f->getCell(ix, iy-1, iz)->getPrimVarHCenter(eos, tau, e0, p, nb, nq, ns, vx0, vy0, vz0) ;
		if(e1>0. && e0>0.){
		ut0=1.0/sqrt(1.0-vx0*vx0-vy0*vy0-vz0*vz0) ;
		ux0=ut0*vx0 ; uy0=ut0*vy0 ; uz0=ut0*vz0 ;
		ut1=1.0/sqrt(1.0-vx1*vx1-vy1*vy1-vz1*vz1) ;
		ux1=ut1*vx1 ; uy1=ut1*vy1 ; uz1=ut1*vz1 ;
		dmu[2][0]=0.25*(ut1*ut1-ut0*ut0)/uuu[0]/dy ; dmu[2][1]=0.25*(ux1*ux1-ux0*ux0)/uuu[1]/dy ;
		dmu[2][2]=0.25*(uy1*uy1-uy0*uy0)/uuu[2]/dy ; dmu[2][3]=0.25*(uz1*uz1-uz0*uz0)/uuu[3]/dy ;
		if(fabs(0.5*(ut1+ut0)/uuu[0])>UDIFF) dmu[2][0]=0.5*(ut1-ut0)/dy ;
		if(fabs(uuu[1])<VMIN || fabs(0.5*(ux1+ux0)/uuu[1])>UDIFF) dmu[2][1]=0.5*(ux1-ux0)/dy ;
		if(fabs(uuu[2])<VMIN || fabs(0.5*(uy1+uy0)/uuu[2])>UDIFF) dmu[2][2]=0.5*(uy1-uy0)/dy ;
		if(fabs(uuu[3])<VMIN || fabs(0.5*(uz1+uz0)/uuu[3])>UDIFF) dmu[2][3]=0.5*(uz1-uz0)/dy ;
		}else{ // matter-vacuum
		dmu[2][0]=dmu[2][1]=dmu[2][2]=dmu[2][3]=0. ;
		}
		// d_z u^\mu
		f->getCell(ix, iy, iz+1)->getPrimVarHCenter(eos, tau, e1, p, nb, nq, ns, vx1, vy1, vz1) ;
		f->getCell(ix, iy, iz-1)->getPrimVarHCenter(eos, tau, e0, p, nb, nq, ns, vx0, vy0, vz0) ;
		if(e1>0. && e0>0.){
		ut0=1.0/sqrt(1.0-vx0*vx0-vy0*vy0-vz0*vz0) ;
		ux0=ut0*vx0 ; uy0=ut0*vy0 ; uz0=ut0*vz0 ;
		ut1=1.0/sqrt(1.0-vx1*vx1-vy1*vy1-vz1*vz1) ;
		ux1=ut1*vx1 ; uy1=ut1*vy1 ; uz1=ut1*vz1 ;
		dmu[3][0]=0.25*(ut1*ut1-ut0*ut0)/uuu[0]/dz/(tau+0.5*dt) ; dmu[3][1]=0.25*(ux1*ux1-ux0*ux0)/uuu[1]/dz/(tau+0.5*dt) ;
		dmu[3][2]=0.25*(uy1*uy1-uy0*uy0)/uuu[2]/dz/(tau+0.5*dt) ; dmu[3][3]=0.25*(uz1*uz1-uz0*uz0)/uuu[3]/dz/(tau+0.5*dt) ;
		if(fabs(0.5*(ut1+ut0)/uuu[0])>UDIFF) dmu[3][0]=0.5*(ut1-ut0)/dz/(tau+0.5*dt) ;
		if(fabs(uuu[1])<VMIN || fabs(0.5*(ux1+ux0)/uuu[1])>UDIFF) dmu[3][1]=0.5*(ux1-ux0)/dz/(tau+0.5*dt) ;
		if(fabs(uuu[2])<VMIN || fabs(0.5*(uy1+uy0)/uuu[2])>UDIFF) dmu[3][2]=0.5*(uy1-uy0)/dz/(tau+0.5*dt) ;
		if(fabs(uuu[3])<VMIN || fabs(0.5*(uz1+uz0)/uuu[3])>UDIFF) dmu[3][3]=0.5*(uz1-uz0)/dz/(tau+0.5*dt) ;
		}else{ // matter-vacuum
		dmu[3][0]=dmu[3][1]=dmu[3][2]=dmu[3][3]=0. ;
		}
	// additional terms from Christoffel symbols :)
	dmu[3][0]+=uuu[3]/(tau-0.5*dt) ; dmu[3][3]+=uuu[0]/(tau-0.5*dt) ;
	// calculation of Z[mu][nu][lambda][rho]
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	for(int k=0; k<4; k++)
	for(int l=0; l<4; l++)
	 Z[i][j][k][l] = 0.0 ;
	// filling Z matrix
	for(int mu=0; mu<4; mu++)
	for(int nu=0; nu<4; nu++)
	for(int lam=0; lam<4; lam++)
	for(int rho=0; rho<4; rho++){
	 if(nu==rho) Z[mu][nu][lam][rho]+=0.5*(gmunu[mu][lam]-uuu[mu]*uuu[lam]) ;
	 if(mu==rho) Z[mu][nu][lam][rho]+=0.5*(gmunu[nu][lam]-uuu[nu]*uuu[lam]) ;
	 if(lam==rho) Z[mu][nu][lam][rho]-=(gmunu[mu][nu]-uuu[mu]*uuu[nu])/3.0 ;
	}
	// calculating sigma[mu][nu]
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++){
	 pi[i][j]=0.0 ;
	for(int k=0; k<4; k++)
	for(int l=0; l<4; l++){
	 pi[i][j]+=Z[i][j][k][l]*dmu[k][l]*2.0*etaS*s/5.068 ;
	}
	}
	Pi=-zetaS*s*(dmu[0][0]+dmu[1][1]+dmu[2][2]+dmu[3][3])/5.068 ; // fm^{-4} --> GeV/fm^3
	du = dmu[0][0]+dmu[1][1]+dmu[2][2]+dmu[3][3] ;
	//--------- debug part: trace check, diag check, transversality check
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++){
	 if(pi[i][j]!=0. && fabs(1.0-pi[j][i]/pi[i][j])>1e-10)
	  cout << "non-diag: " << pi[i][j] << "  " << pi[j][i] << endl ;
	}
	//for(int i=0; i<4; i++){
	//ccccc transversality check for pi[][]
	 //if(fabs(pi[i][0]*uuu[0]-pi[i][1]*uuu[1]-pi[i][2]*uuu[2]-pi[i][3]*tau*tau*uuu[3])/
	 //max(fabs(pi[i][0]),max(fabs(pi[i][1]),max(fabs(pi[i][2]),fabs(pi[i][3]))))>1e-2)
	  //cout << "non-transverse: "<<setw(14)<<pi[i][0]<<setw(14)<<pi[i][1]
	  //<<setw(14)<<pi[i][2]<<setw(14)<<pi[i][3]<<endl
	  //<<setw(14)<<uuu[0]<<setw(14)<<uuu[1]<<setw(14)<<uuu[2]<<setw(14)<<uuu[3]<<endl ;
	  //ccccc transversality check for Delta[][]
	  //double _sum = dmu[i][0]*uuu[0]-dmu[i][1]*uuu[1]-dmu[i][2]*uuu[2]-dmu[i][3]*uuu[3]/tau/tau ;
	  //double _sum = (gmunu[0][i]-uuu[i]*uuu[0])*uuu[0] - (gmunu[1][i]-uuu[i]*uuu[1])*uuu[1]
	  // - (gmunu[2][i]-uuu[i]*uuu[2])*uuu[2] -(gmunu[3][i]-uuu[i]*uuu[3])*uuu[3]/tau/tau ;
	 //if(fabs(_sum)/max(fabs(uuu[0]),max(fabs(uuu[1]),max(fabs(uuu[2]),fabs(uuu[3]))))>1e-5)
	  //cout << "non-transverse: dir= "<<setw(5)<<i<<setw(14)<<_sum<<endl
	  //<<setw(14)<<uuu[0]<<setw(14)<<uuu[1]<<setw(14)<<uuu[2]<<setw(14)<<uuu[3]<<endl ;
	//}
//    if(Pi!=0. && iy==f->getNY()/2 && iz==f->getNZ()/2)
//    cout << setw(14) << "trace pi/pi=" << setw(14) << 
//	(pi[0][0]-pi[1][1]-pi[2][2]-tau*tau*pi[3][3])/(pi[3][3]+1e-50) << endl ;

//----- debug - finding the maximal value of dmu
 //double dmumax=0. ;
 //for(int i=0; i<4; i++)
 //for(int j=0; j<4; j++)
  //if(fabs(dmu[i][j])>dmumax) dmumax=fabs(dmu[i][j]) ;
 //if(dmumax>200){
  //cout<<"dmu>200 triggered; cell # " << ix << " " << iy << " " << iz << ", dmu matrix below\n" ;
  //for(int i=0; i<4; i++){
  //for(int j=0; j<4; j++) cout << setw(14) << dmu[i][j] ;
  //cout << endl ;
  //}
 //}
}


void Hydro::setNSvalues()
{
	double e, p, nb, nq, ns, vx, vy, vz, piNS[4][4], PiNS, dmu[4][4], du ;
	for(int ix=0; ix<f->getNX(); ix++)
	for(int iy=0; iy<f->getNY(); iy++)
	for(int iz=0; iz<f->getNZ(); iz++){
	Cell *c = f->getCell(ix, iy, iz) ;
	c->getPrimVar(eos, tau, e, p, nb, nq, ns, vx, vy, vz) ;
	if(e<=0.) continue ;
	NSquant(ix, iy, iz, piNS, PiNS, dmu, du) ;
	for(int i=0; i<4; i++)
	for(int j=0; j<=i; j++)
		c->setpi(i,j,piNS[i][j]) ;
	c->setPi(PiNS) ;
	}
}



void Hydro::ISformal()
{
	double e, p, nb, nq, ns, vx, vy, vz, T, mub, muq, mus ;
	double piNS[4][4], PiNS, dmu[4][4], du, piH[4][4], PiH ;
	const double gmumu[4] = {1., -1., -1., -1.} ;
	
  // loop #1 (relaxation+source terms)
	for(int ix=0; ix<f->getNX(); ix++)
	for(int iy=0; iy<f->getNY(); iy++)
	for(int iz=0; iz<f->getNZ(); iz++){
	Cell *c = f->getCell(ix, iy, iz) ;
	c->getPrimVarHCenter(eos, tau-dt/2., e, p, nb, nq, ns, vx, vy, vz) ; // instead of getPrimVar()
	if(e<=0.){ // empty cell?
	for(int i=0; i<4; i++)
	for(int j=0; j<=i; j++){
	 c->setpiH0(i,j,0.0); c->setpi0(i,j,0.0) ;
     }
    c->setPiH0(0.0) ;
    c->setPi0(0.0) ;
//############# debug output: expansion rate, gamma derivative
//    if(iz==f->getNZ()/2){
//    foutNS<<setw(14)<<tau<<setw(4)<<ix<<setw(4)<<iy ;
//    for(int p=0; p<7; p++) foutNS<<setw(14)<<0 ;  foutNS<<endl ;}
    }else{ // non-empty cell
  // #### 0) pre-advection
	double xmH =-vx*dt/f->getDx()/2.0 ;
	double ymH =-vy*dt/f->getDy()/2.0 ;
	double zmH =-vz*dt/f->getDz()/(tau-0.5*dt)/2.0 ;
	double wxH [2] = {(1.-fabs(xmH)), fabs(xmH)} ;
	double wyH [2] = {(1.-fabs(ymH)), fabs(ymH)} ;
	double wzH [2] = {(1.-fabs(zmH)), fabs(zmH)} ;
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++){
		piH[i][j]=0.0;
	}
	PiH=0.0 ;
	for(int jx=0; jx<2; jx++)
	for(int jy=0; jy<2; jy++)
	for(int jz=0; jz<2; jz++){
	// pi,Pi-->half-step
	Cell* c1=f->getCell(ix+jx*sign(xmH),iy+jy*sign(ymH),iz+jz*sign(zmH)) ;
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++){
	 piH[i][j]+=wxH[jx]*wyH[jy]*wzH[jz]*c1->getpi(i,j) ;
    }
	PiH+=wxH[jx]*wyH[jy]*wzH[jz]*c1->getPi() ;
    }
  for(int i=0; i<4; i++)
	for(int j=0; j<=i; j++){
		c->setpi0(i,j, piH[i][j]);
	}
	c->setPi0(PiH) ;
	// #### 1) relaxation(pi)+source(pi) terms for half-step
	double gamma=1.0/sqrt(1.0-vx*vx-vy*vy-vz*vz) ;
	NSquant(ix, iy, iz, piNS, PiNS, dmu, du) ;
	eos->eos(e,nb,nq,ns,T,mub,muq,mus,p) ;
//############# debug output: expansion rate, gamma derivative
//    if(iz==f->getNZ()/2){
//    double s = eos->s(e,nb,nq,ns) ;
//    foutNS<<setw(14)<<tau<<setw(4)<<ix<<setw(4)<<iy<<setw(14)<<du<<setw(14)<<dmu[0][0]
//    <<setw(14)<<s<<setw(14)<<piNS[0][0]<<setw(14)<<piNS[0][1]<<setw(14)<<piNS[1][1]<<setw(14)<<piNS[3][3]<<endl ;}
//############# get relaxation times
    double taupi, tauPi ;
    trcoeff->getTau(T, taupi, tauPi) ;
//#############
	// relaxation term, piH,PiH-->half-step
	for(int i=0; i<4; i++)
	for(int j=0; j<=i; j++){
#ifdef FORMAL_SOLUTION
	 c->setpiH0(i,j, (c->getpi0(i,j)-piNS[i][j])*exp(-dt/2.0/gamma/taupi)+piNS[i][j]) ;
#else
	 c->setpiH0(i,j, c->getpi0(i,j)-(c->getpi0(i,j)-piNS[i][j])*dt/2.0/gamma/taupi) ;
#endif
    }
#ifdef FORMAL_SOLUTION
	c->setPiH0( (c->getPi0()-PiNS)*exp(-dt/2.0/gamma/tauPi)+PiNS ) ;
#else
	c->setPiH0( c->getPi0()-(c->getPi0()-PiNS)*dt/2.0/gamma/tauPi ) ;
#endif
    // sources from Christoffel symbols from \dot pi_munu
    double tau1=tau-dt*0.75 ;
    c->addpiH0(0,0, -2.*vz*c->getpi0(0,3)/tau1*dt/2. ) ; // *gamma/gamma
    c->addpiH0(3,3, -(2.*vz/tau1*c->getpi0(0,3))*dt/2. ) ;
    c->addpiH0(3,0, -(vz/tau1*c->getpi0(0,0) + vz/tau1*c->getpi0(3,3))*dt/2. ) ;
    c->addpiH0(1,0, -vz/tau1*c->getpi0(1,3)*dt/2. ) ;
    c->addpiH0(2,0, -vz/tau1*c->getpi0(2,3)*dt/2. ) ;
    c->addpiH0(3,1, -(vz/tau1*c->getpi0(0,1))*dt/2. ) ;
    c->addpiH0(3,2, -(vz/tau1*c->getpi0(0,2))*dt/2. ) ;
    // source from full IS equations (see  draft for the description)
    double u[4] ;
    u[0] = gamma ; u[1] = u[0]*vx ;
    u[2] = u[0]*vy; u[3]=u[0]*vz ;
    for(int i=0; i<4; i++)
    for(int j=0; j<=i; j++){
     c->addpiH0(i,j, -4./3.*c->getpi0(i,j)*du/gamma*dt/2. ) ;
     for(int k=0; k<4; k++) // terms to achieve better transverality to u^mu
     for(int l=0; l<4; l++)
      c->addpiH0(i,j,- c->getpi0(i,k)*u[j]*u[l]*dmu[l][k]*gmumu[k]/gamma*dt/2. 
                     - c->getpi0(j,k)*u[i]*u[l]*dmu[l][k]*gmumu[k]/gamma*dt/2. ) ;
	}
    c->addPiH0( -4./3.*c->getPi0()*du/gamma*dt/2. ) ;
	// 1) #### relaxation(piH)+source(piH) terms for full-step
	for(int i=0; i<4; i++)
	for(int j=0; j<=i; j++){
#ifdef FORMAL_SOLUTION
	 c->setpi0(i,j, (c->getpi0(i,j)-piNS[i][j])*exp(-dt/gamma/taupi)+piNS[i][j]) ;
#else
	 c->setpi0(i,j, c->getpi0(i,j)-(c->getpiH0(i,j)-piNS[i][j])*dt/gamma/taupi) ;
#endif
    }
#ifdef FORMAL_SOLUTION
	c->setPi0( (c->getPi0()-PiNS)*exp(-dt/gamma/tauPi)+PiNS ) ;
#else
	c->setPi0( c->getPi0()-(c->getPiH0()-PiNS)*dt/gamma/tauPi ) ;
#endif
	tau1=tau-dt*0.5 ;
    c->addpi0(0,0, -2.*vz/tau1*c->getpiH0(0,3)*dt ) ; // *gamma/gamma
    c->addpi0(3,3, -(2.*vz/tau1*c->getpiH0(0,3))*dt ) ;
    c->addpi0(3,0, -(vz/tau1*c->getpiH0(0,0) + vz/tau1*c->getpiH0(3,3))*dt ) ;
    c->addpi0(1,0, -vz/tau1*c->getpiH0(1,3)*dt ) ;
    c->addpi0(2,0, -vz/tau1*c->getpiH0(2,3)*dt ) ;
    c->addpi0(3,1, -(vz/tau1*c->getpiH0(0,1))*dt ) ;
    c->addpi0(3,2, -(vz/tau1*c->getpiH0(0,2))*dt ) ;
    // source from full IS equations (see draft for the description)
    for(int i=0; i<4; i++)
    for(int j=0; j<=i; j++){
     c->addpi0(i,j, -4./3.*c->getpiH0(i,j)*du/gamma*dt );
     for(int k=0; k<4; k++) // terms to achieve better transverality to u^mu
     for(int l=0; l<4; l++)
      c->addpi0(i,j,- c->getpiH0(i,k)*u[j]*u[l]*dmu[l][k]*gmumu[k]/gamma*dt 
                    - c->getpiH0(j,k)*u[i]*u[l]*dmu[l][k]*gmumu[k]/gamma*dt) ;
      //if((ix==38 || ix==62)&&iy==f->getNY()/2&&iz==f->getNZ()/2&&i==0&&j==2){
	//	  cout<<setw(5)<<ix<<setw(5)<<i<<setw(5)<<j<<setw(14)<<dpi<<setw(14)<<vx<<setw(14)<<vy<<endl ;
	//  }
	}
    c->addPi0( -4./3.*c->getPiH0()*du/gamma*dt ) ;
    } // end non-empty cell
    } // end loop #1

    // 3) -- post - advection, dt/2 (for full-step values only) ---
  for(int ix=0; ix<f->getNX(); ix++)
	for(int iy=0; iy<f->getNY(); iy++)
	for(int iz=0; iz<f->getNZ(); iz++){
	Cell *c = f->getCell(ix, iy, iz) ;
	c->getPrimVarHCenter(eos, tau-dt/2., e, p, nb, nq, ns, vx, vy, vz) ; // getPrimVar() before
	if(e<=0.) continue ;
	double xmH =-vx*dt/f->getDx()/2.0 ;
	double ymH =-vy*dt/f->getDy()/2.0 ;
	double zmH =-vz*dt/f->getDz()/(tau-0.5*dt)/2.0 ;
	double wxH [2] = {(1.-fabs(xmH)), fabs(xmH)} ;
	double wyH [2] = {(1.-fabs(ymH)), fabs(ymH)} ;
	double wzH [2] = {(1.-fabs(zmH)), fabs(zmH)} ;
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++){
		piH[i][j]=0.0;
	}
	PiH=0.0 ;
	for(int jx=0; jx<2; jx++)
	for(int jy=0; jy<2; jy++)
	for(int jz=0; jz<2; jz++){
	// pi,Pi-->full step, piH,PiH-->half-step
	Cell* c1=f->getCell(ix+jx*sign(xmH),iy+jy*sign(ymH),iz+jz*sign(zmH)) ;
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++){
	 piH[i][j]+=wxH[jx]*wyH[jy]*wzH[jz]*c1->getpi0(i,j) ;
	 //if(ix==45&&iy==32&&iz==0&&i==3&&j==3) // debug
	 //cout<<setw(10)<<ix+jx*sign(xm)<<setw(10)<<iy+jy*sign(ym)<<setw(14)<<iz+jz*sign(zm)
	 //<<c1->pi0[i][j]<<endl ;
    }
	PiH+=wxH[jx]*wyH[jy]*wzH[jz]*c1->getPi0() ;
    }
	//if(ix==45&&iy==32&&iz==0){ // debug
	 //cout<<"advection result (45,32,0):\n" ;
	 //for(int i=0; i<4; i++){
	  //for(int j=0; j<4; j++)
	  //cout << setw(14) << pi[i][j] ; cout<<endl ;
	 //}
	 //cout<<setw(14)<<wx[0]<<setw(14)<<wx[1]<<setw(14)<<wy[0]<<setw(14)<<wy[1]
	 //<<setw(14)<<wz[0]<<setw(14)<<wz[1]<<endl;
	 //cout << "du = " << du << " gamma = " << gamma << endl ;
	//}

	//--------- debug part: trace check, diag check, transversality check
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++){
	 if(piH[i][j]!=0. && fabs(1.0-piH[j][i]/piH[i][j])>1e-10)
	  cout << "non-diag: " << piH[i][j] << "  " << piH[j][i] << endl ;
	}
	//------ end debug
	//======= hydro applicability check (viscous corrections limiter):
	//double maxT0 = max(fabs((e+p)*vx*vx/(1.-vx*vx-vy*vy-vz*vz)+p),
	//   fabs((e+p)*vy*vy/(1.-vx*vx-vy*vy-vz*vz)+p)) ;
	double maxT0 = max((e+p)/(1.-vx*vx-vy*vy-vz*vz)-p,
	(e+p)*(vx*vx+vy*vy+vz*vz)/(1.-vx*vx-vy*vy-vz*vz)+p ) ;
	//double maxpi = max(fabs(pi[1][1]),fabs(pi[2][2])) ;
	double maxpi=0. ;
	for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	 if(fabs(piH[i][j])>maxpi) maxpi=fabs(piH[i][j]) ;
	bool rescaled=false ;
	if(maxT0/maxpi<1.0){
	   //cout<<"hydro_inappl: cell# "<<setw(6)<<ix<<setw(6)<<iy<<setw(6)<<iz<<endl
	   //<< (e+p)*vx*vx/(1.-vx*vx-vy*vy-vz*vz)+p<<"  "<<pi[1][1]<<"  "
	   //<< (e+p)*vy*vy/(1.-vx*vx-vy*vy-vz*vz)+p<<"  "<<pi[2][2]<< endl;
	   for(int i=0; i<4; i++)
	   for(int j=0; j<4; j++){
	    piH[i][j] = 0.1*piH[i][j]*maxT0/maxpi ; }
	    rescaled=true ;
    }
    if(fabs(PiH)>p){
	    if(PiH!=0.) PiH = 0.1*PiH/fabs(PiH)*p ;
	    rescaled=true ;
	}
  if(rescaled) c->setViscCorrCutFlag(maxT0/maxpi) ;
  else c->setViscCorrCutFlag(1.) ;
	//if(ix==22&&iy==25&&iz==0){
	 //cout<<"ISformal (22,25,0):  " << rescaled << endl ;
	 //for(int i=0; i<4; i++){
	  //for(int j=0; j<4; j++)
	  //cout << setw(14) << pi[i][j] ; cout<<endl ;
	 //}
	 //cout<<"ideal part:"<<setw(14)<<e<<setw(14)<<vx<<setw(14)<<vy<<setw(14)<<vz<<endl ;
    //}
    // updating to the new values
    for(int i=0; i<4; i++)
	for(int j=0; j<=i; j++){
		c->setpi(i,j, piH[i][j]);
		c->setpiH(i,j, c->getpiH0(i,j)) ;
	}
	c->setPi(PiH) ; c->setPiH(c->getPiH0()) ;
   } // advection loop (all cells)
}


void Hydro::setQfull()
{
	double uuu [4], Q [7] ; // the 4-velocity
	double e,p,nb,nq,ns,vx,vy,vz ;
	double gmunu [4][4] = {{1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1}} ;
	for(int ix=0; ix<f->getNX(); ix++)
	for(int iy=0; iy<f->getNY(); iy++)
	for(int iz=0; iz<f->getNZ(); iz++){
	 Cell *c = f->getCell(ix, iy, iz) ;
	 c->getQ(Q) ;
	 for(int i=0; i<7; i++) Q[i] = Q[i]/tau ; // remove factor "tau" from Q
	 c->getPrimVar(eos,tau,e,p,nb,nq,ns,vx,vy,vz) ;
	 uuu[0]=1.0/sqrt(1.0-vx*vx-vy*vy-vz*vz) ;
	 uuu[1]=vx*uuu[0];
	 uuu[2]=vy*uuu[0];
	 uuu[3]=vz*uuu[0];
	 for(int i=0; i<4; i++) Q[i] += - c->getPi()*(gmunu[0][i]-uuu[0]*uuu[i]) + c->getpi(0,i) ;
	 c->setQfull(Q) ;
	 //----debug
	 //if(ix==60&&iy==50&&iz==2){
	 //cout<<"setQfull:\n" ;
	 //cout<<setw(14)<<e<<setw(14)<<vx<<setw(14)<<vy<<setw(14)<<vz<<setw(14)<<c->Pi<<endl ;
	 //}
	}
}


// this procedure explicitly uses T_==0, X_==1, Y_==2, Z_==3
void Hydro::visc_flux(Cell *left, Cell *right, int direction)
{
 double flux[4] ;
 int ind2=0 ;
 double dxa ;
 // exit if noth cells are not full with matter
 if(left->getM(direction)<1. && right->getM(direction)<1.) return ;

 if(direction==X_) dxa=f->getDx() ;
 else if(direction==Y_) dxa=f->getDy() ;
 else if(direction==Z_) dxa=f->getDz()*(tau+0.5*dt) ;
 double e, p, nb, nq, ns, vxl, vyl, vzl, vxr, vyr, vzr ;
 // we need to know the velocities at both cell centers at (n+1/2) in order to interpolate to 
 // get the value at the interface
 left->getPrimVarHCenter(eos, tau-dt/2., e, p, nb, nq, ns, vxl, vyl, vzl) ;
 right->getPrimVarHCenter(eos, tau-dt/2., e, p, nb, nq, ns, vxr, vyr, vzr) ;
 vxl=0.5*(vxl+vxr) ; vyl=0.5*(vyl+vyr) ; vzl=0.5*(vzl+vzr) ;
 double v = sqrt(vxl*vxl+vyl*vyl+vzl*vzl) ;
 if(v>1.) { vxl=0.99*vxl/v ; vyl=0.99*vyl/v ; vzl=0.99*vzl/v ; }
 double gamma=1./sqrt(1.-v*v) ;
 double uuu[4] = {gamma, gamma*vxl, gamma*vyl, gamma*vzl} ;
 double gmumu [4] = {1.,-1.,-1.,-1.} ;
 if(direction==X_) ind2=1 ;
 else if(direction==Y_) ind2=2 ;
 else if(direction==Z_) ind2=3 ;
 for(int ind1=0; ind1<4; ind1++){
  flux[ind1] = 0.5*(left->getpiH(ind1,ind2)+right->getpiH(ind1,ind2)) ;
  if(ind1==ind2) flux[ind1] += -0.5*(left->getPiH()+right->getPiH())*gmumu[ind1] ; // gmunu is diagonal
  flux[ind1] += 0.5*(left->getPiH()+right->getPiH())*uuu[ind1]*uuu[ind2] ;
 }
 //---- debug
// if(fabs(flux[Z_])>1e+6){ cout << "flux bug!" << endl ;
 //if(left->getX()==50&&left->getY()==50&&left->getZ()==2){ cout << "flux here:" << endl ;
 //for(int i=0; i<4; i++) cout << setw(14) << flux[i] ; cout<<endl ;
 //for(int i=0; i<4; i++) cout << setw(14) << uuu[i] ; cout<<endl ;
 //for(int i=0; i<4; i++){
	 //for(int j=0; j<4; j++)
	 //cout << setw(14) << right->pi[i][j] ; cout<<endl ;}
 //}
 //---- end debug
 for(int i=0; i<4; i++) flux[i] = flux[i]*dt/dxa ;
	left->addFlux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], 0., 0., 0.) ;
	right->addFlux(flux[T_], flux[X_], flux[Y_], flux[Z_], 0., 0., 0.) ;
}


void Hydro::performStep(void)
{
  debugRiemann = false ; // turn off debug output

	f->updateM(tau, dt) ;

	tau_z = dt/2./log(1+dt/2./tau) ;

	//-----PREDICTOR-ideal
	for(int iy=0; iy<f->getNY(); iy++)
	for(int iz=0; iz<f->getNZ(); iz++)
	for(int ix=0; ix<f->getNX(); ix++){
		Cell *c = f->getCell(ix, iy, iz) ;
		c->saveQprev() ;
		c->clearFlux() ;
	}
		// X dir
	for(int iy=0; iy<f->getNY(); iy++)
	for(int iz=0; iz<f->getNZ(); iz++)
	for(int ix=0; ix<f->getNX()-1; ix++){
		hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix+1, iy, iz), X_, PREDICT) ;
	}
//	cout << "predictor X done\n" ;
		//Y dir
	for(int iz=0; iz<f->getNZ(); iz++)
	for(int ix=0; ix<f->getNX(); ix++)
	for(int iy=0; iy<f->getNY()-1; iy++){
		hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy+1, iz), Y_, PREDICT) ;
	}
//	cout << "predictor Y done\n" ;
			// Z dir
	for(int ix=0; ix<f->getNX(); ix++)
	for(int iy=0; iy<f->getNY(); iy++)
	for(int iz=0; iz<f->getNZ()-1; iz++){
		hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy, iz+1), Z_, PREDICT) ;
	}
//	cout << "predictor Z done\n" ;

	for(int iy=0; iy<f->getNY(); iy++)
	for(int iz=0; iz<f->getNZ(); iz++)
	for(int ix=0; ix<f->getNX(); ix++){
		Cell *c = f->getCell(ix, iy, iz) ;
		source_step(ix, iy, iz, PREDICT) ;
		c->updateQtoQhByFlux() ;
		c->clearFlux() ;
	}

	//----CORRECTOR-ideal

	tau_z = dt/log(1+dt/tau) ;
		// X dir
	for(int iy=0; iy<f->getNY(); iy++)
	for(int iz=0; iz<f->getNZ(); iz++)
	for(int ix=0; ix<f->getNX()-1; ix++){
		hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix+1, iy, iz), X_, CORRECT) ;
	}
//	cout << "corrector X done\n" ;
		//Y dir
	for(int iz=0; iz<f->getNZ(); iz++)
	for(int ix=0; ix<f->getNX(); ix++)
	for(int iy=0; iy<f->getNY()-1; iy++){
		hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy+1, iz), Y_, CORRECT) ;
	}
//	cout << "corrector Y done\n" ;
			// Z dir
	for(int ix=0; ix<f->getNX(); ix++)
	for(int iy=0; iy<f->getNY(); iy++)
	for(int iz=0; iz<f->getNZ()-1; iz++){
		hlle_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy, iz+1), Z_, CORRECT) ;
	}
//	cout << "corrector Z done\n" ;

	for(int iy=0; iy<f->getNY(); iy++)
	for(int iz=0; iz<f->getNZ(); iz++)
	for(int ix=0; ix<f->getNX(); ix++)
	{
		Cell *c = f->getCell(ix, iy, iz) ;
		source_step(ix, iy, iz, CORRECT) ;
		c->updateByFlux() ;
		c->clearFlux() ;
	}
	tau += dt ;
	f->correctImagCells() ;

	//===== viscous part ======
if(trcoeff->isViscous()){
	setQfull() ; // set the values of full Q^\mu including visc corrections
	ISformal() ; // evolution of viscous quantities according to IS equations


		// X dir
	for(int iy=0; iy<f->getNY(); iy++)
	for(int iz=0; iz<f->getNZ(); iz++)
	for(int ix=0; ix<f->getNX()-1; ix++){
		visc_flux(f->getCell(ix, iy, iz), f->getCell(ix+1, iy, iz), X_) ;
	}
//	cout << "visc_flux X done\n" ;
		//Y dir
	for(int iz=0; iz<f->getNZ(); iz++)
	for(int ix=0; ix<f->getNX(); ix++)
	for(int iy=0; iy<f->getNY()-1; iy++){
		visc_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy+1, iz), Y_) ;
	}
//	cout << "visc_flux Y done\n" ;
			// Z dir
	for(int ix=0; ix<f->getNX(); ix++)
	for(int iy=0; iy<f->getNY(); iy++)
	for(int iz=0; iz<f->getNZ()-1; iz++){
		visc_flux(f->getCell(ix, iy, iz), f->getCell(ix, iy, iz+1), Z_) ;
	}

	for(int iy=0; iy<f->getNY(); iy++)
	for(int iz=0; iz<f->getNZ(); iz++)
	for(int ix=0; ix<f->getNX(); ix++){
		visc_source_step(ix, iy, iz) ;
		f->getCell(ix, iy, iz)->updateQfullByFlux() ;
		f->getCell(ix, iy, iz)->correctQideal(eos,tau) ;
		//if(ix==60&&iy==50&&iz==2){
		 //double e,p,nb,nq,ns,vx,vy,vz;
		 //f->getCell(ix, iy, iz)->getPrimVar(eos, e, p, nb, nq, ns, vx, vy, vz) ;
	     //cout<<"correctQideal:\n" ;
	     //cout<<setw(14)<<e<<setw(14)<<vx<<setw(14)<<vy<<setw(14)<<vz<<endl ;
	    //}
	}
}else{ // end viscous part
    
}
	//==== finishing work ====
	f->correctImagCellsFull() ;
}
