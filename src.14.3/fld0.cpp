#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "inc.h"
#include "fld.h"
#include "cll.h"
#include "eos.h"
#include "trancoeff.h"
#include "cornelius.h"

#define OUTPI

// change to hadron EoS (e.g. Laine) to calculate v,T,mu at the surface
#define SWAP_EOS

using namespace std ;

// returns the velocities in cartesian coordinates, fireball rest frame. Y=longitudinal rapidity of fluid
void Fluid::getCMFvariables(Cell *c, double tau, double &e, double &nb, double &nq, double &ns, double &vx, double &vy, double &Y)
{
	double p, vz ;
	c->getPrimVar(eos, tau, e, p, nb, nq, ns, vx, vy, vz) ;
	double eta = getZ(c->getZ()) ;
//	Y = eta + TMath::ATanH(vz) ;
	Y = eta + 1./2.*log((1.+vz)/(1.-vz)) ;
	vx = vx * cosh(Y-eta)/cosh(Y) ;
	vy = vy * cosh(Y-eta)/cosh(Y) ;
}



Fluid::Fluid(EoS *_eos, EoS *_eosH, TransportCoeff *_trcoeff, int _nx, int _ny, int _nz, double _minx, double _maxx, double _miny, double _maxy, double _minz, double _maxz, double _dt, double eCrit)
{
	eos = _eos ;
  eosH = _eosH ;
	trcoeff = _trcoeff ;
	nx = _nx ;
	ny = _ny ;
	nz = _nz ;
	minx = _minx ;
	maxx = _maxx ;
	miny = _miny ;
	maxy = _maxy ;
	minz = _minz ;
	maxz = _maxz ;
	dx = (maxx-minx)/(nx-1) ;
	dy = (maxy-miny)/(ny-1) ;
	dz = (maxz-minz)/(nz-1) ;
  dt = _dt ;

	cell = new Cell [nx*ny*nz] ;

	cell0 = new Cell ;
	cell0 ->setPrimVar(eos, 1.0, 0., 0., 0., 0., 0., 0., 0.) ; // tau is not important here, since *0
	cell0 ->setAllM(0.) ;

	for(int ix=0; ix<nx; ix++)
	for(int iy=0; iy<ny; iy++)
	for(int iz=0; iz<nz; iz++)
	{
		getCell(ix,iy,iz)->setPrev(X_, getCell(ix-1, iy, iz)) ;
		getCell(ix,iy,iz)->setNext(X_, getCell(ix+1, iy, iz)) ;
		getCell(ix,iy,iz)->setPrev(Y_, getCell(ix, iy-1, iz)) ;
		getCell(ix,iy,iz)->setNext(Y_, getCell(ix, iy+1, iz)) ;
		getCell(ix,iy,iz)->setPrev(Z_, getCell(ix, iy, iz-1)) ;
		getCell(ix,iy,iz)->setNext(Z_, getCell(ix, iy, iz+1)) ;
		getCell(ix,iy,iz)->setPos(ix, iy, iz) ;
	}

	output_nt = 0 ;
	output_nx = 0 ;
	output_ny = 0 ;
  
//---- Cornelius init
  double arrayDx [4] = {dt, dx, dy, dz} ;
  cornelius = new Cornelius ;
  cornelius->init(4, eCrit, arrayDx) ;
  vEff = 0. ;
  EtotSurf = 0.0 ;
}


Fluid::~Fluid()
{
	delete [] cell ;
	delete cell0 ;
}


void Fluid::initOutput(char *dir, int maxstep, double tau0, int cmpr2dOut)
{
//    directory = dir ;
      compress2dOut = cmpr2dOut ;
    cout << "maxstep=" << maxstep << endl ;
      char command [255] ;
      sprintf(command, "mkdir -p %s",dir) ;
      system(command) ;
      string outx = dir ; outx.append("/outx.dat");
      string outxvisc = dir ; outxvisc.append("/outx.visc.dat");
      string outyvisc = dir ; outyvisc.append("/outy.visc.dat");
      string outdiagvisc = dir ; outdiagvisc.append("/diag.visc.dat");
      string outy = dir ; outy.append("/outy.dat");
      string outdiag = dir ; outdiag.append("/outdiag.dat");
      string outz = dir ; outz.append("/outz.dat");
      string outaniz = dir ; outaniz.append("/out.aniz.dat");
      string out2d = dir ; out2d.append("/out2D.dat");
      string outfreeze = dir ; outfreeze.append("/freezeout.dat");
      foutx.open(outx.c_str()) ; fouty.open(outy.c_str()); foutz.open(outz.c_str()) ;
      foutdiag.open(outdiag.c_str()) ;
      fout2d.open(out2d.c_str()) ;
      foutxvisc.open(outxvisc.c_str()) ;
      foutyvisc.open(outyvisc.c_str()) ;
      foutdiagvisc.open(outdiagvisc.c_str()) ;
      fout_aniz.open(outaniz.c_str()) ;
      ffreeze.open(outfreeze.c_str()) ;
      //################################################################
      // important remark. for correct diagonal output, nx=ny must hold.
      //################################################################
      foutxvisc << maxstep + 1 << "  " << getNX() << endl ;
      foutyvisc << maxstep + 1 << "  " << getNY() << endl ;
      foutdiagvisc << maxstep + 1 << "  " << getNX() << endl ;
      foutx << "# " << maxstep + 1 << "  " << getNX() << endl ;
      fouty << "# "  << maxstep + 1 << "  " << getNY() << endl ;
      foutdiag << "# "  << maxstep + 1 << "  " << getNX() << endl ;
      foutz << "# "  << maxstep + 1 << "  " << getNZ() << endl ;
      fout2d << " " << maxstep+1  << "  " << (getNX()-5)+1 << "  " << (getNY()-5)+1 << endl ;
      fout2d << tau0 << "  " << tau0+0.05*maxstep << "  " << getX(2) << "  " << getX(getNX()-3) << "  "
	  << getY(2) << "  " << getY(getNY()-3) << endl ;
      outputPDirections(tau0);
}


void Fluid::correctImagCells(void)
{
	double Q [7] ;
// Z
	for(int ix=0; ix<nx; ix++)
	for(int iy=0; iy<ny; iy++){
		// left boundary
		getCell(ix,iy,2)->getQ(Q) ;
		getCell(ix,iy,1)->setQ(Q) ;
		getCell(ix,iy,0)->setQ(Q) ;
		// right boundary
		getCell(ix,iy,nz-3)->getQ(Q) ;
		getCell(ix,iy,nz-2)->setQ(Q) ;
		getCell(ix,iy,nz-1)->setQ(Q) ;
	}
// Y
	for(int ix=0; ix<nx; ix++)
	for(int iz=0; iz<nz; iz++){
		// left boundary
		getCell(ix,2,iz)->getQ(Q) ;
		getCell(ix,1,iz)->setQ(Q) ;
		getCell(ix,0,iz)->setQ(Q) ;
		// right boundary
		getCell(ix,ny-3,iz)->getQ(Q) ;
		getCell(ix,ny-2,iz)->setQ(Q) ;
		getCell(ix,ny-1,iz)->setQ(Q) ;
	}
// X
	for(int iy=0; iy<ny; iy++)
	for(int iz=0; iz<nz; iz++){
		// left boundary
		getCell(2,iy,iz)->getQ(Q) ;
		getCell(1,iy,iz)->setQ(Q) ;
		getCell(0,iy,iz)->setQ(Q) ;
		// right boundary
		getCell(nx-3,iy,iz)->getQ(Q) ;
		getCell(nx-2,iy,iz)->setQ(Q) ;
		getCell(nx-1,iy,iz)->setQ(Q) ;
	}
}


void Fluid::correctImagCellsFull(void)
{
	double Q [7], _pi[4][4], _Pi ;
// Z
	for(int ix=0; ix<nx; ix++)
	for(int iy=0; iy<ny; iy++){
		// left boundary
		getCell(ix,iy,2)->getQ(Q) ;
		getCell(ix,iy,1)->setQ(Q) ;
		getCell(ix,iy,0)->setQ(Q) ;
		for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
			_pi[i][j]=getCell(ix,iy,2)->getpi(i,j) ;
		_Pi=getCell(ix,iy,2)->getPi() ;
		
		for(int i=0; i<4; i++)
		for(int j=0; j<=i; j++){
			getCell(ix,iy,0)->setpi(i,j,_pi[i][j]) ;
			getCell(ix,iy,1)->setpi(i,j,_pi[i][j]) ;
		}
		getCell(ix,iy,0)->setPi(_Pi) ;
		getCell(ix,iy,1)->setPi(_Pi) ;
		// right boundary
		getCell(ix,iy,nz-3)->getQ(Q) ;
		getCell(ix,iy,nz-2)->setQ(Q) ;
		getCell(ix,iy,nz-1)->setQ(Q) ;
		for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
			_pi[i][j]=getCell(ix,iy,nz-3)->getpi(i,j) ;
		_Pi=getCell(ix,iy,nz-3)->getPi() ;
		
		for(int i=0; i<4; i++)
		for(int j=0; j<=i; j++){
			getCell(ix,iy,nz-2)->setpi(i,j,_pi[i][j]) ;
			getCell(ix,iy,nz-1)->setpi(i,j,_pi[i][j]) ;
		}
		getCell(ix,iy,nz-2)->setPi(_Pi) ;
		getCell(ix,iy,nz-1)->setPi(_Pi) ;
	}
// Y
	for(int ix=0; ix<nx; ix++)
	for(int iz=0; iz<nz; iz++){
		// left boundary
		getCell(ix,2,iz)->getQ(Q) ;
		getCell(ix,1,iz)->setQ(Q) ;
		getCell(ix,0,iz)->setQ(Q) ;
		for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
			_pi[i][j]=getCell(ix,2,iz)->getpi(i,j) ;
		_Pi=getCell(ix,2,iz)->getPi() ;
		
		for(int i=0; i<4; i++)
		for(int j=0; j<=i; j++){
			getCell(ix,0,iz)->setpi(i,j,_pi[i][j]) ;
			getCell(ix,1,iz)->setpi(i,j,_pi[i][j]) ;
		}
		getCell(ix,0,iz)->setPi(_Pi) ;
		getCell(ix,1,iz)->setPi(_Pi) ;
		// right boundary
		getCell(ix,ny-3,iz)->getQ(Q) ;
		getCell(ix,ny-2,iz)->setQ(Q) ;
		getCell(ix,ny-1,iz)->setQ(Q) ;
		for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
			_pi[i][j]=getCell(ix,ny-3,iz)->getpi(i,j) ;
		_Pi=getCell(ix,ny-3,iz)->getPi() ;
		
		for(int i=0; i<4; i++)
		for(int j=0; j<=i; j++){
			getCell(ix,ny-2,iz)->setpi(i,j,_pi[i][j]) ;
			getCell(ix,ny-1,iz)->setpi(i,j,_pi[i][j]) ;
		}
		getCell(ix,ny-2,iz)->setPi(_Pi) ;
		getCell(ix,ny-1,iz)->setPi(_Pi) ;
	}
// X
	for(int iy=0; iy<ny; iy++)
	for(int iz=0; iz<nz; iz++){
		// left boundary
		getCell(2,iy,iz)->getQ(Q) ;
		getCell(1,iy,iz)->setQ(Q) ;
		getCell(0,iy,iz)->setQ(Q) ;
		for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
			_pi[i][j]=getCell(2,iy,iz)->getpi(i,j) ;
		_Pi=getCell(2,iy,iz)->getPi() ;
		
		for(int i=0; i<4; i++)
		for(int j=0; j<=i; j++){
			getCell(0,iy,iz)->setpi(i,j,_pi[i][j]) ;
			getCell(1,iy,iz)->setpi(i,j,_pi[i][j]) ;
		}
		getCell(0,iy,iz)->setPi(_Pi) ;
		getCell(1,iy,iz)->setPi(_Pi) ;
		// right boundary
		getCell(nx-3,iy,iz)->getQ(Q) ;
		getCell(nx-2,iy,iz)->setQ(Q) ;
		getCell(nx-1,iy,iz)->setQ(Q) ;
		for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
			_pi[i][j]=getCell(nx-3,iy,iz)->getpi(i,j) ;
		_Pi=getCell(nx-3,iy,iz)->getPi() ;
		
		for(int i=0; i<4; i++)
		for(int j=0; j<=i; j++){
			getCell(nx-2,iy,iz)->setpi(i,j,_pi[i][j]) ;
			getCell(nx-1,iy,iz)->setpi(i,j,_pi[i][j]) ;
		}
		getCell(nx-2,iy,iz)->setPi(_Pi) ;
		getCell(nx-1,iy,iz)->setPi(_Pi) ;
	}
}

/*
void Fluid::updateM(double tau, double dt)
{
	for(int iy=0; iy<getNY(); iy++)
	for(int iz=0; iz<getNZ(); iz++)
	for(int ix=0; ix<getNX(); ix++){
//		if(getCell(ix,iy,iz)->getM(X_)<1. && (getCell(ix-1,iy,iz)->getM(X_)==1. || getCell(ix+1,iy,iz)->getM(X_)==1.))
		if(getCell(ix,iy,iz)->getM(X_)<1. && (getCell(ix-1,iy,iz)->getLM()==1. || getCell(ix+1,iy,iz)->getLM()==1.))
		{ getCell(ix,iy,iz)->addM(X_, dt/dx) ; ix++; }
	}

	for(int iz=0; iz<getNZ(); iz++)
	for(int ix=0; ix<getNX(); ix++)
	for(int iy=0; iy<getNY(); iy++){
//		if(getCell(ix,iy,iz)->getM(Y_)<1. && (getCell(ix,iy-1,iz)->getM(Y_)==1. || getCell(ix,iy+1,iz)->getM(Y_)==1.))
		if(getCell(ix,iy,iz)->getM(Y_)<1. && (getCell(ix,iy-1,iz)->getLM()==1. || getCell(ix,iy+1,iz)->getLM()==1.))
		{ getCell(ix,iy,iz)->addM(Y_, dt/dy) ; iy++; }
	}


	for(int ix=0; ix<getNX(); ix++)
	for(int iy=0; iy<getNY(); iy++)
	for(int iz=0; iz<getNZ(); iz++){
//		if(getCell(ix,iy,iz)->getM(Z_)<1. && (getCell(ix,iy,iz-1)->getM(Z_)==1. || getCell(ix,iy,iz+1)->getM(Z_)==1.))
		if(getCell(ix,iy,iz)->getM(Z_)<1. && (getCell(ix,iy,iz-1)->getLM()==1. || getCell(ix,iy,iz+1)->getLM()==1.))
		{ getCell(ix,iy,iz)->addM(Z_, dt/dz/tau) ; iz++; }
	}
}
*/

void Fluid::updateM(double tau, double dt)
{
	for(int ix=0; ix<getNX(); ix++)
	for(int iy=0; iy<getNY(); iy++)
	for(int iz=0; iz<getNZ(); iz++){
	Cell* c = getCell(ix,iy,iz) ;
	c->setDM(X_, 0.) ;
	c->setDM(Y_, 0.) ;
	c->setDM(Z_, 0.) ;
	if(getCell(ix,iy,iz)->getLM()<1.){
		if(getCell(ix+1,iy,iz)->getM(X_)>=1. || getCell(ix-1,iy,iz)->getM(X_)>=1.) c->setDM(X_, dt/dx) ;
		if(getCell(ix,iy+1,iz)->getM(Y_)>=1. || getCell(ix,iy-1,iz)->getM(Y_)>=1.) c->setDM(Y_, dt/dy) ;
		if(getCell(ix,iy,iz+1)->getM(Z_)>=1. || getCell(ix,iy,iz-1)->getM(Z_)>=1.) c->setDM(Z_, dt/dz/tau) ;

		if(c->getDM(X_)==0. && c->getDM(Y_)==0.){
			if(getCell(ix+1,iy+1,iz)->getLM()>=1. || getCell(ix+1,iy-1,iz)->getLM()>=1. ||
			getCell(ix-1,iy+1,iz)->getLM()>=1. || getCell(ix-1,iy-1,iz)->getLM()>=1.){
				c->setDM(X_, 0.707*dt/dx) ;
				c->setDM(Y_, 0.707*dt/dy) ;
			}
		}
	} //if
	}

	for(int ix=0; ix<getNX(); ix++)
	for(int iy=0; iy<getNY(); iy++)
	for(int iz=0; iz<getNZ(); iz++){
		Cell* c = getCell(ix,iy,iz) ;
		c->addM(X_,c->getDM(X_)) ;
		c->addM(Y_,c->getDM(Y_)) ;
		c->addM(Z_,c->getDM(Z_)) ;
	}
}


void Fluid::outputPDirections(double tau)
{
  outputGnuplot(tau) ;
  return ;
  
	double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz ;

	// X direction
	foutx << tau << endl ;
	foutxvisc << tau << endl ;
	for(int ix=0; ix<nx; ix++){
	double x = getX(ix) ;
	Cell *c = getCell(ix,ny/2,nz/2) ;
	getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz) ;
	eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
	foutx << setw(6) << x << setw(14) << vx << setw(14) << vy << setw(14) << e << setw(14) << nb <<  setw(14) << t <<  setw(14) << mub << endl ;
	//foutx << setw(6) << " " << setw(14) << (e+p)*vx*vx/(1.-vx*vx-vy*vy-vz*vz)+p << setw(14) << c->pi[1][1] << endl ;
	//if(ix==36){
	//c->getQ(Q) ;
	//foutx << "(36,50,2)Qideal  " ; for(int i=0; i<4; i++) foutx<<setw(14)<<Q[i]; foutx<<endl ;
	//c->getQfull(Q) ;
	//foutx << "(36,50,2)Qfull  " ; for(int i=0; i<4; i++) foutx<<setw(14)<<Q[i]; foutx<<endl ;
	//}
	foutxvisc << setw(6) << x ;
	foutxvisc << setw(14) << c->getpi(0,0) << setw(14) << c->getpi(0,1) << setw(14) << c->getpi(0,2); 
	foutxvisc << setw(14) <<c->getpi(0,3) << setw(14) << c->getpi(1,1) << setw(14) << c->getpi(1,2); 
	foutxvisc << setw(14) <<c->getpi(1,3) << setw(14) << c->getpi(2,2) << setw(14) << c->getpi(2,3); 
	foutxvisc << setw(14) <<c->getpi(3,3) << setw(14) << c->getPi() << endl;
	foutxvisc << setw(20) << "anomalies tr pi/pi=" << setw(14) << 
	(c->getpi(0,0)-c->getpi(1,1)-c->getpi(2,2)-c->getpi(3,3))/(c->getpi(1,1)+c->getpi(2,2)+1e-50) <<
	setw(10) << "transv." <<
	(c->getpi(0,0)-c->getpi(0,1)*vx-c->getpi(0,2)*vy-c->getpi(0,3)*vz)/(fabs(c->getpi(3,3))+1e-50) << endl ;
//	foutxvisc << setw(14) << "tmunu_id/corr" << setw(14) << (e+p)/(1.-vx*vx-vy*vy-vz*vz)-p
//	 << setw(14) << c->pi[3][3] << setw(14) << c->Pi*(1.-1./(1.-vx*vx-vy*vy-vz*vz)) << endl ;
	}
	foutx << endl ;
	foutxvisc << endl ;

	// Y direction
	fouty << tau << endl ;
  foutyvisc << tau << endl ;
	for(int iy=0; iy<ny; iy++){
	double y = getY(iy) ;
  Cell *c = getCell(nx/2,iy,nz/2) ;
	getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz) ;
	eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
	fouty << setw(6) << y << setw(14) << vy << setw(14) << vx << setw(14) << e << setw(14) << nb <<  setw(14) << t <<  setw(14) << mub << endl ;
//  foutyvisc << setw(6) << y ;
//	foutyvisc << setw(14) << c->pi[0][0] << setw(14) << c->pi[0][1] << setw(14) << c->pi[0][2]; 
//	foutyvisc << setw(14) <<c->pi[0][3] << setw(14) << c->pi[1][1] << setw(14) << c->pi[1][2]; 
//	foutyvisc << setw(14) <<c->pi[1][3] << setw(14) << c->pi[2][2] << setw(14) << c->pi[2][3]; 
//	foutyvisc << setw(14) <<c->pi[3][3] << setw(14) << c->Pi << endl;
//	foutyvisc << setw(20) << "anomalies tr pi/pi=" << setw(14) << 
//	(c->pi[0][0]-c->pi[1][1]-c->pi[2][2]-c->pi[3][3])/(c->pi[1][1]+c->pi[2][2]+1e-50) <<
//	setw(10) << "transv." <<
//	(c->pi[0][0]-c->pi[0][1]*vx-c->pi[0][2]*vy-c->pi[0][3]*vz)/(fabs(c->pi[3][3])+1e-50) << endl ;
	}
	fouty << endl ;
  foutyvisc << endl ;
  
  // diagonal
	foutdiag << tau << endl ;
  foutdiagvisc << tau << endl ;
	for(int ix=0; ix<nx; ix++){
	double x = getY(ix) ;
  Cell *c = getCell(ix,ix,nz/2) ;
	getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz) ;
	eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
	foutdiag << setw(14) << sqrt(2.)*x << setw(14) << vx << setw(14) << vy << setw(14) << e << setw(14) << nb <<  setw(14) << t <<  setw(14) << mub << endl ;
//  foutdiagvisc << setw(14) << sqrt(2.)*x ;
//	foutdiagvisc << setw(14) << c->pi[0][0] << setw(14) << c->pi[0][1] << setw(14) << c->pi[0][2]; 
//	foutdiagvisc << setw(14) <<c->pi[0][3] << setw(14) << c->pi[1][1] << setw(14) << c->pi[1][2]; 
//	foutdiagvisc << setw(14) <<c->pi[1][3] << setw(14) << c->pi[2][2] << setw(14) << c->pi[2][3]; 
//	foutdiagvisc << setw(14) <<c->pi[3][3] << setw(14) << c->Pi << endl;
//	foutdiagvisc << setw(20) << "anomalies tr pi/pi=" << setw(14) << 
//	(c->pi[0][0]-c->pi[1][1]-c->pi[2][2]-c->pi[3][3])/(c->pi[1][1]+c->pi[2][2]+1e-50) <<
//	setw(10) << "transv." <<
//	(c->pi[0][0]-c->pi[0][1]*vx-c->pi[0][2]*vy-c->pi[0][3]*vz)/(fabs(c->pi[3][3])+1e-50) << endl ;
	}
	foutdiag << endl ;
  foutdiagvisc << endl ;

	// Z direction
	foutz << tau << endl ;
	for(int iz=0; iz<nz; iz++){
	double z = getZ(iz) ;
	getCMFvariables(getCell(nx/2,ny/2,iz), tau, e, nb, nq, ns, vx, vy, vz) ;
	eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
	foutz << setw(6) << z << setw(14) << vz << setw(14) << e << setw(14) << nb <<  setw(14) << t <<  setw(14) << mub << endl ;
	}
	foutz << endl ;

	// 2d output - for escape4 calculation
	for(int ix=2; ix<nx-2; ix++)
	for(int iy=2; iy<ny-2; iy++){
	double x = getX(ix) ;
	double y = getY(iy) ;
	getCMFvariables(getCell(ix,iy,nz/2), tau, e, nb, nq, ns, vx, vy, vz) ;
  eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
	fout2d << " " << vx << " " << vy << " " << e << " " << t << " " << 0.0 << endl ;
	}
}


void Fluid::outputGnuplot(double tau)
{
	double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz ;

	// X direction
	for(int ix=0; ix<nx; ix++){
	double x = getX(ix) ;
	Cell *c = getCell(ix,ny/2,nz/2) ;
	getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz) ;
	eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
	foutx << setw(14) << tau << setw(14) << x << setw(14) << vx << setw(14) << vy << setw(14) << e
  << setw(14) << nb <<  setw(14) << t <<  setw(14) << mub ;
	foutx << setw(14) << c->getpi(0,0) << setw(14) << c->getpi(0,1) << setw(14) << c->getpi(0,2); 
	foutx << setw(14) <<c->getpi(0,3) << setw(14) << c->getpi(1,1) << setw(14) << c->getpi(1,2); 
	foutx << setw(14) <<c->getpi(1,3) << setw(14) << c->getpi(2,2) << setw(14) << c->getpi(2,3); 
	foutx << setw(14) <<c->getpi(3,3) << setw(14) << c->getPi() << endl;
	}
	foutx << endl ;

	// Y direction
	for(int iy=0; iy<ny; iy++){
	double y = getY(iy) ;
  Cell *c = getCell(nx/2,iy,nz/2) ;
	getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz) ;
	eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
	fouty << setw(14) << tau << setw(14) << y << setw(14) << vy << setw(14) << vx << setw(14) << e
  << setw(14) << nb <<  setw(14) << t <<  setw(14) << mub ;
	fouty << setw(14) << c->getpi(0,0) << setw(14) << c->getpi(0,1) << setw(14) << c->getpi(0,2); 
	fouty << setw(14) <<c->getpi(0,3) << setw(14) << c->getpi(1,1) << setw(14) << c->getpi(1,2); 
	fouty << setw(14) <<c->getpi(1,3) << setw(14) << c->getpi(2,2) << setw(14) << c->getpi(2,3); 
	fouty << setw(14) <<c->getpi(3,3) << setw(14) << c->getPi() << endl;
	}
	fouty << endl ;

  
  // diagonal
	for(int ix=0; ix<nx; ix++){
	double x = getY(ix) ;
  Cell *c = getCell(ix,ix,nz/2) ;
	getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz) ;
	eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
	foutdiag << setw(14) << tau << setw(14) << sqrt(2.)*x << setw(14) << vx 
  << setw(14) << vy << setw(14) << e << setw(14) << nb <<  setw(14) << t <<  setw(14) << mub << endl ;
	foutdiag << setw(14) << c->getpi(0,0) << setw(14) << c->getpi(0,1) << setw(14) << c->getpi(0,2); 
	foutdiag << setw(14) <<c->getpi(0,3) << setw(14) << c->getpi(1,1) << setw(14) << c->getpi(1,2); 
	foutdiag << setw(14) <<c->getpi(1,3) << setw(14) << c->getpi(2,2) << setw(14) << c->getpi(2,3); 
	foutdiag << setw(14) <<c->getpi(3,3) << setw(14) << c->getPi() << endl;
	}
	foutdiag << endl ;

	// Z direction
	for(int iz=0; iz<nz; iz++){
	double z = getZ(iz) ;
	Cell *c = getCell(nx/2,ny/2,iz) ;
	getCMFvariables(getCell(nx/2,ny/2,iz), tau, e, nb, nq, ns, vx, vy, vz) ;
	eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
	foutz << setw(14) << tau << setw(14) << z << setw(14) << vz << setw(14) << vx << setw(14) << e
  << setw(14) << nb <<  setw(14) << t <<  setw(14) << mub ;
	foutz << setw(14) << c->getpi(0,0) << setw(14) << c->getpi(0,1) << setw(14) << c->getpi(0,2); 
	foutz << setw(14) <<c->getpi(0,3) << setw(14) << c->getpi(1,1) << setw(14) << c->getpi(1,2); 
	foutz << setw(14) <<c->getpi(1,3) << setw(14) << c->getpi(2,2) << setw(14) << c->getpi(2,3); 
	foutz << setw(14) <<c->getpi(3,3) << setw(14) << c->getPi() << endl;
	}
	foutz << endl ;

}


// unput: geom. rapidity + velocities in Bjorken frame, --> output: velocities in lab.frame
void transformToLab(double eta, double &vx, double &vy, double &vz)
{
	const double Y = eta + 1./2.*log((1.+vz)/(1.-vz)) ;
	vx = vx * cosh(Y-eta)/cosh(Y) ;
	vy = vy * cosh(Y-eta)/cosh(Y) ;
  vz = tanh(Y) ;
}


void Fluid::outputSurface(double tau)
{
 double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz, Q[7] ;
 double E = 0., Efull = 0., Px=0., vt_num=0., vt_den=0., vxvy_num=0., vxvy_den=0., pi0x_num=0., pi0x_den=0.,
 txxyy_num=0., txxyy_den=0. ;
 double eta=0 ;
//-- Cornelius: allocating memory for corner points
  double ****ccube = new double***[2];
  for (int i1=0; i1 < 2; i1++) {
    ccube[i1] = new double**[2];
    for (int i2=0; i2 < 2; i2++) {
      ccube[i1][i2] = new double*[2];
      for (int i3=0; i3 < 2; i3++) {
        ccube[i1][i2][i3] = new double[2];
      }
    }
  }
//----end Cornelius
#ifdef SWAP_EOS
 swap(eos, eosH) ;
#endif
 fout2d << endl ;
 for(int ix=2; ix<nx-2; ix++)
 for(int iy=2; iy<ny-2; iy++)
 for(int iz=2; iz<nz-2; iz++){
  Cell *c = getCell(ix,iy,iz) ;
  getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz) ;
  c->getQ(Q) ;
  eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
  eta=getZ(iz) ;
  E += tau*(e+p)/(1.-vx*vx-vy*vy-tanh(vz)*tanh(vz))*(cosh(eta)-tanh(vz)*sinh(eta)) - tau*p*cosh(eta) ;
//---- inf check
  if(isinf(E)){
   cout<<"EEinf"<<setw(14)<<e<<setw(14)<<p<<setw(14)<<vx<<setw(14)<<vy<<setw(14)<<vz<<endl ;
   exit(1) ;
  }
//--------------
  Efull += tau*(e+p)/(1.-vx*vx-vy*vy-tanh(vz)*tanh(vz))*(cosh(eta)-tanh(vz)*sinh(eta)) - tau*p*cosh(eta) ;
  if(trcoeff->isViscous()) Efull += tau*c->getpi(0,0)*cosh(eta)+tau*c->getpi(0,3)*sinh(eta);
  Px += tau*(e+p)*vx/(1.-vx*vx-vy*vy-tanh(vz)*tanh(vz)) ;
  vt_num += e/sqrt(1.-vx*vx-vy*vy)*sqrt(vx*vx+vy*vy) ;
  vt_den += e/sqrt(1.-vx*vx-vy*vy) ;
  vxvy_num += e*(fabs(vx)-fabs(vy)) ;
  vxvy_den += e ;
  txxyy_num += (e+p)/(1.-vx*vx-vy*vy-tanh(vz)*tanh(vz))*(vx*vx-vy*vy) ;
  txxyy_den += (e+p)/(1.-vx*vx-vy*vy-tanh(vz)*tanh(vz))*(vx*vx+vy*vy)+2.*p ;
  pi0x_num += e/(1.-vx*vx-vy*vy-tanh(vz)*tanh(vz))*fabs(c->getpi(0,1)) ;
  pi0x_den += e/(1.-vx*vx-vy*vy-tanh(vz)*tanh(vz)) ;

  /*if(ix%compress2dOut==0&&iy%compress2dOut==0)*/
  //fout2d<<setw(10)<<tau<<setw(4)<<iz<<setw(4)<<ix<<setw(4)<<iy<<setw(14)<<e<<setw(14)<<vx<<setw(14)<<vy<<setw(14)<<vz
  //<<setw(14)<<c->getViscCorrCutFlag()<<setw(14)<<c->getpi(0,0)<<setw(14)<<c->getpi(0,1)<<setw(14)<<c->getpi(0,2)
  //<<setw(14)<<c->getpi(0,3)<<setw(14)<<c->getpi(1,1)<<setw(14)<<c->getpi(1,2)<<setw(14)<<c->getpi(1,3)<<setw(14)
  //<<c->getpi(2,2)<<setw(14)<<c->getpi(2,3)<<setw(14)<<c->getpi(3,3)<<setw(14)<<Q[0]<<setw(14)<<Q[1]<<setw(14)<<Q[2]<<setw(14)<<Q[3]
  //<<endl ;
//----- Cornelius stuff
  double vxCube[2][2][2][2], vyCube[2][2][2][2], vzCube[2][2][2][2] ;
  double TCube[2][2][2][2], mubCube[2][2][2][2], muqCube[2][2][2][2], musCube[2][2][2][2] ;
  double eCube[2][2][2][2], pCube[2][2][2][2] ;
  double piSquare[2][2][2][10], PiSquare[2][2][2] ;
  for(int jx=0; jx<2; jx++)
  for(int jy=0; jy<2; jy++)
  for(int jz=0; jz<2; jz++){
	  Cell* cc = getCell(ix+jx,iy+jy,iz+jz) ;
	  cc->getPrimVar(eos,tau,e,pCube[1][jx][jy][jz],nb,nq,ns,vxCube[1][jx][jy][jz],vyCube[1][jx][jy][jz],vzCube[1][jx][jy][jz]) ;
	  ccube[1][jx][jy][jz] = e ;
    eCube[1][jx][jy][jz] = e ;
    eos->eos(e,nb,nq,ns,TCube[1][jx][jy][jz],mubCube[1][jx][jy][jz],muqCube[1][jx][jy][jz],musCube[1][jx][jy][jz],p) ;
	  cc->getPrimVarPrev(eos,tau-dt,e,pCube[0][jx][jy][jz],nb,nq,ns,vxCube[0][jx][jy][jz],vyCube[0][jx][jy][jz],vzCube[0][jx][jy][jz]) ;
	  ccube[0][jx][jy][jz] = e ;
    eCube[0][jx][jy][jz] = e ;
    eos->eos(e,nb,nq,ns,TCube[0][jx][jy][jz],mubCube[0][jx][jy][jz],muqCube[0][jx][jy][jz],musCube[0][jx][jy][jz],p) ;
	  // ---- get viscous tensor
	  for(int ii=0; ii<4; ii++)
	  for(int jj=0; jj<=ii; jj++)
      piSquare[jx][jy][jz][index44(ii,jj)] = cc->getpi(ii,jj) ;
    PiSquare[jx][jy][jz] = cc->getPi() ;
  }
  cornelius->find_surface_4d(ccube);
  if(cornelius->get_Nelements()==1)
  for(int isegm=0; isegm<1; isegm++){
    //ffreeze<<"cell  "<<ix<<"  "<<iy<<"  "<<iz<<endl ;
    //for(int jx=0; jx<2; jx++)
    //for(int jy=0; jy<2; jy++)
    //for(int jz=0; jz<2; jz++)
    //ffreeze<<"["<<jx<<","<<jy<<","<<jz<<"] --> "<<setw(14)<<ccube[0][jx][jy][jz]<<setw(14)<<ccube[1][jx][jy][jz]<<endl ;
          ffreeze.precision(15) ;
	  ffreeze<<setw(24)<<tau+cornelius->get_centroid_elem(isegm,0)<<setw(24)<<getX(ix)+cornelius->get_centroid_elem(isegm,1)
	  <<setw(24)<<getY(iy)+cornelius->get_centroid_elem(isegm,2)<<setw(24)<<getZ(iz)+cornelius->get_centroid_elem(isegm,3) ;
	  //for(int m=0; m<4; m++) ffreeze<<setw(14)<<cornelius->get_centroid_elem(isegm,m) ;
	  //for(int m=0; m<4; m++) ffreeze<<setw(14)<<cornelius->get_normal_elem(isegm,m) ;
    //ffreeze<<endl ;
    // ---- interpolation procedure
    double vxC=0., vyC=0., vzC=0., TC=0., mubC=0., muqC=0., musC=0., piC[10], PiC=0. ; // values at the centre, to be interpolated
    double eC=0., pC=0. ;
    for(int ii=0; ii<10; ii++) piC[ii] = 0.0 ;
    double wCenT[2] = {1. - cornelius->get_centroid_elem(isegm,0)/dt, cornelius->get_centroid_elem(isegm,0)/dt} ;
    double wCenX[2] = {1. - cornelius->get_centroid_elem(isegm,1)/dx, cornelius->get_centroid_elem(isegm,1)/dx} ;
    double wCenY[2] = {1. - cornelius->get_centroid_elem(isegm,2)/dy, cornelius->get_centroid_elem(isegm,2)/dy} ;
    double wCenZ[2] = {1. - cornelius->get_centroid_elem(isegm,3)/dz, cornelius->get_centroid_elem(isegm,3)/dz} ;
    for(int jt=0; jt<2; jt++)
    for(int jx=0; jx<2; jx++)
    for(int jy=0; jy<2; jy++)
    for(int jz=0; jz<2; jz++){
      vxC += vxCube[jt][jx][jy][jz]*wCenT[jt]*wCenX[jx]*wCenY[jy]*wCenZ[jz] ;
      vyC += vyCube[jt][jx][jy][jz]*wCenT[jt]*wCenX[jx]*wCenY[jy]*wCenZ[jz] ;
      vzC += vzCube[jt][jx][jy][jz]*wCenT[jt]*wCenX[jx]*wCenY[jy]*wCenZ[jz] ;
       TC +=  TCube[jt][jx][jy][jz]*wCenT[jt]*wCenX[jx]*wCenY[jy]*wCenZ[jz] ;
     mubC +=mubCube[jt][jx][jy][jz]*wCenT[jt]*wCenX[jx]*wCenY[jy]*wCenZ[jz] ;
     muqC +=muqCube[jt][jx][jy][jz]*wCenT[jt]*wCenX[jx]*wCenY[jy]*wCenZ[jz] ;
     musC +=musCube[jt][jx][jy][jz]*wCenT[jt]*wCenX[jx]*wCenY[jy]*wCenZ[jz] ;
       eC +=  eCube[jt][jx][jy][jz]*wCenT[jt]*wCenX[jx]*wCenY[jy]*wCenZ[jz] ;
       pC +=  pCube[jt][jx][jy][jz]*wCenT[jt]*wCenX[jx]*wCenY[jy]*wCenZ[jz] ;
    }
    //double Teos ;
    //eos->eos(eC,0.,0.,0.,Teos, mub, muq, mus, p) ; // now temperature from EoS
    //if(fabs(eC-0.5)>0.01) {cout<<"+++ eC= "<<setw(14)<<eC<<", cell "<<ix<<" "<<iy<<" "<<iz<<endl;}
    for(int jx=0; jx<2; jx++)
    for(int jy=0; jy<2; jy++)
    for(int jz=0; jz<2; jz++){
      for(int ii=0; ii<10; ii++)
        piC[ii] += piSquare[jx][jy][jz][ii]*wCenX[jx]*wCenY[jy]*wCenZ[jz] ;
      PiC += PiSquare[jx][jy][jz]*wCenX[jx]*wCenY[jy]*wCenZ[jz] ;
    }
    double v2C = vxC*vxC+vyC*vyC+vzC*vzC ;
    if(v2C>1.){
      vxC *= sqrt(0.99/v2C) ;
      vyC *= sqrt(0.99/v2C) ;
      vzC *= sqrt(0.99/v2C) ;
      v2C = 0.99 ;
    }
    double etaC = getZ(iz) + cornelius->get_centroid_elem(isegm,3) ;
    transformToLab(etaC, vxC, vyC, vzC) ; // viC is now in lab.frame!
    double gammaC = 1./sqrt(1.-vxC*vxC-vyC*vyC-vzC*vzC) ;
    //ffreeze<<setw(14)<<ccube[0][0][0][0]<<setw(14)<<ccube[0][1][0][0]<<endl<<setw(14)<<ccube[1][0][0][0]<<setw(14)<<ccube[1][1][0][0]<<endl ;
    double uC [4] = {gammaC, gammaC*vxC, gammaC*vyC, gammaC*vzC} ;
    const double tauC = tau+cornelius->get_centroid_elem(isegm,0) ;
    double dsigma [4] ;
    // for(int ii=0; ii<4; ii++) dsigma[ii] = cornelius->get_normal_elem(0,ii) ;
    // ---- transform dsigma to lab.frame :
    const double ch = cosh(etaC) ;
    const double sh = sinh(etaC) ;
    dsigma [0] = tauC*( ch*cornelius->get_normal_elem(0,0) - sh/tauC*cornelius->get_normal_elem(0,3) ) ;
    dsigma [3] = tauC*(-sh*cornelius->get_normal_elem(0,0) + ch/tauC*cornelius->get_normal_elem(0,3) );
    dsigma [1] = tauC*cornelius->get_normal_elem(0,1) ;
    dsigma [2] = tauC*cornelius->get_normal_elem(0,2) ;
    double dVEff = 0.0 ;
    for(int ii=0; ii<4; ii++) dVEff += dsigma[ii]*uC[ii] ; // normalize for Delta eta=1
    vEff += dVEff ;
    for(int ii=0; ii<4; ii++) ffreeze<<setw(24)<<dsigma[ii] ;
    for(int ii=0; ii<4; ii++) ffreeze<<setw(24)<<uC[ii] ;
    ffreeze<<setw(24)<<TC<<setw(24)<<mubC<<setw(24)<<muqC<<setw(24)<<musC ;
#ifdef OUTPI
    double picart[10] ;
    /*pi00*/ picart[index44(0,0)] = ch*ch*piC[index44(0,0)]+2.*ch*sh*piC[index44(0,3)]+sh*sh*piC[index44(3,3)] ;
    /*pi01*/ picart[index44(0,1)] = ch*piC[index44(0,1)]+sh*piC[index44(3,1)] ;
    /*pi02*/ picart[index44(0,2)] = ch*piC[index44(0,2)]+sh*piC[index44(3,2)] ;
    /*pi03*/ picart[index44(0,3)] = ch*sh*(piC[index44(0,0)]+piC[index44(3,3)])+(ch*ch+sh*sh)*piC[index44(0,3)] ;
    /*pi11*/ picart[index44(1,1)] = piC[index44(1,1)] ;
    /*pi12*/ picart[index44(1,2)] = piC[index44(1,2)] ;
    /*pi13*/ picart[index44(1,3)] = sh*piC[index44(0,1)]+ch*piC[index44(3,1)] ;
    /*pi22*/ picart[index44(2,2)] = piC[index44(2,2)] ;
    /*pi23*/ picart[index44(2,3)] = sh*piC[index44(0,2)]+ch*piC[index44(3,2)] ;
    /*pi33*/ picart[index44(3,3)] = sh*sh*piC[index44(0,0)]+ch*ch*piC[index44(3,3)]+2.*sh*ch*piC[index44(0,3)] ;
    for(int ii=0; ii<10; ii++) ffreeze<<setw(24)<<picart[ii] ;
    ffreeze<<setw(24)<<PiC<<endl ;
#else
    ffreeze<<setw(24)<<dVEff<<endl ;
#endif
    double dEsurfVisc = 0. ;
    for(int i=0; i<4; i++) dEsurfVisc += picart[index44(0,i)]*dsigma[i] ;
    EtotSurf += (eC+pC)*uC[0]*dVEff - pC*dsigma[0] + dEsurfVisc ;
    /*if(fabs(eC-0.5)>0.01){ cout<<"Warning: eC = "<<setw(14)<<eC
    <<setw(14)<<getX(ix)<<setw(14)<<getY(iy)<<setw(14)<<getZ(iz)<<endl ;
    for(int i=0; i<2; i++)
    for(int j=0; j<2; j++){
      cout<<"cube: " ;
      for(int k=0; k<2; k++)
      for(int l=0; l<2; l++)
        cout<<setw(14)<<eCube[i][j][k][l] ;
      cout<<endl ;
    }
    } */ // end eC check
    //EtotSurf += (0.5+0.0795)*uC[0]*dVEff - 0.0795*dsigma[0] + dEsurfVisc ;
  // root search for EoS switching
/*  const double _p = eos->p(0.5,0.,0.,0.) ;
  const double alpha = 0.159 ;
  const double sigma2 = dsigma[0]*dsigma[0]-dsigma[1]*dsigma[1]-dsigma[2]*dsigma[2]-dsigma[3]*dsigma[3] ;
  const double AA = (0.5 + _p)*(0.5 - _p)*dVEff*dVEff + _p*_p*sigma2 ;
  const double Asigma = (0.5 + _p)*dVEff*dVEff - _p*sigma2 ;
  //ffreeze<<setw(5)<<"DDD"<<setw(24)<< AA ;
  //ffreeze<<setw(24)<< Asigma << setw(24) << 0.5*_p*dsigma2 << endl ;
  const double Det = sqrt((1.0-alpha)*(1.0-alpha)*Asigma*Asigma+4.0*alpha*sigma2*AA) ;
  ffreeze<<setw(5)<<"DDD"<<setw(24)<< (-(1.0-alpha)*Asigma - Det)/(2.0*alpha*sigma2)
    <<setw(24)<< (-(1.0-alpha)*Asigma + Det)/(2.0*alpha*sigma2) <<endl ;
*/
  }
  if(cornelius->get_Nelements()>1) cout << "oops, Nelements>1\n" ;
//----- end Cornelius
 }
 E=E*dx*dy*dz ;
 Efull=Efull*dx*dy*dz ;
 fout_aniz << setw(12) << tau << setw(14) << vt_num/vt_den <<
 setw(14) << vxvy_num/vxvy_den << setw(14) << pi0x_num/pi0x_den << endl ;
 cout << endl << setw(12) << "E = " << setw(14) << E << "  Efull = " << Efull << endl ;
 cout << setw(12) << "Px = " << setw(14) << Px << "  vEff = " << vEff << "  Esurf = " <<setw(14)<<EtotSurf << endl ;
//-- Cornelius: all done, let's free memory
  for (int i1=0; i1 < 2; i1++) {
    for (int i2=0; i2 < 2; i2++) {
      for (int i3=0; i3 < 2; i3++) {
        delete[] ccube[i1][i2][i3];
      }
      delete[] ccube[i1][i2];
    }
    delete[] ccube[i1];
  }
  delete[] ccube;
//----end Cornelius
#ifdef SWAP_EOS
  swap(eos, eosH) ;
#endif
}
