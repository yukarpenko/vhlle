#pragma once
#include <iosfwd>
#include <algorithm>
#include "inc.h"
class EoS ;

int index44(const int &i, const int &j) ;

class Cell
{
private :
	double Q[7] ;  // here, Q, Qh, Qprev etc. ~tau*T^{0i}, like Hirano'01
	double Qh[7] ; // values at (n+1/2) timestep
	double Qprev [7] ; // values at the end of previous timestep
	double Qfull [7] ; // full T^{0\mu} with viscous terms, WITHOUT tau factor
  double pi[10], piH[10], Pi, PiH ;     //viscous, WITHOUT tau factor
  double pi0[10], piH0[10], Pi0, PiH0 ; //viscous, WITHOUT tau factor
	double flux[7] ;
	Cell *next [3] ;
	Cell *prev [3] ;
	double m [3] ;
	double dm [3] ;
	int ix, iy, iz ;
  double viscCorrCut ; // flag if the viscous corrections are cut for this cell
public :
	Cell() ;
	~Cell() {} ;
	inline void setPos(int iix, int iiy, int iiz) { ix = iix; iy = iiy; iz = iiz ; }
	inline int getX(void) { return ix ; }
	inline int getY(void) { return iy ; }
	inline int getZ(void) { return iz ; }

	inline void setQ(double *_Q){
		for(int i=0; i<7; i++) Q[i] = _Q[i] ;
		if(Q[T_]<0.){ for(int i=0; i<7; i++) Q[i]=0. ; }
	}
	inline void setQh(double *_Qh){
		for(int i=0; i<7; i++) Qh[i] = _Qh[i] ;
		if(Qh[T_]<0.){ for(int i=0; i<7; i++) Qh[i]=0. ; }
	}
	inline void setQfull(double *_Qf){
		for(int i=0; i<7; i++) Qfull[i] = _Qf[i] ;
		if(Qfull[T_]<0.){ for(int i=0; i<7; i++) Qfull[i]=0. ; 
		/*cout<<"Cell: Qfull[0]<0 "; for(int j=0;j<4;j++)cout<<"  "<<Qfull[j] ;
		cout<<endl ;*/ }
	}

  inline double getpi(const int &i, const int &j){ return pi[index44(i,j)] ; }
  inline double getpiH(const int &i, const int &j){ return piH[index44(i,j)] ; }
  inline double getpi0(const int &i, const int &j){ return pi0[index44(i,j)] ; }
  inline double getpiH0(const int &i, const int &j){ return piH0[index44(i,j)] ; }
  inline double getPi(void){ return Pi ; }
  inline double getPiH(void){ return PiH ; }
  inline double getPi0(void){ return Pi0 ; }
  inline double getPiH0(void){ return PiH0 ; }
 
  inline void setpi(const int &i, const int &j, const double &val){ pi[index44(i,j)] = val ; }
  inline void setpiH(const int &i, const int &j, const double &val){ piH[index44(i,j)] = val ; }
  inline void setpi0(const int &i, const int &j, const double &val){ pi0[index44(i,j)] = val ; }
  inline void setpiH0(const int &i, const int &j, const double &val){ piH0[index44(i,j)] = val ; }
  inline void addpi0(const int &i, const int &j, const double &val){ pi0[index44(i,j)] += val ; }
  inline void addpiH0(const int &i, const int &j, const double &val){ piH0[index44(i,j)] += val ; }
  inline void setPi(const double &val){ Pi = val ; }
  inline void setPiH(const double &val){ PiH = val ; }
  inline void setPi0(const double &val){ Pi0 = val ; }
  inline void setPiH0(const double &val){ PiH0 = val ; }
  inline void addPi0(const double &val){ Pi0 += val ; }
  inline void addPiH0(const double &val){ PiH0 += val ; }

	inline void getQ(double *_Q) { for(int i=0; i<7; i++) _Q[i] = Q[i] ; }
	inline void getQh(double *_Qh) { for(int i=0; i<7; i++) _Qh[i] = Qh[i] ; }
  inline void getQprev(double *_Qp) { for(int i=0; i<7; i++) _Qp[i] = Qprev[i] ; }
	inline void getQfull(double *_Qf) { for(int i=0; i<7; i++) _Qf[i] = Qfull[i] ; }
	inline void saveQprev(void) { for(int i=0; i<7; i++) Qprev[i] = Q[i] ; }
	inline void setNext(int i, Cell* c) { next[i-1] = c ; }
	inline void setPrev(int i, Cell* c) { prev[i-1] = c ; }
	inline Cell* getNext(int i) { return next[i-1] ; }
	inline Cell* getPrev(int i) { return prev[i-1] ; }

	inline void setAllM(double value) { m[0] = m[1] = m[2] = value ; }
	inline void addM(int dir, double inc) { m[dir-1] += inc; if(m[dir-1]>0.9) for(int i=0; i<3; i++) m[i]=1. ; }
	inline double getM(int dir) { return m[dir-1] ; }
	inline double getLM(void) { return std::max(m[0],std::max(m[1],m[2])) ; }
	inline void setDM(int dir, double value) { dm[dir-1] = value ; }
	inline double getDM(int dir) { return dm[dir-1] ; }
	
	inline void setpi0(double values[4][4])
	{ for(int i=0; i<4; i++) for(int j=0; j<4; j++) pi0[index44(i,j)] = values[i][j] ; }
	inline void setpiH0(double values[4][4])
	{ for(int i=0; i<4; i++) for(int j=0; j<4; j++) piH0[index44(i,j)] = values[i][j] ; }

	void getPrimVar(EoS *eos, double tau, double &_e, double &_p, double &_nb, double &_nq, double &_ns, double &_vx, double &_vy, double &_vz) ;

	void getPrimVarLeft(EoS *eos, double tau, double &_e, double &_p, double &_nb, double &_nq, double &_ns, double &_vx, double &_vy, double &_vz, int dir) ;
	void getPrimVarRight(EoS *eos, double tau, double &_e, double &_p, double &_nb, double &_nq, double &_ns, double &_vx, double &_vy, double &_vz, int dir) ;

	void getPrimVarHLeft(EoS *eos, double tau, double &_e, double &_p, double &_nb, double &_nq, double &_ns, double &_vx, double &_vy, double &_vz, int dir) ;
	void getPrimVarHRight(EoS *eos, double tau, double &_e, double &_p, double &_nb, double &_nq, double &_ns, double &_vx, double &_vy, double &_vz, int dir) ;
	void getPrimVarHCenter(EoS *eos, double tau, double &_e, double &_p, double &_nb, double &_nq, double &_ns, double &_vx, double &_vy, double &_vz) ;
	void getPrimVarPrev(EoS *eos, double tau, double &_e, double &_p, double &_nb, double &_nq, double &_ns, double &_vx, double &_vy, double &_vz) ;
	void getPrimVarFull(EoS *eos, double &_e, double &_p, double &_nb, double &_nq, double &_ns, double &_vx, double &_vy, double &_vz) ;
	void setPrimVar(EoS *eos, double tau, double _e, double _nb, double _nq, double _ns, double _vx, double _vy, double _vz) ;

	inline void addFlux(double Ft, double Fx, double Fy, double Fz, double Fnb, double Fnq, double Fns) 
	{ flux[T_] += Ft ; flux[X_] += Fx ; flux[Y_] += Fy ; flux[Z_] += Fz ; flux[NB_] += Fnb ; flux[NQ_] += Fnq ; flux[NS_] += Fns ; }
	inline void clearFlux(void) { for(int i=0; i<7; i++) flux[i] = 0. ; }
	void updateByFlux() ;
	void updateQtoQhByFlux() ;
	void updateQfullByFlux() ;
	void correctQideal(EoS *eos, double tau) ;
  inline void setViscCorrCutFlag(double value) { viscCorrCut=value ; }
  inline double getViscCorrCutFlag(void) { return viscCorrCut ; }
	void Dump(double tau) ;
};
