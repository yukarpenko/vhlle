//#include <TError.h>
//#include <TApplication.h>
//#include <TGraph.h>
//#include <TCanvas.h>
//#include <TMath.h>
//#include <TGraph.h>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <cstdlib>

#include <iostream>
#include <fstream>
//#include <TF1.h>

#include "eos.h"
#include "eoHadron.h"

using namespace std ;


EoSHadron::EoSHadron(char *filename)
{
  ifstream fin(filename) ;
  if(!fin.good()){ cout<<"I/O error with "<<filename<<endl ;  exit(1) ; }
  double eminhead, emaxhead, nbminhead, nbmaxhead, nqminhead, nqmaxhead ; // these are read from the header line
  fin >> ne >> nnb >> nnq >> eminhead >> emaxhead >> nbminhead >> nbmaxhead >> nqminhead >> nqmaxhead ;

  ptab = new double [ne*nnb*nnq] ;
  Ttab = new double [ne*nnb*nnq] ;
mubtab = new double [ne*nnb*nnq] ;
muqtab = new double [ne*nnb*nnq] ;
mustab = new double [ne*nnb*nnq] ;
statustab = new int [ne*nnb*nnq] ;

  double* e = new double [ne] ;
  double* nb = new double [nnb] ;
  double* nq = new double [nnq] ;


  for(int ie=0; ie<ne; ie++)
  for(int inb=0; inb<nnb; inb++)
  for(int inq=0; inq<nnq; inq++){
    // e  nb  p  T  mub  muq  mus  status
    fin >> e[ie] >> nb[inb] >> nq[inq] >> ptab[index3(ie,inb,inq)] >> Ttab[index3(ie,inb,inq)]
     >> mubtab[index3(ie,inb,inq)] >> muqtab[index3(ie,inb,inq)] >> mustab[index3(ie,inb,inq)] 
     >> statustab[index3(ie,inb,inq)] ;
  }
  emin = e[0] ;
  emax = e[ne-1] ;
  nbmin = nb[0] ;
  nbmax = nb[nnb-1] ;
  nqmin = nq[0] ;
  nqmax = nq[nnq-1] ;
  if(fabs(eminhead-emin)>1e-4 || fabs(emaxhead-emax)>1e-4 || fabs(nbminhead-nbmin)>1e-4 || fabs(nbmaxhead-nbmax)>1e-4)
    { cout << "wrong eps or nb range: "<<setw(14)<<emin<<setw(14)<<emax<<setw(14)<<nbmin<<setw(14)<<nbmax<<endl ; exit(1) ; }
  cout<<"EoHadron: table "<<filename<<" read, [emin,emax,nmin,nmax] = "<<emin<<"  "<<emax<<"  "<<nbmin<<"  "<<nbmax
  <<"  "<<nqmin<<"  "<<nqmax<<endl ;
  delete [] e ;
  delete [] nb ;
  delete [] nq ;
}


EoSHadron::~EoSHadron()
{
    delete [] ptab ;
    delete [] Ttab ;
    delete [] mubtab ;
    delete [] muqtab ;
    delete [] mustab ;
    delete [] statustab ;
}


void EoSHadron::eos(double e, double nb, double nq, double ns,
		double &T, double &mub, double &muq, double &mus, double &p)
{
	if(e<0.) { T = mub = muq = mus = p = 0. ; return ; }
  const double de = (emax-emin)/(ne-1) ;
  const double dnb = (nbmax-nbmin)/(nnb-1) ;
  const double dnq = (nqmax-nqmin)/(nnq-1) ;
  int ie = (int)((e-emin)/de) ;
  int inb = (int)((nb-nbmin)/dnb);
  int inq = (int)((nq-nqmin)/dnq);
  if(ie<0) ie=0 ;
  if(inb<0) inb=0 ;
  if(inq<0) inq=0 ;
  if(ie>ne-2) ie=ne-2 ;
  if(inb>nnb-2) inb=nnb-2 ;
  if(inq>nnq-2) inq=nnq-2 ;
  const double em =  e-emin-ie*de ;
  const double nbm = nb-nbmin-inb*dnb ;
  const double nqm = nq-nqmin-inq*dnq ;
  
  if(statustab[index3(ie,inb+1,inq+1)]==1) { T = mub = muq = mus = p = 0. ; return ; }
                
	double we [2] = {1.-em/de, em/de} ;
	double wnb [2] = {1.-nbm/dnb, nbm/dnb} ;
  double wnq [2] = {1.-nqm/dnq, nqm/dnq} ;

	T = mub = muq = mus = p = 0.0 ;
	for(int je=0; je<2; je++)
	for(int jnb=0; jnb<2; jnb++)
  for(int jnq=0; jnq<2; jnq++){
		p += we[je]*wnb[jnb]*wnq[jnq]*ptab[index3(ie+je,inb+jnb,inq+jnq)] ;
		T += we[je]*wnb[jnb]*wnq[jnq]*Ttab[index3(ie+je,inb+jnb,inq+jnq)] ;
		mub += we[je]*wnb[jnb]*wnq[jnq]*mubtab[index3(ie+je,inb+jnb,inq+jnq)] ;
    muq += we[je]*wnb[jnb]*wnq[jnq]*muqtab[index3(ie+je,inb+jnb,inq+jnq)] ;
    mus += we[je]*wnb[jnb]*wnq[jnq]*mustab[index3(ie+je,inb+jnb,inq+jnq)] ;
	}
	if(p<0.0) p = 0.0 ;
        // cout <<  e <<" "<< nb <<" "<< nq <<" "<< ns <<" "<< _T<<endl;
}


double EoSHadron::p(double e, double nb, double nq, double ns)
{
	if(e<0.) return 0.0 ;
  const double de = (emax-emin)/(ne-1) ;
  const double dnb = (nbmax-nbmin)/(nnb-1) ;
  const double dnq = (nqmax-nqmin)/(nnq-1) ;
  int ie = (int)((e-emin)/de) ;
  int inb = (int)((nb-nbmin)/dnb);
  int inq = (int)((nq-nqmin)/dnq);
  if(ie<0) ie=0 ;
  if(inb<0) inb=0 ;
  if(inq<0) inq=0 ;
  if(ie>ne-2) ie=ne-2 ;
  if(inb>nnb-2) inb=nnb-2 ;
  if(inq>nnq-2) inq=nnq-2 ;
  const double em =  e-emin-ie*de ;
  const double nbm = nb-nbmin-inb*dnb ;
  const double nqm = nq-nqmin-inq*dnq ;
  
  if(statustab[index3(ie,inb+1,inq+1)]==1)  return 0.0 ;
                
	double we [2] = {1.-em/de, em/de} ;
	double wnb [2] = {1.-nbm/dnb, nbm/dnb} ;
  double wnq [2] = {1.-nqm/dnq, nqm/dnq} ;

	double p = 0.0 ;
	for(int je=0; je<2; je++)
	for(int jnb=0; jnb<2; jnb++)
  for(int jnq=0; jnq<2; jnq++){
		p += we[je]*wnb[jnb]*wnq[jnq]*ptab[index3(ie+je,inb+jnb,inq+jnq)] ;
	}
	if(p<0.0) p = 0.0 ;
  return p ;
}
