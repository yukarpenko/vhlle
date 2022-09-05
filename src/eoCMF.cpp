#include <math.h>
#include <iomanip>
#include <ctime>
#include <cstdlib>

#include <iostream>
#include <fstream>

#include "eos.h"
#include "eoCMF.h"

using namespace std;

double EPS_CUTOFF = 1e-5;

EoSCMF::EoSCMF() : EoSCMF("eos/CMF_eos.dat", 1001, 401){
};

EoSCMF::EoSCMF(std::string filename, int N_T, int N_nB) {

 this->N_T = N_T;
 this->N_nB = N_nB;
 ptab = new double*[N_T];
 etab = new double*[N_T];
 mubtab = new double*[N_T];
 mustab = new double*[N_T];
 stab = new double*[N_T];
 qfractab = new double*[N_T];
 for (int i = 0; i < N_T; i++) {
  ptab[i] = new double[N_nB];
  etab[i] = new double[N_nB];
  mubtab[i] = new double[N_nB];
  mustab[i] = new double[N_nB];
  stab[i] = new double[N_nB];
  qfractab[i] = new double[N_nB];
 }
 double* temp = new double[N_T];
 double* n = new double[N_nB];

 ifstream eosfile(filename);
 if (!eosfile.good()) {
  cout << "I/O error with " << filename << endl;
  exit(1);
 }

 dT = 0.5;
 dN = 0.1;
 EPS_0 = 146.51751415742*0.87326281;
 N_0 = 0.15891*0.87272727;

 int i_nB, i_T;
 double nB, T;
 nmax = 0;
 for(int i=0;i<N_T*N_nB;i++) {
  eosfile >> T >>
          nB;
  i_nB =(nB+0.5*dN)/dN;
  i_T = (T+0.5*dT)/dT;
  n[i_nB] = nB;
  temp[i_T] = T;
  eosfile >> etab[i_T][i_nB] >>
          ptab[i_T][i_nB] >>
          stab[i_T][i_nB] >>
          mubtab[i_T][i_nB] >>
          mustab[i_T][i_nB] >>
          qfractab[i_T][i_nB];

  etab[i_T][i_nB] *= EPS_0/1000.0;    // --> e[GeV/fm3]
  mubtab[i_T][i_nB] /= 1000.0;  // --> mub[GeV]
  mustab[i_T][i_nB] /= 1000.0;  // --> mus[GeV]
  ptab[i_T][i_nB] *= EPS_0/1000.0;     // --> p[GeV/fm3]
  stab[i_T][i_nB] *= N_0;      // --> s[1/fm3]
 }


 Tmin = temp[0]/1000.;
 Tmax = temp[N_T - 2]/1000.;
 nmin = n[0] * N_0;
 nmax = n[N_nB - 1] * N_0;
 cout << "EoSCMF: table " << filename
     << " read, [Tmin,Tmax,nmin,nmax] = " << Tmin << "  " << Tmax << "  "
     << nmin << "  " << nmax << endl;
 delete[] n;
 delete[] temp;
}

EoSCMF::~EoSCMF() {
 for (int i = 0; i < N_T; i++) {
  delete[] ptab[i];
  delete[] etab[i];
  delete[] mubtab[i];
  delete[] mustab[i];
  delete[] stab[i];
 }
 delete ptab;
 delete etab;
 delete mubtab;
 delete mustab;
 delete stab;
}

void EoSCMF::get(double e, double nb, double& p, double& T, double& mub,
                 double& mus) {

 if (e < EPS_CUTOFF) {
  T = mub = mus = p = 0.;
  return;
 }

 T = EoSCMF::get_temp(e, nb);
 if(T < 10.){
   T = mub = mus = p = 0.;
   return;
 }

 int i_T = (int)((T - Tmin + dT*0.5) / dT);
 int i_nB = (int)((nb - nmin + dN*0.5) / dN);
 if (i_T < 0) i_T = 0;
 if (i_nB < 0) i_nB = 0;
 if (i_T > i_T - 2) i_T = i_T - 2;
 if (i_nB > i_nB - 2) i_nB = i_nB - 2;
 const double Tm = T - Tmin - i_T * dT;
 const double nm = nb - nmin - i_nB * dN;

 double we[2] = {1. - Tm / dT, Tm / dT};
 double wn[2] = {1. - nm / dN, nm / dN};

 mub = mus = p = 0.0;
 for (int je = 0; je < 2; je++)
  for (int jn = 0; jn < 2; jn++) {
   p += we[je] * wn[jn] * ptab[i_T + je][i_nB + jn];
   mub += we[je] * wn[jn] * mubtab[i_T + je][i_nB + jn];
   mus += we[je] * wn[jn] * mustab[i_T + je][i_nB + jn];
  }
 if (p < 0.0) p = 0.0;
}

double EoSCMF::p(double e, double nb, double nq, double ns) {
// if (e < 146. && nb < 6.) {
   if (e < 1.63 && nb < 6.) {
     if (e < EPS_CUTOFF) return 0.0;

     double T = EoSCMF::get_temp(e, nb);

     int i_T = (int)((T - Tmin + dT*0.5) / dT);
     int i_nB = (int)((nb - nmin+ dN*0.5) / dN);
     if (i_T < 0) i_T = 0;
     if (i_nB < 0) i_nB = 0;
     if (i_T > N_T - 2) i_T = i_T - 2;
     if (i_nB > N_nB - 2) i_nB = i_nB - 2;
     const double Tm = T - Tmin - i_T * dT;
     const double nm = nb - nmin - i_nB * dN;

     double we[2] = {1. - Tm / dT, Tm / dT};
     double wn[2] = {1. - nm / dN, nm / dN};

     double p = 0.0;
     for (int je = 0; je < 2; je++)
         for (int jn = 0; jn < 2; jn++) p += we[je] * wn[jn] * ptab[i_T + je][i_nB + jn];

     if (p < 0.0) p = 0.0;
     return p;
 }else
     return 0.2964 * e;
}

void EoSCMF::eos(double e, double nb, double nq, double ns, double& T,
                    double& mub, double& muq, double& mus, double& p) {
 if (e < EPS_CUTOFF) {
   p = 0.;
   T = 0.;
   mub = mus = 0.0;
   muq = 0.0;
   return;
 }

 if (EPS_CUTOFF < e && e < 146. &&  nb < 6.)
     this->get(e, nb, p, T, mub, mus);
 else {
  p = 0.2964 * e;
  T = 0.15120476935 * pow(e, 0.25);
  mub = mus = 0.0;
 }
 muq = 0.0;  // generally it's not zero
}

int EoSCMF::binary_search(double *a, int first, int last, double search_num) {
//        binary search to find index of closest element of array a[]
//        to the  values of searh_num
//        the function is executed iteratively
 int middle;

 if (last > first) {
  if(abs(last-first) == 1){
   if(fabs(a[last]-search_num) < fabs(a[first]-search_num))
    return last;
   else
    return first;
  }
  middle = (first + last) / 2;
  //Checking if the element is present at middle loc
  if (a[middle] == search_num) {
   return middle;
  }

   //Checking if the search element is present in greater half
  else{ if (a[middle] < search_num) {
    return binary_search(a, middle + 1, last, search_num);
   }
    //Checking if the search element is present in lower half
   else{
    return binary_search(a, first, middle - 1, search_num);
   }
  }

 }
 else if (last == first)
  return last;
 return -1;
}

int EoSCMF::get_ind(double eps, double n) {

 int nB_ind = (n + 0.5*dN) / dN;
 if (nB_ind < 0) nB_ind = 0;

 double a[N_T];

 for (int i = 0; i < N_T; i++) {
  a[i] = etab[i][nB_ind];
 }
 
 int t_ind = binary_search(a, 0, N_T, eps);
 if (t_ind >= N_T-1) t_ind = N_T-2;

 return t_ind;
}

double EoSCMF::get_temp(double eps, double n) {
 double el1, el2, eh1, eh2, elow, ehigh, etabval, ediff1, ediff2, temp;
 int n_ind = n / dN;
 if (n_ind < 0) n_ind = 0;

 int t_ind = get_ind(eps, n);

 el1 = etab[t_ind][n_ind];
 el2 = etab[t_ind][n_ind + 1];
 eh1 = etab[t_ind + 1][n_ind];
 eh2 = etab[t_ind + 1][n_ind + 1];

 elow = (el2 - el1) / dN * (n - (double) n_ind * dN) + el1;
 ehigh = (eh2 - eh1) / dN * (n - (double) n_ind * dN) + eh1;

 etabval = ehigh - elow;
 ediff1 = eps - ehigh;
 ediff2 = eps - elow;
//    if (ediff1 *ediff2 < 0.0f) {
 temp = ((double) t_ind + (ediff2 / etabval))*dT;

 if(temp > 500.){
   temp = 499.;
 }

 return temp;
}
