#include <cmath>
#include <iostream>
#include <iomanip>
#include "smearingKernel.h"
 
using namespace std;

GaussSmear::GaussSmear(int nSmear)
{
 N = nSmear;
 s = new double [(2*N+1)*(2*N+1)];
 // special case: N=0
 if(N==0) {
  s[0] = 1.;
  return;
 }
 // calculating the smearing matrix
 double sum = 0.;
 for(int i=-N; i<=N; i++)
 for(int j=-N; j<=N; j++)
 {
  double distr = exp(- i*i/(double)(N*N) - j*j/(double)(N*N));
  s[(i+N)*(2*N+1)+j+N] = distr;
  sum += distr;
 }
 for(int i=-N; i<=N; i++)
 for(int j=-N; j<=N; j++)
 {
  s[(i+N)*(2*N+1)+j+N] /= sum;
 }
}

GaussSmear::~GaussSmear()
{
 delete [] s;
}
