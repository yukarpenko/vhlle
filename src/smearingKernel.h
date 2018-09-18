
class GaussSmear {
 int N;  // range of smearing
 double* s;  // 2D array of the coefficients [-N..N][-N..N]
public:
 GaussSmear(int nSmear);
 ~GaussSmear();
 int getN(void) {return N; }
 double get(int i, int j) {
  if(i<-N or i>N or j<-N or j>N) return 0.;
  else return s[(i+N)*(2*N+1)+j+N];
 }
};
