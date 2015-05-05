#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "jets.h"

using namespace std;

// --------- individual jet ---------

void Jet::propagate(double tau, double dtau)
{
 double pt = sqrt(px*px+py*py);
 double boost = cosh(rap_jet)/cosh(rap_jet-eta);
 x += px / pt * vt_jet * boost * dtau;
 y += py / pt * vt_jet * boost * dtau;
 eta += tanh(rap_jet - eta) / tau * dtau;
}


void Jet::addEnergyMom(double dE, double dPx, double dPy, double dPz)
{
 double E  = mt*cosh(rap);
 double Pz = mt*sinh(rap);
 E  += dE;
 px += dPx;
 py += dPy;
 Pz += dPz;
 if(E>0. && fabs(Pz)<E){
 rap = atanh(Pz/E);
 mt = sqrt(E*E - Pz*Pz);
 }else{
  mt = -1.0;
  cout<<"Warning: jet lost all its energy\n";
 }
}

// --------- Jets -------------------

Jets::Jets(char* filename)
{
 ifstream fin(filename);
 double tau0;
 if(!fin){ cout<<"cannot open "<<filename<<endl; exit(1); }
 while (!fin.eof()) {
  double x, y, eta, px, py, rap, mt, vt_jet, rap_jet;
  fin >> tau0 >> x >> y >> eta >> mt >> px >> py >> rap >> vt_jet >> rap_jet;
  if(!fin.fail())
   jets.push_back(Jet(x, y, eta, px, py, rap, mt, vt_jet, rap_jet));
 }
 cout<<"Jets:  tau0 = "<<tau0<<endl;
 cout<<"Jets:  "<<jets.size()<<"  jets read in.\n";
}


void Jets::propagate(double tau, double dtau)
{
 for(vector<Jet>::iterator it=jets.begin(); it!=jets.end(); it++){
  it->propagate(tau, dtau);
 }
}


void Jets::checkJets(void)
{
  for(vector<Jet>::iterator it=jets.begin(); it!=jets.end(); it++){
  if(it->Mt()<0.){
   cout<<"Jets: removing jet\n";
   jets.erase(it);
  }
 }
}
