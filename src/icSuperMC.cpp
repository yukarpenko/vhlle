#include <cfloat>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstring>
#include <numeric>

#include "eos.h"
#include "eoChiral.h"
#include "fld.h"
#include "icSuperMC.h"
#include "rmn.h"
#include "s95p.h"

using namespace std;

ICSuperMC::ICSuperMC(string filename, double tau0, string setup_data) : name_file(filename), tau0(tau0) {
  initialize_superMC_parameters(setup_data);
 }


ICSuperMC::~ICSuperMC(){};

void ICSuperMC::setIC(Fluid* f,  EoS *eos) {
 std::cout << "loading SuperMC\n";
 int nx = f->getNX();
 int ny = f->getNY();
 int nz = f->getNZ();

 double ybeam = acosh(sNN/(2.*mN));

 double dummyX, dummyY; //dummy variable to read the file
 int dummyId;
 double xcdmA=0,ycdmA=0,maxRA=0; //position of the center of mass and maximum distance from it for nucleus A
 double xcdmB=0,ycdmB=0,maxRB=0; //position of the center of mass and maximum distance from it for nucleus B
 int prevId=1; //used to distinguish events
 size_t size_grid = nx*ny*nz;
 vector<double> T00_xyz(size_grid,0),T0z_xyz(size_grid,0), baryon_density(size_grid,0);
 vector<double> eventxA, eventxB, eventyA, eventyB; //vector of wounded nuclei for single event 

 ifstream fin(name_file);
 if (!fin.good()) {
  std::cout << "I/O error with " << name_file << endl;
  exit(1);
 }

 int nevents = 0;  // events counter
 string line;
 istringstream instream;

 double ci = C(eta0, sigmaeta); 
 
 while (!fin.eof()) {
  getline(fin, line);
  instream.str(line);
  instream.seekg(0);
  instream.clear();
  // Read line
  instream >> dummyX >> dummyY >> dummyId;
  
  if(prevId == 2 && dummyId == 1){ //if event changes, compute quantities
    xcdmA = accumulate(eventxA.begin(),eventxA.end(),0.0)/eventxA.size();
    ycdmA = accumulate(eventyA.begin(),eventyA.end(),0.0)/eventyA.size();
    xcdmB = accumulate(eventxB.begin(),eventxB.end(),0.0)/eventxB.size();
    ycdmB = accumulate(eventyB.begin(),eventyB.end(),0.0)/eventyB.size();

    maxRA = max_distance_from_center(eventxA,eventyA,xcdmA,ycdmA);
    maxRB = max_distance_from_center(eventxB,eventyB,xcdmB,ycdmB);
    
    for (int ix = 0; ix < nx; ix++){
      for (int iy = 0; iy < ny; iy++){
        double thickA  = T(f->getX(ix),f->getY(iy),eventxA,eventyA, xcdmA, ycdmA, maxRA);
        double thickB  = T(f->getX(ix),f->getY(iy),eventxB,eventyB, xcdmB, ycdmB, maxRB);
        double y_CM = yCM(thickA, thickB, ybeam);
        double mass= Mass(thickA, thickB, ybeam, mN);
        double norm = normalization(mass, ci, sigmaeta, eta0);
  
        for(int iz = 0; iz < nz; iz++){
          double En = energy(norm, f->getZ(iz), (1-eff)*y_CM, eta0, sigmaeta,ybeam);
          T00_xyz[iz + nz*iy + nz*ny*ix]+= En*cosh(eff*y_CM);
          T0z_xyz[iz + nz*iy + nz*ny*ix]+= En*sinh(eff*y_CM); //not needed to divide by tau0 because Q[Z_] already takes that into account      

          baryon_density[iz + nz*iy + nz*ny*ix] += (thickA*baryon_density_profile(f->getZ(iz), etaB, sigmaIN, sigmaOUT, 1) +
                                              thickB*baryon_density_profile(f->getZ(iz), etaB, sigmaIN, sigmaOUT, 2)); //taken from 2106.08125
        }
      }
    }
    
    eventxA.clear();
    eventxB.clear();
    eventyA.clear();
    eventyB.clear();
    nevents++;
    xcdmA = 0;
    ycdmA = 0;
    maxRA = 0;
    xcdmB = 0;
    ycdmB = 0;
    maxRB = 0;
  }

  if(dummyId == 1){ //check nucleon id to set variables (single event)
    eventxA.push_back(dummyX);
    eventyA.push_back(dummyY);
  }
  if(dummyId == 2){
    eventxB.push_back(dummyX);
    eventyB.push_back(dummyY);
  }
  prevId = dummyId;                                              
 }

  xcdmA = accumulate(eventxA.begin(),eventxA.end(),0.0)/eventxA.size();
  ycdmA = accumulate(eventyA.begin(),eventyA.end(),0.0)/eventyA.size();
  xcdmB = accumulate(eventxB.begin(),eventxB.end(),0.0)/eventxB.size();
  ycdmB = accumulate(eventyB.begin(),eventyB.end(),0.0)/eventyB.size();

  maxRA = max_distance_from_center(eventxA,eventyA,xcdmA,ycdmA);
  maxRB = max_distance_from_center(eventxB,eventyB,xcdmB,ycdmB);
  
  for (int ix = 0; ix < nx; ix++){
    for (int iy = 0; iy < ny; iy++){
      double thickA  = T(f->getX(ix),f->getY(iy),eventxA,eventyA, xcdmA, ycdmA, maxRA);
      double thickB  = T(f->getX(ix),f->getY(iy),eventxB,eventyB, xcdmB, ycdmB, maxRB);
      double y_CM = yCM(thickA, thickB, ybeam);
      double mass= Mass(thickA, thickB, ybeam, mN);
      double norm = normalization(mass, ci, sigmaeta, eta0);

      for(int iz = 0; iz < nz; iz++){
        double En = energy(norm, f->getZ(iz), (1-eff)*y_CM, eta0, sigmaeta,ybeam);
        T00_xyz[iz + nz*iy + nz*ny*ix]+= En*cosh(eff*y_CM);
        T0z_xyz[iz + nz*iy + nz*ny*ix]+= En*sinh(eff*y_CM); //not needed to divide by tau0 because Q[Z_] already takes that into account       

        baryon_density[iz + nz*iy + nz*ny*ix] += (thickA*baryon_density_profile(f->getZ(iz), etaB, sigmaIN, sigmaOUT, 1) +
                                            thickB*baryon_density_profile(f->getZ(iz), etaB, sigmaIN, sigmaOUT, 2));
      }
    }
  }
  
  eventxA.clear();
  eventxB.clear();
  eventyA.clear();
  eventyB.clear();
  nevents++;
  xcdmA = 0;
  ycdmA = 0;
  maxRA = 0;
  xcdmB = 0;
  ycdmB = 0;
  maxRB = 0;

  fin.close(); //done reading


double e, p, nb, nq, ns, vx, vy, vz;
double EMAX = 0;
double Q[7]={0,0,0,0,0,0,0};
double gamma = 0;

for (int ix = 0; ix < nx; ix++){
  for (int iy = 0; iy < ny; iy++){
   for (int iz = 0; iz < nz; iz++) {        
    Q[T_] = T00_xyz[iz + nz*iy + nz*ny*ix]/nevents;
    Q[X_] = 0.;
    Q[Y_] = 0.;
    Q[Z_] = T0z_xyz[iz + nz*iy + nz*ny*ix]/nevents;
    
    transformPV(eos, Q, e, p, nb, nq, ns, vx, vy, vz); //we need the gamma factor to initialize the currents
    gamma = 1./sqrt(1-vx*vx-vy*vy-vz*vz);
    
    Q[NB_] = gamma*baryon_density[iz + nz*iy + nz*ny*ix]/nevents;
    Q[NQ_] = gamma*ZoverA*baryon_density[iz + nz*iy + nz*ny*ix]/nevents;
    Q[NS_] = 0.;
    
    transformPV(eos, Q, e, p, nb, nq, ns, vx, vy, vz);
    if(e>EMAX) EMAX = e;
    // if (e < 1e-7 || fabs(f->getX(ix)) > 10. || fabs(f->getY(iy)) > 10. ||
    //     fabs(f->getZ(iz)) > 5.) {
    //  e = nb = nq = 0.0;
    //  vx = vy = vz = 0.0;
    // }
    Cell* c = f->getCell(ix, iy, iz);
    c->setPrimVar(eos, tau0, e, nb, nq, ns, vx, vy, vz);
    if (e > 1e-7) c->setAllM(1.);
   }
  }
}
cout<<"The max Energy density in the LRF (in fm^-4) is: "<<EMAX<<endl;
}

double ICSuperMC::T(double x, double y, vector<double> &xa, vector<double> &ya, double xcdm, double ycdm, double maxR){
    double T=0;
    if(sqrt(pow((x-xcdm),2)+pow((y-ycdm),2)) < maxR +4*w){ //thickness is set to zero if the point of the grid is far from the distribution on participants
      for(size_t partA=0; partA<xa.size();partA++){
          T +=1./(2.*M_PI*w*w)*exp(-1./(2.*w*w)*(pow(x-xa[partA],2)+pow(y-ya[partA],2)));
      }
    }
    return T;
}

double ICSuperMC::yCM(double TA, double TB, double ybeam){
    double ycm = atanh((TA-TB)/(TA+TB+1e-8) *tanh(ybeam));
    return ycm;
}

double ICSuperMC::Mass(double TA, double TB, double ybeam, double mN){
    double M = mN*sqrt(TA*TA+TB*TB+2.*TA*TB*cosh(2.*ybeam));
    //double gevtofm = 5.067728853; //this way the mass is in 1/fm^3 units
    return M;
}

double ICSuperMC::normalization(double M, double C, double sigmaeta, double eta0){
    double N = M/(2.*sinh(eta0) + sqrt(0.5*M_PI) *sigmaeta *exp(sigmaeta*sigmaeta*0.5) *C );
    return N;
}

double ICSuperMC::C(double eta0, double sigmaeta){
    double c = exp(eta0) *erfc(-sqrt(0.5)*sigmaeta) + exp(-eta0) *erfc(sqrt(0.5)*sigmaeta);
    return c;
}

double ICSuperMC::energy(double N, double etas, double ycombo, double eta0, double sigmaeta, double ybeam){
    double eta0_new = std::min(eta0,ybeam - ycombo);
    double absarg = abs(etas-ycombo)-eta0_new; //ycombo = ycm-yl=(1-f)ycm
    if(absarg>=0){
        return N*exp(-1./(2*sigmaeta*sigmaeta) *absarg*absarg);
    }
    return N;
}

double ICSuperMC::baryon_density_profile(double etas, double etaB, double sigmaIN, double sigmaOUT, int id){
  // refer to 2003.05852 for the notation
    double result = 0;
    double norm = 1./(sqrt(M_PI_2)*(sigmaIN+sigmaOUT)); //normalization computed explicitly from 1804.10557
    if(id==1){
      double variable = etas -etaB;
      if(variable>=0){
          result = norm*exp(-1./(2*sigmaOUT*sigmaOUT) *variable*variable);
      }
      else{
        result = norm*exp(-1./(2*sigmaIN*sigmaIN) *variable*variable);
      }
    }
    else if(id==2){
      double variable = etas + etaB;
      if(variable>=0){
        result = norm*exp(-1./(2*sigmaIN*sigmaIN) *variable*variable);
      }    
      else{
        result = norm*exp(-1./(2*sigmaOUT*sigmaOUT) *variable*variable);
      }  
    }
    else{
      cout<<"problem with particle ID!"<<endl;
      exit(1);
      return 0;
    }
    return result;
}

double ICSuperMC::max_distance_from_center(vector<double> x, vector<double> y, double xcdm, double ycdm){
  double max_dist = 0;
  for(int i=0; i<x.size();i++){
      if( max_dist<sqrt(pow((x[i]-xcdm),2) + pow((y[i]-ycdm),2)) ){
        max_dist = sqrt(pow((x[i]-xcdm),2) + pow((y[i]-ycdm),2));
      }
    }
  return max_dist;
}

void ICSuperMC::initialize_superMC_parameters(std::string setup_data){
  ifstream fin(setup_data);
  if (!fin.is_open()) {
    cout << "cannot open superMC file " << setup_data << endl;
    exit(1);
  }
  cout << "Reading superMC parameters from " << setup_data << endl;
  char parName[255], parValue[255];
  while (fin.good()) {
    string line;
    getline(fin, line);
    istringstream sline(line);
    sline >> parName >> parValue;
    if (strcmp(parName, "sNN") == 0)
    sNN = atof(parValue);
    else if (strcmp(parName, "eta0") == 0)
    eta0 = atof(parValue);
    else if (strcmp(parName, "sigmaeta") == 0)
    sigmaeta = atof(parValue);
    else if (strcmp(parName, "w") == 0)
    w = atof(parValue);
    else if (strcmp(parName, "eff") == 0)
    eff = atof(parValue);
    else if (strcmp(parName, "etaB") == 0)
    etaB = atof(parValue);
    else if (strcmp(parName, "sigmaIN") == 0)
    sigmaIN = atof(parValue);
    else if (strcmp(parName, "sigmaOUT") == 0)
    sigmaOUT = atof(parValue);  
    }

    cout<< "superMC parameters initialized:"<<endl;
    cout<<"sNN: " << sNN << endl;
    cout<<"eta0: " << eta0 << endl;
    cout<<"sigmaeta: " << sigmaeta << endl;
    cout<<"w: "<< w << endl;
    cout<<"eff: " << eff << endl;
    cout<<"etaB: " << etaB << endl;
    cout<<"sigmaIN: " << sigmaIN << endl;
    cout<<"sigmaOUT: " << sigmaOUT<<endl;
  }
  


