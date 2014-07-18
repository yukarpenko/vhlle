#include "eos.h"
#include "inc.h"

void allshock(double el, double vl, double nl, double er, double vr, double nr, EoS *eos,
			  double &e1, double &v1, double &n1, double &e2, double &n2, double &s1, double &s3) ;

void transformPV(EoS *eos, double Q[7], double &e, double &p, double &nb, double &nq, double &ns, double &vx, double &vy, double &vz) ;
void transformPVBulk(EoS *eos, double Pi, double Q[7], double &e, double &p, double &nb, double &nq, double &ns, double &vx, double &vy, double &vz) ;

void transformCV(double e, double p, double nb, double nq, double ns, double vx, double vy, double vz, double Q []) ;
