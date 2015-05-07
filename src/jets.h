#include <vector>

// class for individual jet
class Jet {
  double x, y, eta;
  double px, py, rap, mt, vt_jet, rap_jet;

 public:
  Jet(double _x, double _y, double _eta, double _px, double _py, double _rap,
      double _mt, double _vt_jet, double _rap_jet)
      : x(_x),
        y(_y),
        eta(_eta),
        px(_px),
        py(_py),
        rap(_rap),
        mt(_mt),
        vt_jet(_vt_jet),
        rap_jet(_rap_jet){};
  double X(void) { return x; }
  double Y(void) { return y; }
  double Eta(void) { return eta; }
  double Px(void) { return px; }
  double Py(void) { return py; }
  double Rap(void) { return rap; }
  double Mt(void) { return mt; }
  double Vt_jet(void) { return vt_jet; }
  double Rap_jet(void) { return rap_jet; }
  // meaningful functions
  void propagate(double tau, double dtau);
  void addEnergyMom(double dE, double dPx, double dPy, double dPz);
};

// handles all the jets in the system
class Jets {
  std::vector<Jet> jets;

 public:
  Jets(char* filename);
  void propagate(double tau, double dtau);
  void checkJets(void);
  Jet* jet(int i) { return &jets[i]; }
  int nJets(void) { return jets.size(); }
};
