
class Nucleon {
 double x0, y0, eta0;
 public:
  double x, y, eta, rap, charge;
  Nucleon(double x, double y, double eta, double rap, double charge): x(x), x0(x), y(y), y0(y), eta(eta), eta0(eta), rap(rap), charge(charge) {}
  double getX0() { return x0; }
  double getY0() { return y0; }
  double getEta0() { return eta0; }
};
