class TGraph;

class CrossSections {
 TGraph *gSigmaT, *gSigmaE, *gSigmaPiN;
public:
 CrossSections(void);
 ~CrossSections(void);
 void NN(double Ekin, double& sigmaT, double& sigmaE, double& sigmaP);
 double piN(double s);
};
