#ifndef NUMERICAL_INTEGRATION_H
#define NUMERICAL_INTEGRATION_H
#include <vector>

void GetCoefs2DLaguerre32Legendre32(double ay, double by, std::vector<double> &xlag, std::vector<double> &wlag, std::vector<double> &xleg, std::vector<double> &wleg);

void GetCoefs2DLegendre32Legendre32(double ay, double by, double a2y, double b2y, std::vector<double> &xlag, std::vector<double> &wlag, std::vector<double> &xleg, std::vector<double> &wleg);

void GetCoefsIntegrateLegendre32(double a, double b, std::vector<double> &x, std::vector<double> &w);

void GetCoefsIntegrateLegendre10(double a, double b, std::vector<double> &x, std::vector<double> &w);

void GetCoefsIntegrateLegendre5(double a, double b, std::vector<double> &x, std::vector<double> &w);

void GetCoefsIntegrateLegendre40(double a, double b, std::vector<double> &x, std::vector<double> &w);

void GetCoefsIntegrateLaguerre15(std::vector<double> &x, std::vector<double> &w);

void GetCoefsIntegrateLaguerre32(std::vector<double> &x, std::vector<double> &w);

#endif
