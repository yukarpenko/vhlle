#include <iostream>
#include <algorithm>
#include <vector>
#include <stdexcept> 
#include <string>

#include "eos.h"
#include "fld.h"
#include "hdo.h"

using Matrix = std::vector<std::vector<double>>;
using Cube = std::vector<std::vector<std::vector<Matrix>>>;

class Vorticity {
 public:
  // Standard constructor sets all elements of dbeta_ and dbetaCube_ to zero
  Vorticity()
      : dbeta_(4, std::vector<double>(4, 0.0)),
        dbetaCube_(2, std::vector<std::vector<Matrix>>(
                          2, std::vector<Matrix>(
                                 2, Matrix(4, std::vector<double>(4, 0.0))))),
        dbetaAveraged_(4, std::vector<double>(4, 0.0)) {}

  // Get the thermal vorticity tensor
  Matrix dbeta() const { return dbeta_; }

  // Get a single element of the thermal vorticity tensor
  double dbeta(int i, int j) const {
    if (i >= 0 && i < 4 && j >= 0 && j < 4) {
      return dbeta_[i][j];
    } else {
      throw std::invalid_argument("dbeta index out of range");
    }
  }

  // Set the vorticity tensor for a single cell
  void set_dbeta(const Matrix& dbeta) {
    if (dbeta.size() == 4 && dbeta[0].size() == 4) {
      for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
          if (dbeta[i][j] > 100) {
            std::cerr
                << "Warning: dbeta component greater than 100, setting to 100."
                << std::endl;
            dbeta_[i][j] = 100;
          } else if (dbeta[i][j] < -100) {
            std::cerr
                << "Warning: dbeta component less than -100, setting to -100."
                << std::endl;
            dbeta_[i][j] = -100;
          } else {
            dbeta_[i][j] = dbeta[i][j];
          }
        }
      }
    } else {
      throw std::invalid_argument("dbeta must have shape (4, 4)");
    }
  }

  void calculateDbeta(const int column_index, Fluid& f, EoS& eos,
                      const FluidState& current, const FluidState& previous,
                      const double dt, const double tau) {
    // current cell
    double T, mub, muq, mus, p;
    eos.eos(current.eosData.e, current.eosData.nb, current.eosData.nq,
            current.eosData.ns, T, mub, muq, mus, p);

    // previous cell
    double T_0, mub_0, muq_0, mus_0, p_0;
    eos.eos(previous.eosData.e, previous.eosData.nb, previous.eosData.nq,
            previous.eosData.ns, T_0, mub_0, muq_0, mus_0, p_0);

    FourVector dX = {dt, f.getDx(), f.getDy(), f.getDz()};

    switch (column_index) {
      case 0:
        // calculate dbeta_[0][...] components
        dbeta_[0][0] = (current.u[0] / T - previous.u[0] / T_0) / dX[0];
        dbeta_[0][1] = (current.u[1] / T - previous.u[1] / T_0) / dX[0];
        dbeta_[0][2] = (current.u[2] / T - previous.u[2] / T_0) / dX[0];
        dbeta_[0][3] = (current.u[3] / T - previous.u[3] / T_0) / dX[0];
        break;

      case 1:
        // calculate dbeta_[1][...] components
        dbeta_[1][0] = 0.5 * (current.u[0] / T - previous.u[0] / T_0) / dX[1];
        dbeta_[1][1] = 0.5 * (current.u[1] / T - previous.u[1] / T_0) / dX[1];
        dbeta_[1][2] = 0.5 * (current.u[2] / T - previous.u[2] / T_0) / dX[1];
        dbeta_[1][3] = 0.5 * (current.u[3] / T - previous.u[3] / T_0) / dX[1];
        break;

      case 2:
        // calculate dbeta_[2][...] components
        dbeta_[2][0] = 0.5 * (current.u[0] / T - previous.u[0] / T_0) / dX[2];
        dbeta_[2][1] = 0.5 * (current.u[1] / T - previous.u[1] / T_0) / dX[2];
        dbeta_[2][2] = 0.5 * (current.u[2] / T - previous.u[2] / T_0) / dX[2];
        dbeta_[2][3] = 0.5 * (current.u[3] / T - previous.u[3] / T_0) / dX[2];
        break;

      case 3:
        // calculate dbeta_[2][...] components
        dbeta_[3][0] = 0.5 * (current.u[0] / T - previous.u[0] / T_0) /
                       (dX[3] * (tau + 0.5 * dt));
        dbeta_[3][1] = 0.5 * (current.u[1] / T - previous.u[1] / T_0) /
                       (dX[3] * (tau + 0.5 * dt));
        dbeta_[3][2] = 0.5 * (current.u[2] / T - previous.u[2] / T_0) /
                       (dX[3] * (tau + 0.5 * dt));
        dbeta_[3][3] = 0.5 * (current.u[3] / T - previous.u[3] / T_0) /
                       (dX[3] * (tau + 0.5 * dt));
        break;

      default:
        // If dbeta_index is not 0, 1, 2, or 3, raise an error.
        throw std::invalid_argument("Invalid dbeta_index: " +
                                    std::to_string(column_index));
    }
  }

  // Getter for the dbeta_ tensor of the current cell.
  const Matrix& getDbeta() const { return dbeta_; }

  // Set the dbeta vector to zero
  void resetDbeta() {
    for (auto& row : dbeta_) {
      std::fill(row.begin(), row.end(), 0.0);  // Efficiently reset each row
    }
  }

  // Calculates the dbetaCube_ for a 2x2x2 cube of adjacent cells.
  void calculateDbetaCube(
      const std::vector<std::vector<std::vector<Matrix>>>& hydroCells);

  // Interpolates the dbetaCube_ values using weighting factors and stores the
  // result in dbetaCentral_.
  void averageDbetaCube(const std::vector<double>& weightingFactors);

  // Print the vorticity tensor to file
  void print_to_file(const std::string& filename);

 private:
  // Thermal vorticity tensor for a single cell.
  Matrix dbeta_;

  // Thermal vorticity tensor for a 2x2x2 cube of adjacent cells.
  Cube dbetaCube_;

  // Averaged dbeta from dbetaCube_
  Matrix dbetaAveraged_;
};