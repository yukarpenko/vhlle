#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "eos.h"
#include "fld.h"
#include "hdo.h"
#include "vtk.h"

namespace {

/// Helper function to set the precision for VTK ASCII output.
void set_precision_vtk_ascii_output(std::ostream &file) {
  file << std::defaultfloat << std::setprecision(6);
}

/**
 * Magnitudes below this (in code units: GeV, fm and products thereof) are
 * round-off noise or physically empty tails and are written as exact 0.
 */
constexpr double vtk_output_zero_threshold = 1e-6;

/// Sets tiny values (and -0) to 0; NaN and inf pass through deliberately.
double set_tiny_values_to_zero(double value) {
  return std::abs(value) < vtk_output_zero_threshold ? 0.0 : value;
}

/// Helper function to split the string of quantities.
std::vector<std::string> split(const std::string &s, const char delim) {
  std::vector<std::string> result{};
  std::stringstream ss{s};
  std::string item{};

  while (std::getline(ss, item, delim)) {
    const auto first = item.find_first_not_of(" \t\n\r\f\v");
    const auto last = item.find_last_not_of(" \t\n\r\f\v");
    if (first == std::string::npos) {
      result.emplace_back();
    } else {
      result.emplace_back(item.substr(first, last - first + 1));
    }
  }
  return result;
}

}  // namespace

void VtkOutput::write_header(std::ofstream &file, Fluid &fluid,
                             const std::string &quantity) {
  const std::size_t point_count = num_of_cells_x_direction_ *
                                  num_of_cells_y_direction_ *
                                  num_of_cells_eta_direction_;
  file << "# vtk DataFile Version 2.0\n"
    << quantity << "\n"
    << "ASCII\n"
    << "DATASET STRUCTURED_POINTS\n"
    << "DIMENSIONS " << num_of_cells_x_direction_ << " "
      << num_of_cells_y_direction_ << " " << num_of_cells_eta_direction_ << "\n"
    << "SPACING " << fluid.getDx() << " " << fluid.getDy()
                  << " " << fluid.getDz() << "\n"
    << "ORIGIN " << xmin_ << " " << ymin_ << " " << etamin_ << "\n"
    << "POINT_DATA " <<  point_count << "\n";
}

std::filesystem::path VtkOutput::make_filename(const std::string &quantity,
                                               const int counter) const {
  std::ostringstream filename{};
  filename << quantity
           << "_taustep"
           << std::setw(5)
           << std::setfill('0')
           << counter
           << ".vtk";
  return output_path_ / filename.str();
}

void VtkOutput::write_vtk_scalar(std::ofstream &file, Fluid &fluid,
                                 const double tau,
                                 const std::string &quantity) {
  file << "SCALARS " << quantity << " double 1\n"
       << "LOOKUP_TABLE default\n";

  for (std::size_t ieta = 0; ieta < num_of_cells_eta_direction_; ieta++) {
    for (std::size_t iy = 0; iy < num_of_cells_y_direction_; iy++) {
      for (std::size_t ix = 0; ix < num_of_cells_x_direction_; ix++) {
        double e, nb, nq, ns, p, vx, vy, vz;
        Cell* cell = fluid.getCell(ix, iy, ieta);
        cell->getPrimVar(&eos_, tau, e, p, nb, nq, ns, vx, vy, vz);
        double q = 0;
        // scalar quantities
        if (quantity == "eps") {
          q = e;
        } else if (quantity == "nb") {
          q = nb;
        } else if (quantity == "nq") {
          q = nq;
        } else if (quantity == "ns") {
          q = ns;
        } else if (quantity == "p") {
          q = p;
        } else if (quantity == "Pi") {
          q = cell->getPi();
        // scalar quantities that need eos()
        } else if (quantity == "mub" || quantity == "muq" || quantity == "mus"
                   || quantity == "T") {
          double mub, muq, mus, T;
          eos_.eos(e, nb, nq, ns, T, mub, muq, mus, p);
          if (quantity == "mub") {
            q = mub;
          } else if (quantity == "muq") {
            q = muq;
          } else if (quantity == "mus") {
            q = mus;
          } else if (quantity == "T") {
            q = T;
          }
        } else {
          throw std::logic_error(
              "Unhandled scalar VTK quantity: '" + quantity + "'");
        }
        if (ix > 0) {
          file.put(' ');
        }
        file << set_tiny_values_to_zero(q);
      }
      file.put('\n');
    }
  }
}

void VtkOutput::write_vtk_vector(std::ofstream &file, Fluid &fluid,
                                 const double tau,
                                 const std::string &quantity) {
  // The only implemented vector quantity is velocity.
  if (quantity != "v") {
    throw std::logic_error("Unhandled vector VTK quantity: '" + quantity + "'");
  }
  file << "VECTORS " << quantity << " double\n";

  for (std::size_t ieta = 0; ieta < num_of_cells_eta_direction_; ieta++) {
    for (std::size_t iy = 0; iy < num_of_cells_y_direction_; iy++) {
      for (std::size_t ix = 0; ix < num_of_cells_x_direction_; ix++) {
        double e, p, nb, nq, ns, vx, vy, vz;
        Cell* cell = fluid.getCell(ix, iy, ieta);
        cell->getPrimVar(&eos_, tau, e, p, nb, nq, ns, vx, vy, vz);
        file << set_tiny_values_to_zero(vx) << ' '
             << set_tiny_values_to_zero(vy) << ' '
             << set_tiny_values_to_zero(vz) << '\n';
      }
    }
  }
}

void VtkOutput::write_vtk_tensor(std::ofstream &file, Fluid &fluid,
                                 const std::string &quantity) {
  // The only implemented tensor quantity is the shear stress tensor.
  if (quantity != "pi") {
    throw std::logic_error("Unhandled tensor VTK quantity: '" + quantity + "'");
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      file << "SCALARS " << quantity << i << j << " double 1\n"
           << "LOOKUP_TABLE default\n";

      for (std::size_t ieta = 0; ieta < num_of_cells_eta_direction_; ieta++) {
        for (std::size_t iy = 0; iy < num_of_cells_y_direction_; iy++) {
          for (std::size_t ix = 0; ix < num_of_cells_x_direction_; ix++) {
            Cell* cell = fluid.getCell(ix, iy, ieta);
            if (ix > 0) {
              file.put(' ');
            }
            file << set_tiny_values_to_zero(cell->getpi(i, j));
          }
          file.put('\n');
        }
      }
    }
  }
}

void VtkOutput::write(const Hydro &h, const std::string &quantities) {
  if (quantities.empty()) {
    std::cerr << "VTK quantities in configuration file are empty. "
      "No VTK output will be created.\n";
    return;
  }

  Fluid &fluid = *h.getFluid();
  const double tau = h.getTau();
  num_of_cells_x_direction_ = fluid.getNX();
  num_of_cells_y_direction_ = fluid.getNY();
  num_of_cells_eta_direction_ = fluid.getNZ();
  const std::vector<std::string> quantities_list = split(quantities, ',');

  for (const auto &q : quantities_list) {
    const auto quantity_it = valid_quantities_.find(q);
    if (quantity_it == valid_quantities_.end()) {
      std::cerr << "Given quantity '" << q << "' is not an "
        "implemented VTK quantity. This entry will be skipped.\n";
      continue;
    }

    const auto filename = make_filename(q, vtk_output_counter_);
    std::ofstream file{filename};
    if (!file) {
      throw std::runtime_error(
        "Unable to open VTK output file: " + filename.string());
    }

    set_precision_vtk_ascii_output(file);
    write_header(file, fluid, q);
    switch (quantity_it->second) {
    case QuantityType::Scalar:
      write_vtk_scalar(file, fluid, tau, q);
      break;
    case QuantityType::Vector:
      write_vtk_vector(file, fluid, tau, q);
      break;
    case QuantityType::Tensor:
      write_vtk_tensor(file, fluid, q);
      break;
    default:
      throw std::logic_error("Unhandled VTK quantity type '" + q + "'.");
    }

    file.flush();
    if (!file) {
      throw std::runtime_error(
          "Failed while writing VTK output file: " + filename.string());
    }
  }

  ++vtk_output_counter_;
}
