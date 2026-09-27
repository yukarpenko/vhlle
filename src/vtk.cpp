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

void VtkOutput::write_header(std::ofstream &file, const Hydro &h,
                             const std::string &quantity) {
  num_of_cells_x_direction_ = h.getFluid()->getNX();
  num_of_cells_y_direction_ = h.getFluid()->getNY();
  num_of_cells_eta_direction_ = h.getFluid()->getNZ();
  const std::size_t point_count = num_of_cells_x_direction_ *
                                  num_of_cells_y_direction_ *
                                  num_of_cells_eta_direction_;
  file << "# vtk DataFile Version 2.0\n"
    << quantity << "\n"
    << "ASCII\n"
    << "DATASET STRUCTURED_POINTS\n"
    << "DIMENSIONS " << num_of_cells_x_direction_ << " "
      << num_of_cells_y_direction_ << " " << num_of_cells_eta_direction_ << "\n"
    << "SPACING " << h.getFluid()->getDx() << " " << h.getFluid()->getDy()
                  << " " << h.getFluid()->getDz() << "\n"
    << "ORIGIN " << xmin_ << " " << ymin_ << " " << etamin_ << "\n"
    << "POINT_DATA " <<  point_count << "\n";
}

std::string VtkOutput::make_filename(const std::string &quantity,
                                     const int counter) const {
  char suffix[24];
  std::snprintf(suffix, sizeof(suffix), "_taustep%05d.vtk", counter);
  return path_ + "/" + quantity + suffix;
}

void VtkOutput::write_vtk_scalar(std::ofstream &file, const Hydro &h,
                                 const std::string &quantity) {
  file << "SCALARS " << quantity << " double 1\n"
       << "LOOKUP_TABLE default\n";
  file << std::setprecision(3);
  file << std::fixed;

  for (std::size_t ieta = 0; ieta < num_of_cells_eta_direction_; ieta++) {
    for (std::size_t iy = 0; iy < num_of_cells_y_direction_; iy++) {
      for (std::size_t ix = 0; ix < num_of_cells_x_direction_; ix++) {
        double e, nb, nq, ns, p, vx, vy, vz;
        Cell* cell = h.getFluid()->getCell(ix, iy, ieta);
        if (cartesian_) {
         cell->getPrimVar(&eos_, 1.0, e, p, nb, nq, ns, vx, vy, vz);
        } else {
         cell->getPrimVar(&eos_, h.getTau(), e, p, nb, nq, ns, vx, vy, vz);
        }
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
        file << q << " ";
      }
      file << "\n";
    }
  }
}

void VtkOutput::write_vtk_vector(std::ofstream &file, const Hydro &h,
                                 const std::string &quantity) {
  // The only implemented vector quantity is velocity.
  if (quantity != "v") {
    throw std::logic_error("Unhandled vector VTK quantity: '" + quantity + "'");
  }
  file << "VECTORS " << quantity << " double\n";
  file << std::setprecision(3);
  file << std::fixed;

  for (std::size_t ieta = 0; ieta < num_of_cells_eta_direction_; ieta++) {
    for (std::size_t iy = 0; iy < num_of_cells_y_direction_; iy++) {
      for (std::size_t ix = 0; ix < num_of_cells_x_direction_; ix++) {
        double e, p, nb, nq, ns, vx, vy, vz;
        Cell* cell = h.getFluid()->getCell(ix, iy, ieta);
        if (cartesian_) {
         cell->getPrimVar(&eos_, 1.0, e, p, nb, nq, ns, vx, vy, vz);
        } else {
         cell->getPrimVar(&eos_, h.getTau(), e, p, nb, nq, ns, vx, vy, vz);
        }
        file << vx << ' ' << vy << ' ' << vz << '\n';
      }
    }
  }
}

void VtkOutput::write_vtk_tensor(std::ofstream &file, const Hydro &h,
                                 const std::string &quantity) {
  // The only implemented tensor quantity is the shear stress tensor.
  if (quantity != "pi") {
    throw std::logic_error("Unhandled tensor VTK quantity: '" + quantity + "'");
  }
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      file << "SCALARS " << quantity << i << j << " double 1\n"
           << "LOOKUP_TABLE default\n";
      file << std::setprecision(3);
      file << std::fixed;

      for (std::size_t ieta = 0; ieta < num_of_cells_eta_direction_; ieta++) {
        for (std::size_t iy = 0; iy < num_of_cells_y_direction_; iy++) {
          for (std::size_t ix = 0; ix < num_of_cells_x_direction_; ix++) {
            Cell* cell = h.getFluid()->getCell(ix, iy, ieta);
            file << cell->getpi(i, j) << " ";
          }
          file << "\n";
        }
      }
    }
  }
}

std::vector<std::string> split(const std::string &s, const char delim) {
  std::vector<std::string> result;
  std::stringstream ss{s};
  std::string item;

  while (std::getline(ss, item, delim)) {
    result.push_back(item);
  }

  return result;
}

void VtkOutput::write(const Hydro &h, const std::string &quantities) {
  if (quantities.empty()) {
    std::cerr << "VTK quantities in configuration file are empty. "
      "No VTK output will be created.\n";
    return;
  }
  const std::vector<std::string> quantities_list = split(quantities, ',');
  for (const auto &q : quantities_list) {
    const auto quantity_it = valid_quantities_.find(q);
    if (quantity_it == valid_quantities_.end()) {
      std::cerr << "Given quantity '" << q << "' is not an "
        "implemented VTK quantity. This entry will be skipped.\n";
      continue;
    }
    std::ofstream file{};

    file.open(make_filename(q, vtk_output_counter_), std::ios::out);
    write_header(file, h, q);
    switch (quantity_it->second) {
    case QuantityType::Scalar:
      write_vtk_scalar(file, h, q);
      break;
    case QuantityType::Vector:
      write_vtk_vector(file, h, q);
      break;
    case QuantityType::Tensor:
      write_vtk_tensor(file, h, q);
      break;
    default:
      throw std::logic_error("Unhandled VTK quantity type '" + q + "'.");
    }
    file.close();
  }

  vtk_output_counter_++;
}
