#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "eos.h"
#include "fld.h"
#include "hdo.h"
#include "vtk.h"

void VtkOutput::write_header(std::ofstream &file, const Hydro h,
                             const std::string &description) {
  num_of_cells_x_direction_ = h.getFluid()->getNX();
  num_of_cells_y_direction_ = h.getFluid()->getNY();
  num_of_cells_eta_direction_ = h.getFluid()->getNZ();
  file << "# vtk DataFile Version 2.0\n"
    << description << "\n"
    << "ASCII\n"
    << "DATASET STRUCTURED_POINTS\n"
    << "DIMENSIONS " << num_of_cells_x_direction_ << " "
      << num_of_cells_y_direction_ << " " << num_of_cells_eta_direction_ << "\n"
    << "SPACING " << h.getFluid()->getDx() << " " << h.getFluid()->getDy()
                  << " " << h.getFluid()->getDz() << "\n"
    << "ORIGIN " << xmin_ << " " << ymin_ << " " << etamin_ << "\n"
    << "POINT_DATA " << num_of_cells_x_direction_ * num_of_cells_y_direction_
                        * num_of_cells_eta_direction_ << "\n";

  return;
}

std::string VtkOutput::make_filename(const std::string &descr, int counter) {
  char suffix[24];
  snprintf(suffix, sizeof(suffix), "_taustep%05i.vtk",counter);
  return path_ + std::string("/") + descr + std::string(suffix);
}

void VtkOutput::write_vtk_scalar(std::ofstream &file, const Hydro h,
                                 const std::string &quantity) {
  file << "SCALARS " << quantity << " double 1\n"
       << "LOOKUP_TABLE default\n";
  file << std::setprecision(3);
  file << std::fixed;

  for (int ieta = 0; ieta < num_of_cells_eta_direction_; ieta++) {
    for (int iy = 0; iy < num_of_cells_y_direction_; iy++) {
      for (int ix = 0; ix < num_of_cells_x_direction_; ix++) {
        double e, nb, nq, ns, p, vx, vy, vz;
        Cell* cell = h.getFluid()->getCell(ix,iy,ieta);
        cell->getPrimVar(eos_, h.getTau(), e, p, nb, nq, ns, vx, vy, vz);
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
          eos_->eos(e, nb, nq, ns, T, mub, muq, mus, p);
          if (quantity == "mub") {
            q = mub;
          } else if (quantity == "muq") {
            q = muq;
          } else if (quantity == "mus") {
            q = mus;
          } else if (quantity == "T") {
            q = T;
          }
        }
        file << q << " ";
      }
      file << "\n";
    }
  }
}

void VtkOutput::write_vtk_vector(std::ofstream &file, const Hydro h,
                                 const std::string &quantity) {
  file << "VECTORS " << quantity << " double\n";
  file << std::setprecision(3);
  file << std::fixed;

  for (int ieta = 0; ieta < num_of_cells_eta_direction_; ieta++) {
    for (int iy = 0; iy < num_of_cells_y_direction_; iy++) {
      for (int ix = 0; ix < num_of_cells_x_direction_; ix++) {
        double e, p, nb, nq, ns, vx, vy, vz;
        Cell* cell = h.getFluid()->getCell(ix,iy,ieta);
        cell->getPrimVar(eos_, h.getTau(), e, p, nb, nq, ns, vx, vy, vz);
        std::vector<double> q = {0.,0.,0.};
        if (quantity == "v") {
          q = {vx, vy, vz};
        }
        file << q.at(0) << " " << q.at(1) << " " << q.at(2) << "\n";
      }
    }
  }
}

void VtkOutput::write_vtk_tensor(std::ofstream &file, const Hydro h,
                                 const std::string &quantity) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++ ) {
      file << "SCALARS " << quantity << std::to_string(i) << std::to_string(j)
           << " double 1\n"
           << "LOOKUP_TABLE default\n";
      file << std::setprecision(3);
      file << std::fixed;

      for (int ieta = 0; ieta < num_of_cells_eta_direction_; ieta++) {
        for (int iy = 0; iy < num_of_cells_y_direction_; iy++) {
          for (int ix = 0; ix < num_of_cells_x_direction_; ix++) {
            Cell* cell = h.getFluid()->getCell(ix,iy,ieta);
            double q = 0;
            if (quantity == "pi") {
              q = cell->getpi(i,j);
            }
            file << q << " ";
          }
          file << "\n";
        }
      }
    }
  }
}

std::vector<std::string> split (const std::string &s, char delim) {
  std::vector<std::string> result;
  std::stringstream ss(s);
  std::string item;

  while (getline(ss, item, delim)) {
    result.push_back(item);
  }

  return result;
}

bool VtkOutput::is_quantity_implemented(const std::string &quantity) {
  bool quantity_is_valid = (valid_quantities_.find(quantity)
                            != valid_quantities_.end());
  return quantity_is_valid;
}

void VtkOutput::write(const Hydro h, const std::string &quantities) {
  std::vector<std::string> quantities_list = split(quantities,',');
  for (std::string q : quantities_list){
    if (!is_quantity_implemented(q)) {
      std::cout << "Given quantity '" << q << "' is not an "
        "implemented VTK quantity. This entry will be skipped." << std::endl;
      continue;
    }
    std::ofstream file;

    file.open(make_filename(q, vtk_output_counter_), std::ios::out);
    write_header(file, h, q);
    if (valid_quantities_.at(q) == "scalar") {
      write_vtk_scalar(file, h, q);
    } else if (valid_quantities_.at(q) == "vector") {
      write_vtk_vector(file, h, q);
    } else if (valid_quantities_.at(q) == "tensor") {
      write_vtk_tensor(file, h, q);
    } else {
      std::cout << "Quantity '" << q << "' is neither stated to be a scalar, "
        "nor a vector, nor a tensor. Skipping this quantity. Please check the "
        "map in file src/vtk.h." << std::endl;
    }
    file.close();
  }

  vtk_output_counter_++;
  return;
}
