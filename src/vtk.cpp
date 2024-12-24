#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "eos.h"
#include "fld.h"
#include "hdo.h"
#include "vtk.h"

void VtkOutput::write_header(std::ofstream &file, const Hydro h, const std::string &description) {
  if (cartesian_) {
    file << "# vtk DataFile Version 2.0\n"
      << description << "\n"
      << "ASCII\n"
      << "DATASET STRUCTURED_POINTS\n"
      << "DIMENSIONS " << h.getFluid()->getNX() << " " << h.getFluid()->getNY() << " " << (h.getFluid()->getNZ())*20 << "\n"
      << "SPACING " << h.getFluid()->getDx() << " " << h.getFluid()->getDy() << " " << h.getTau()*std::sinh(h.getFluid()->getDz()*(h.getFluid()->getNZ()/2))/(20*h.getFluid()->getNZ()/2) << "\n"
      << "ORIGIN " << xmin_ << " " << ymin_ << " " << etamin_ << "\n"
      << "POINT_DATA " << h.getFluid()->getNX()*h.getFluid()->getNY()*(h.getFluid()->getNZ())*20 << "\n";
  } else {
    file << "# vtk DataFile Version 2.0\n"
      << description << "\n"
      << "ASCII\n"
      << "DATASET STRUCTURED_POINTS\n"
      << "DIMENSIONS " << h.getFluid()->getNX() << " " << h.getFluid()->getNY() << " " << h.getFluid()->getNZ() << "\n"
      << "SPACING " << h.getFluid()->getDx() << " " << h.getFluid()->getDy() << " " << h.getFluid()->getDz() << "\n"
      << "ORIGIN " << xmin_ << " " << ymin_ << " " << etamin_ << "\n"
      << "POINT_DATA " << h.getFluid()->getNX()*h.getFluid()->getNY()*h.getFluid()->getNZ() << "\n";
  }

  return;
}

std::string VtkOutput::make_filename(const std::string &descr, int counter) {
  char suffix[24];
  snprintf(suffix, sizeof(suffix), "_taustep%05i.vtk",counter);
  return path_ + std::string("/") + descr + std::string(suffix);
}

std::vector<double> VtkOutput::smearing_factor_and_poseta(const Hydro h, const int iz, const int z_length) {
  double poseta = iz;
  double factor = 1;
  if (cartesian_) {
    double total_length = h.getTau()*std::sinh(h.getFluid()->getDz()*(h.getFluid()->getNZ()/2));
    double pos = 0;
    if (iz < z_length/2) {
      pos = -total_length+total_length/(z_length/2)*iz;
      for (int ieta = 0; ieta < h.getFluid()->getNZ()/2; ieta++) {
        if (h.getTau()*std::sinh(h.getFluid()->getDz()*(ieta-h.getFluid()->getNZ()/2)) > pos) {
          factor = fabs((h.getTau()*std::sinh(h.getFluid()->getDz()*(ieta-h.getFluid()->getNZ()/2))-pos)/pos);
          poseta = ieta;
          break;
        }
      }
    } else {
      pos = total_length/(z_length/2)*(iz-z_length/2);
      for (int ieta = h.getFluid()->getNZ()/2; ieta < h.getFluid()->getNZ(); ieta++) {
        if (h.getTau()*std::sinh(h.getFluid()->getDz()*(ieta-h.getFluid()->getNZ()/2)) > pos) {
          factor = fabs((h.getTau()*std::sinh(h.getFluid()->getDz()*(ieta-h.getFluid()->getNZ()/2))-pos)/pos);
          poseta = ieta;
          break;
        }
      }
    }
    if (pos == 0) {
      factor = 1;
    }
  }
  std::vector<double> factor_and_poseta = {factor, poseta};
  return factor_and_poseta;
}

void VtkOutput::write_vtk_scalar(std::ofstream &file, const Hydro h,
                                 std::string &quantity) {
  file << "SCALARS " << quantity << " double 1\n"
       << "LOOKUP_TABLE default\n";
  file << std::setprecision(3);
  file << std::fixed;
  int z_length = h.getFluid()->getNZ();
  if (cartesian_) {
    z_length = 20*z_length;
  }

  for (int iz = 0; iz < z_length; iz++) {
    std::vector<double> factor_and_poseta = smearing_factor_and_poseta(h, iz, z_length);
    double factor = factor_and_poseta.at(0);
    double poseta = factor_and_poseta.at(1);
    for (int iy = 0; iy < h.getFluid()->getNY(); iy++) {
      for (int ix = 0; ix < h.getFluid()->getNX(); ix++) {
        double e, mub, muq, mus, nb, nq, ns, p, T, vx, vy, vz;
        double e2, mub2, muq2, mus2, nb2, nq2, ns2, p2, T2, vx2, vy2, vz2;
        Cell* cell = h.getFluid()->getCell(ix,iy,poseta);
        Cell* cell2;
        if (poseta > 0) {
          cell2 = h.getFluid()->getCell(ix,iy,poseta-1);
        } else {
          cell2 = cell;
        }
        cell->getPrimVar(eos_, h.getTau(), e, p, nb, nq, ns, vx, vy, vz);
        cell2->getPrimVar(eos_, h.getTau(), e2, p2, nb2, nq2, ns2, vx2, vy2, vz2);
        double q = 0;
        if (quantity == "eps") {
          q = factor*e+(factor-1)*e2;
        } else if (quantity == "mub") {
          eos_->eos(e, nb, nq, ns, T, mub, muq, mus, p);
          eos_->eos(e2, nb2, nq2, ns2, T2, mub2, muq2, mus2, p2);
          q = factor*mub+(factor-1)*mub2;
        } else if (quantity == "muq") {
          eos_->eos(e, nb, nq, ns, T, mub, muq, mus, p);
          eos_->eos(e2, nb2, nq2, ns2, T2, mub2, muq2, mus2, p2);
          q = factor*muq+(factor-1)*muq2;
        } else if (quantity == "mus") {
          eos_->eos(e, nb, nq, ns, T, mub, muq, mus, p);
          eos_->eos(e2, nb2, nq2, ns2, T2, mub2, muq2, mus2, p2);
          q = factor*mus+(factor-1)*mus2;
        } else if (quantity == "nb") {
          q = factor*nb+(factor-1)*nb2;
        } else if (quantity == "nq") {
          q = factor*nq+(factor-1)*nq2;
        } else if (quantity == "ns") {
          q = factor*ns+(factor-1)*ns2;
        } else if (quantity == "p") {
          q = factor*p+(factor-1)*p2;
        } else if (quantity == "Pi") {
          double Pi = cell->getPi();
          double Pi2 = cell2->getPi();
          q = factor*Pi+(factor-1)*Pi2;
        } else if (quantity == "T") {
          eos_->eos(e, nb, nq, ns, T, mub, muq, mus, p);
          eos_->eos(e2, nb2, nq2, ns2, T2, mub2, muq2, mus2, p2);
          q = factor*T+(factor-1)*T2;
        }
        file << q << " ";
      }
      file << "\n";
    }
  }
}

void VtkOutput::write_vtk_vector(std::ofstream &file, const Hydro h,
                                  std::string &quantity) {
  file << "VECTORS " << quantity << " double\n";
  file << std::setprecision(3);
  file << std::fixed;
  int z_length = h.getFluid()->getNZ();
  if (cartesian_) {
    z_length = 20*z_length;
  }

  for (int iz = 0; iz < z_length; iz++) {
    std::vector<double> factor_and_poseta = smearing_factor_and_poseta(h, iz, z_length);
    double factor = factor_and_poseta.at(0);
    double poseta = factor_and_poseta.at(1);
    for (int iy = 0; iy < h.getFluid()->getNY(); iy++) {
      for (int ix = 0; ix < h.getFluid()->getNX(); ix++) {
        double e, p, nb, nq, ns, vx, vy, vz, T, mub, muq, mus;
        double e2, p2, nb2, nq2, ns2, vx2, vy2, vz2, T2, mub2, muq2, mus2;
        Cell* cell = h.getFluid()->getCell(ix,iy,poseta);
        Cell* cell2;
        if (poseta > 0) {
          cell2 = h.getFluid()->getCell(ix,iy,poseta-1);
        } else {
          cell2 = cell;
        }
        cell->getPrimVar(eos_, h.getTau(), e, p, nb, nq, ns, vx, vy, vz);
        cell2->getPrimVar(eos_, h.getTau(), e2, p2, nb2, nq2, ns2, vx2, vy2, vz2);
        std::vector<double> q = {0.,0.,0.};
        if (quantity == "v") {
          q = {factor*vx+(factor-1)*vx2, factor*vy+(factor-1)*vy2, factor*vz+(factor-1)*vz2};
        }
        file << q.at(0) << " " << q.at(1) << " " << q.at(2) << "\n";
      }
    }
  }
}

void VtkOutput::write_vtk_tensor(std::ofstream &file, const Hydro h,
                                  std::string &quantity) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++ ) {
      file << "SCALARS " << quantity << std::to_string(i) << std::to_string(j)
           << " double 1\n"
           << "LOOKUP_TABLE default\n";
      file << std::setprecision(3);
      file << std::fixed;
      int z_length = h.getFluid()->getNZ();
      if (cartesian_) {
        z_length = 20*z_length;
      }

      for (int iz = 0; iz < z_length; iz++) {
        std::vector<double> factor_and_poseta = smearing_factor_and_poseta(h, iz, z_length);
        double factor = factor_and_poseta.at(0);
        double poseta = factor_and_poseta.at(1);
        for (int iy = 0; iy < h.getFluid()->getNY(); iy++) {
          for (int ix = 0; ix < h.getFluid()->getNX(); ix++) {
            Cell* cell = h.getFluid()->getCell(ix,iy,poseta);
            Cell* cell2;
            if (poseta > 0) {
              cell2 = h.getFluid()->getCell(ix,iy,poseta-1);
            } else {
              cell2 = cell;
            }
            double q = 0;
            if (quantity == "pi") {
              double pi_ij = cell->getpi(i,j);
              double pi2_ij = cell2->getpi(i,j);
              q = factor*pi_ij+(factor-1)*pi2_ij;
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
  std::stringstream ss (s);
  std::string item;

  while (getline(ss, item, delim)) {
    result.push_back(item);
  }

  return result;
}

bool VtkOutput::is_quantity_implemented(std::string &quantity) {
  bool quantity_is_valid = (valid_quantities_.find(quantity)
                            != valid_quantities_.end());
  return quantity_is_valid;
}

void VtkOutput::write(const Hydro h, std::string &quantities) {
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
