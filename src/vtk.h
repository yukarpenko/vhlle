#pragma once

#include <cstddef>
#include <filesystem>
#include <iosfwd>
#include <map>
#include <string>

class EoS;
class Fluid;
class Hydro;

class VtkOutput {
 private:
  enum class QuantityType {
    Scalar,
    Vector,
    Tensor
  };

  std::filesystem::path output_path_;
  EoS &eos_;
  double xmin_, ymin_, etamin_;
  std::size_t num_of_cells_x_direction_{0}, num_of_cells_y_direction_{0},
              num_of_cells_eta_direction_{0};
  int vtk_output_counter_{0};  // Number of vtk output in current event

  /* The following map contains all currently supported VTK quantities. If
   * multiple quantities are desired, the delimiter in the config file has to be
   * a comma. For example: `VTK_output_values eps,mub,nq,T,v,pi`. Any whitespace
   * will be stripped, i.e. `VTK_output_values T, eps` would also be possible.
   * For vector quantities only the corresponding three vector will be written
   * to the output file. */
  inline static const std::map<std::string, QuantityType> valid_quantities_ = {
    {"eps", QuantityType::Scalar},        // energy density
    {"mub", QuantityType::Scalar},        // baryon chemical potential
    {"muq", QuantityType::Scalar},        // electric chemical potential
    {"mus", QuantityType::Scalar},        // strangeness chemical potential
    {"nb", QuantityType::Scalar},         // baryon density
    {"nq", QuantityType::Scalar},         // charge density
    {"ns", QuantityType::Scalar},         // strangeness density
    {"p", QuantityType::Scalar},          // pressure
    {"Pi", QuantityType::Scalar},         // bulk pressure
    {"pi", QuantityType::Tensor},         // shear stress tensor
    {"T", QuantityType::Scalar},          // temperature
    {"v", QuantityType::Vector}           // velocity
  };

  void write_header(std::ofstream &file, Fluid &fluid,
                    const std::string &quantity);
  void write_vtk_scalar(std::ofstream &file, Fluid &fluid, const double tau,
                        const std::string &quantity);
  void write_vtk_vector(std::ofstream &file, Fluid &fluid, const double tau,
                        const std::string &quantity);
  void write_vtk_tensor(std::ofstream &file, Fluid &fluid,
                        const std::string &quantity);
  std::filesystem::path make_filename(const std::string &quantity,
                                      int counter) const;

 public:
  /**
   * Create a new VTK output.
   *
   * \param output_path Directory in which VTK files are written.
   * \param eos The equation of state.
   * \param xmin x coordinate of the first cell (center of fluid cell).
   * \param ymin y coordinate of the first cell (center of fluid cell).
   * \param etamin eta coordinate of the first cell (center of fluid cell).
   */
  VtkOutput(std::string output_path, EoS &eos, double xmin, double ymin,
            double etamin):
              output_path_(output_path),
              eos_(eos),
              xmin_(xmin),
              ymin_(ymin),
              etamin_(etamin)
            {}

  void write(const Hydro &h, const std::string &quantities);
};
