#pragma once

#include <cstddef>
#include <iosfwd>
#include <map>
#include <string>

class EoS;
class Hydro;

class VtkOutput {
 private:
  std::string path_;
  EoS* eos_;
  double xmin_, ymin_, etamin_;
  std::size_t num_of_cells_x_direction_{0}, num_of_cells_y_direction_{0},
              num_of_cells_eta_direction_{0};
  bool cartesian_;
  int vtk_output_counter_{0};  // Number of vtk output in current event

  /* The following map contains all currently supported VTK quantities.
   * If multiple quantities are desired, the delimiter in the config file has
   * to be a comma without any whitespaces in between, for example:
   * `VTK_output_valus eps,mub,nq,T,v,pi`
   * For vector quantities only the corresponding three vector will be written
   * to the output file. */
  const std::map<std::string, std::string> valid_quantities_ = {
    {"eps", "scalar"},        // energy density
    {"mub", "scalar"},        // baryon chemical potential
    {"muq", "scalar"},        // electric chemical potential
    {"mus", "scalar"},        // strangeness chemical potential
    {"nb", "scalar"},         // baryon density
    {"nq", "scalar"},         // charge density
    {"ns", "scalar"},         // strangeness density
    {"p", "scalar"},          // pressure
    {"Pi", "scalar"},         // bulk pressure
    {"pi", "tensor"},         // shear stress tensor
    {"T", "scalar"},          // temperature
    {"v", "vector"}           // velocity
  };

  void write_header(std::ofstream &file, const Hydro &h,
                    const std::string &quantity);
  void write_vtk_scalar(std::ofstream &file, const Hydro &h,
                        const std::string &quantity);
  void write_vtk_vector(std::ofstream &file, const Hydro &h,
                        const std::string &quantity);
  void write_vtk_tensor(std::ofstream &file, const Hydro &h,
                        const std::string &quantity);
  std::string make_filename(const std::string &quantity, int counter) const;

 public:
  /**
   * Create a new VTK output.
   *
   * \param path Path to the output file.
   * \param eos The equation of state.
   * \param xmin x coordinate of the first cell (center of fluid cell).
   * \param ymin y coordinate of the first cell (center of fluid cell).
   * \param etamin eta coordinate of the first cell (center of fluid cell).
   */
  VtkOutput(std::string path, EoS* eos, double xmin, double ymin,
            double etamin, bool cartesian):
              path_(path),
              eos_(eos),
              xmin_(xmin),
              ymin_(ymin),
              etamin_(etamin),
              cartesian_(cartesian)
            {}

  void write(const Hydro &h, const std::string &quantities);
};
