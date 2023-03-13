#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>

class Cell;
class Fluid;
class EoS;
class Hydro;

class VtkOutput {

 private:

 std::string  path_;
 /// Number of vtk output in current event
 int vtk_output_counter_ = -1;
 EoS* eos_;
 bool cartesian_;
 double xmin_, ymin_,zmin_;


 public:
  /**
   * Create a new VTK output.
   *
   * \param path Path to the output file.
   * \param name Name of the output.
   * \param out_par Additional information on the configured output.
   */
  VtkOutput( std::string  path, EoS* eos, double xmin,double ymin,double zmin, bool cartesian=false):
              path_(path),
              eos_(eos),
              cartesian_(cartesian),
              xmin_(xmin),
              ymin_(ymin),
              zmin_(zmin)
            {}

  void write(const Hydro h, std::string &quantity);
  void write_header(std::ofstream &file, const Hydro h,  const std::string &description);
  void write_vtk_scalar(std::ofstream &file, const Hydro h, std::string &quantity);
  std::string make_filename (const std::string &descr, int counter);

  
};