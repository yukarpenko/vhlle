#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>

class Hydro;



class VtkOutput {
 public:
  /**
   * Create a new VTK output.
   *
   * \param path Path to the output file.
   * \param name Name of the output.
   * \param out_par Additional information on the configured output.
   */
  VtkOutput(const std::string &path, const std::string &name,
            const std::string &out_par);
  ~VtkOutput();

  void write(Hydro &h);
  void write_header(std::ofstream &file, Hydro &h, const std::string &description);
  std::string make_filename (const std::string &descr, int counter);

  private:

  const std::string base_path_;
  const std::string name;
  /// Number of vtk output in current event
  int vtk_output_counter_ = -1;
};