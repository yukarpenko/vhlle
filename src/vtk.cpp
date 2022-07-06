#include "vtk.h"
#include <stdexcept>

/*struct FileDeleter {
  /// The class has no members, so this is a noop.
  constexpr FileDeleter() = default;

  
  void operator()(std::FILE* f) const {
    if (f == nullptr) {
      return;
    }
    if (0 != std::fclose(f)) {
      throw std::runtime_error(std::strerror(errno));
    }
  }
};

using FilePtr = std::unique_ptr<std::FILE, FileDeleter>;*/

VtkOutput::VtkOutput(const std::string &path, const std::string &name,
                     const std::string &out_par)
    : name(name),
      base_path_(path){}

VtkOutput::~VtkOutput() {}


void VtkOutput::write_header(std::ofstream &file, Hydro &h, const std::string &description){
    file << "# vtk DataFile Version 2.0\n"
       << description << "\n";
       //<< "ASCII\n"
       //<< "DATASET STRUCTURED_POINTS\n"
       //<< "DIMENSIONS " << dim[0] << " " << dim[1] << " " << dim[2] << "\n"
       //<< "SPACING " << cs[0] << " " << cs[1] << " " << cs[2] << "\n"
       //<< "ORIGIN " << orig[0] << " " << orig[1] << " " << orig[2] << "\n"
       //<< "POINT_DATA " << lattice.size() << "\n";
    return;
}

std::string VtkOutput::make_filename(const std::string &descr, int counter) {
  char suffix[22];
  snprintf(suffix, sizeof(suffix), "_tstep%05i.vtk",counter);
  return base_path_ + std::string("/") + descr + std::string(suffix);
}

void VtkOutput::write(Hydro &h){
    char filename[32];
    std::string t="test";
    std::ofstream file;
    file.open(make_filename(t, vtk_output_counter_++), std::ios::out);
    write_header(file, h, t);
    return;
}
