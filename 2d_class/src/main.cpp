#include "fdtd_2d.h"
#include <json/json.h>
#include <fstream>

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <control_file.scf>" << std::endl;
    return 1;
  }

  std::ifstream file(argv[1], std::ifstream::binary);
  if (!file) {
    std::cerr << "Error: could not open file " << argv[1] << std::endl;
    return 1;
  }

  Json::Value root;
  file >> root;

  FDTD_2D solver(root);
  // solver.run_TEz();
  solver.run_TMz();
  
  return 0;
}
