#include <iostream> 
#include <vector> 
#include <cmath>
#include <string>
#include <fstream>
#include "FDTD_2D.h"
#include <jsoncpp/json/value.h>
#include <json/json.h>

using namespace std;

int main(int argc, char* argv[]) {
  
  
  // const int Nx = 200;
  // const int Ny = 200;
  // const double t0 = 40.0;
  // const double spread = 12.0;
  
  FDTD_2D sim(argv[1]);
  sim.init();
  sim.print_initial_state();
  sim.fdtd_2d_basic_TEz_PML();

  return 0;
}