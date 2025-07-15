#include <iostream> 
#include <vector> 
#include <cmath>
#include <string>
#include <fstream>
#include <FDTD_1D.h>

using namespace std;


// Main function
int main(int argc, char* argv[]) {
  const int Nt = 1000;
  const int Nx = 200;
  const double t0 = 40.0;
  const double spread = 12.0;
  
  FDTD_1D sim(argv[1]);
  sim.init();
  // sim.fdtd_1d_basic();
  sim.fdtd_1d_graded_pml();
  // sim.fdtd_1d_graded_Cpml();
  
  
  return 0;
}
