#include <iostream> 
#include <vector> 
#include <cmath>
#include <string>
#include <fstream>
#include <FDTD_1D.h>

using namespace std;

int main(int argc, char* argv[]) {
  
  FDTD_1D sim(argv[1]);
  sim.init();
  // sim.fdtd_1d_basic();
  // sim.fdtd_1d_graded_pml();
  sim.fdtd_1d_graded_Cpml();
    
  return 0;
}
