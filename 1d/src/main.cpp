

#include <iostream> 
#include <vector> 
#include <cmath>
#include <string>
#include <fstream>
#include <FDTD_1D.h>

using namespace std;


// Main function
int main() {
  const int Nt = 1000;
  const int Nx = 200;
  const double t0 = 40.0;
  const double spread = 12.0;

  
  FDTD_1D sim;
  sim.set_spacing(Nx);
  sim.set_time_steps(Nt);
  sim.set_excitation(t0, spread);
  // sim.fdtd_1d_basic();
  sim.fdtd_1d_abc();


  return 0;
}
