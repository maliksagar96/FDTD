#include <iostream> 
#include <vector> 
#include <cmath>
#include <string>
#include <fstream>
#include "FDTD_2D.h"

using namespace std;

int main() {

const int Nt = 1000;
  const int Nx = 200;
  const int Ny = 200;
  const double t0 = 40.0;
  const double spread = 12.0;
  
  FDTD_2D sim;
  sim.set_spacing(Nx, Ny);
  sim.set_time_steps(Nt);
  sim.set_excitation(t0, spread);
  sim.fdtd_2d_basic();



  return 0;
}