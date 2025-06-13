/*

This simulation is done using the first chapter of the book Electromagntic simulation by Dennis M. Sullivan.

Question - Why does Hy travel in this way?

*/


#include <iostream> 
#include <vector> 
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

void writeArrayToFile(const std::vector<double>& arr, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error opening file: " << filename << "\n";
        return;
    }

    for (double val : arr) {
        outFile << val << "\n";
    }

    outFile.close();
}


int main() {

  int Nt = 1000;

  int Nx = 200;
  int Nc = Nx/2;


  double t0 = 40.0, spread = 12, T = 0;
  int NSTEPS = 1;
  
  const double epsilon = 8.854e-12;   // Permittivity of free space
  const double mu = 4 * M_PI * 1e-7;  // Permeability of free space

  vector<double> Ex(Nx, 0);
  vector<double> Hy(Nx, 0);

  for(int t = 0;t<Nt;t++) {
    for(int i = 0;i<Nx;i++) {
      Ex[i] += 0.5*(Hy[i-1] - Hy[i]);
    }
  
    Ex[Nc] = exp(-0.5*(pow((t0 - t)/spread, 2)));

    for(int i = 0;i<Nx;i++) {
      Hy[i] += 0.5*(Ex[i] - Ex[i+1]);
    }

    string filename = "data/Hy_" + to_string(t) + ".txt";
    writeArrayToFile(Hy, filename);

    

  }



  return 0;
}