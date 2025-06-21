/*

This simulation is done using the first chapter of the book Electromagntic simulation by Dennis M. Sullivan.

We are doing the 

*/

#include <iostream> 
#include <vector> 
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

void fdtd_1d_basic(vector<double>& ex, vector<double>& hy) {
  
  for(int t = 0;t<Nt;t++) {
    for(int i = 0;i<Nx;i++) {
      ex[i] += 0.5*(hy[i-1] - hy[i]);
    }
  
    ex[Nc] = exp(-0.5*(pow((t0 - t)/spread, 2)));

    for(int i = 0;i<Nx;i++) {
      hy[i] += 0.5*(ex[i] - ex[i+1]);
    }

    string filename = "../data/Hy_1D_" + to_string(t) + ".txt";
    writeArrayToFile(hy, filename);
    
  }

}

class 1D_FDTD{
  public:
    1D_FDTD();
    ~1D_FDTD();
    void set_time_steps();
    void fdtd_1d_basic();
    void writeArrayToFile(const std::vector<double>& arr, const std::string& filename);
    void set_excitation();

  private:
    vector<double> ex, vector<double> hy;
    double t0, spread;
    int nt;
}

void 1D_FDTD::set_time_steps(int nt) {
  this.nt = nt;
}

void 1D_FDTD::fdtd_1d_basic(){
  for(int t = 0;t<Nt;t++) {
    for(int i = 0;i<Nx;i++) {
      ex[i] += 0.5*(hy[i-1] - hy[i]);
    }
  
    ex[Nc] = exp(-0.5*(pow((t0 - t)/spread, 2)));

    for(int i = 0;i<Nx;i++) {
      hy[i] += 0.5*(ex[i] - ex[i+1]);
    }    
  }    
}

void 1D_FDTD::writeArrayToFile(const std::vector<double>& arr, const std::string& filename) {
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

void 1D_FDTD::set_excitation(double t0, double spread) {
  
}

int main() {

  int Nt = 1000;

  int Nx = 200;
  int Nc = Nx/2;

  double t0 = 40.0, spread = 12, T = 0;
  int NSTEPS = 1;
  
  const double epsilon = 8.854e-12;   // Permittivity of free space
  const double mu = 4 * M_PI * 1e-7;  // Permeability of free space

  vector<double> ex(Nx, 0);
  vector<double> hy(Nx, 0);

  fdtd_1d_basic(ex, hy);

  return 0;
}