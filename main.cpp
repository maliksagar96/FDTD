#include <iostream>
#include <vector>
#include <fstream> 
#include <cmath>
#include <string>

using namespace std;

void writeMatrixToFile(const std::vector<std::vector<double>>& field, const std::string& filename) {
  std::ofstream fout(filename);
  if (!fout.is_open()) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
  }

  for (const auto& row : field) {
      for (const auto& val : row) {
          fout << val << " ";
      }
      fout << "\n";
  }

  fout.close();
}

int main() {

//Defining cell limits
double x_min = 0.0, x_max = 2.0;
double y_min = 0.0, y_max = 2.0;

//Number of cells in x and y directions
int Nx = 200; 
int Ny = 200; 

//x and y spacing based on x and y range and Number of grid points in x and y direction. 
double dx = (x_max - x_min) / Nx;
double dy = (y_max - y_min) / Ny;

// Yee grid field allocations
std::vector<std::vector<double>> Ex(Nx, std::vector<double>(Ny, 0.0));
std::vector<std::vector<double>> Ey(Nx, std::vector<double>(Ny, 0.0));
std::vector<std::vector<double>> Hz(Nx, std::vector<double>(Ny, 0.0));

// for (int j = 0; j < Ny; ++j) {
//   Ey[10][j] = 10*exp(-pow((j - Ny/2.0) / 10.0, 2));  // Gaussian in Y at x=1
// }

writeMatrixToFile(Ey, "Ey_output.txt");

const double epsilon = 8.854e-12;   // Permittivity of free space
const double mu = 4 * M_PI * 1e-7;  // Permeability of free space

double c = 1.0 / sqrt(mu * epsilon);
double dt = 1.0 / (c * sqrt(1.0/(dx*dx) + 1.0/(dy*dy)));

dt *= 0.5;
int Nt = 10000;
double T0 = 20;
double tau = 6;
//Time marching.
Ey[Nx/2][Ny/2] = 0.5;
for(int t = 0;t<Nt;t++) {
  //Update Hz
  
  for (int i = 0; i < Nx-1; ++i){
    for (int j = 0; j < Ny-1; ++j){
      Hz[i][j] = Hz[i][j] - (dt / (mu * dx)) * (Ey[i+1][j] - Ey[i][j]) + (dt / (mu * dy)) * (Ex[i][j+1] - Ex[i][j]);
    }      
  }
    
  for (int i = 0; i < Nx; ++i) {
    for (int j = 1; j < Ny; ++j) {  // j starts at 1
        Ex[i][j] = Ex[i][j] + (dt / (epsilon * dy)) * (Hz[i][j] - Hz[i][j-1]);
    }
}

for (int i = 1; i < Nx; ++i) {  // i starts at 1
    for (int j = 0; j < Ny; ++j) {
        Ey[i][j] = Ey[i][j] - (dt / (epsilon * dx)) * (Hz[i][j] - Hz[i-1][j]);
    }
}

writeMatrixToFile(Ey, "data/Ey_" + std::to_string(t) + ".txt");


}



  return 0;
}