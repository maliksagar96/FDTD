#include <iostream> 
#include <vector> 
#include <cmath>
#include <string>
#include <fstream>

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

  int nt = 1000;

  int nx = 200;
  int ny = 200;


  double t0 = 20.0, spread = 12, T = 0;
  
  const double epsilon = 8.854e-12;   // Permittivity of free space
  const double mu = 4 * M_PI * 1e-7;  // Permeability of free space
  const double omega = 1;
  
std::vector<std::vector<double>> ez(nx, std::vector<double>(nx, 0.0));
std::vector<std::vector<double>> dz(nx, std::vector<double>(nx, 0.0));
std::vector<std::vector<double>> hx(nx, std::vector<double>(nx, 0.0));
std::vector<std::vector<double>> hy(nx, std::vector<double>(nx, 0.0));
std::vector<std::vector<double>> epsilon_r(nx, std::vector<double>(nx, 1.0));


  for(int t = 0;t<nt;t++) {
    //ez 
    for(int i = 1;i<nx;i++) { 
      for(int j = 1;j<nx;j++) {
        dz[i][j] += 0.5*(hy[i][j] - hy[i-1][j] - hx[i][j] + hx[i][j-1]);    
      }
    }
  
    for(int i = 1;i<nx;i++) { 
      for(int j = 1;j<nx;j++) {
        ez[i][j] = epsilon_r[i][j] * dz[i][j];    
      }
    }

    //ez = pulse
    ez[nx/2][nx/2] = exp(-0.5*(pow((t0 - t)/spread, 2)));
    // ez[nx/2][nx/2] = sin(omega * t);

    //hy and Hz
    for(int i = 0;i<nx-1;i++) {
      for(int j = 0;j<nx-1;j++) {
        hx[i][j] += 0.5*(ez[i][j] - ez[i][j+1]);
      }
    }

    
    for(int i = 0;i<nx-1;i++) {
      for(int j = 0;j<nx-1;j++) {
        hy[i][j] += 0.5*(ez[i+1][j] - ez[i][j]);
      }
    }

    writeMatrixToFile(ez, "../data/Ez/Ez_" + std::to_string(t) + ".txt");

  }



  return 0;
}