#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

double EPSILON_0 = 8.8541878128e-12;
double MU_0      = 1.256637062e-6;
double c0        = 2.99792458e8;

double frequency = 1e9;
double center_wavelength = c0/frequency;
double simulation_size = 20 * center_wavelength;
int    cells_per_wavelength = 100;

double simulation_time = 1e-12;	
double pulse_width = 1e-14;
double pulse_delay = 4e-14;
double amplitude   = 1.0;

double step_size = center_wavelength / cells_per_wavelength;
double dx = step_size;
double dy = step_size;
int    N_space_cells = int(simulation_size / step_size);
double dt =  0.75 / ( c0 * sqrt( (1.0/(dx*dx)) + (1.0/(dy*dy)) ));
// int    N_time_steps = int(simulation_time / dt);

int    N_time_steps = 2000;
int data_capture_interval = 50;

vector<vector<double>> Ex, Ey, Hz;

int main() {

  cout<<"N_space_cells = "<<N_space_cells<<endl;
  cout<<"N_time_steps = "<<N_time_steps<<endl;
  cout<<"Step_size = "<<dx<<endl;

  int Nx = N_space_cells;
  int Ny = N_space_cells;

  Ex.resize(Nx, vector<double>(Ny, 0.0));
  Ey.resize(Nx, vector<double>(Ny, 0.0));
  Hz.resize(Nx, vector<double>(Ny, 0.0));

  // source position (center of domain)
  int src_i = Nx/2;
  int src_j = Ny/2;

  for(int t = 0; t < N_time_steps; t++) {

    // --- Update Hz ---
    for(int i = 0; i < Nx-1; i++) {
      for(int j = 0; j < Ny-1; j++) {
        Hz[i][j] -= (dt/MU_0) * ((Ey[i+1][j] - Ey[i][j]) / dx - (Ex[i][j+1] - Ex[i][j]) / dy);
      }
    }

    // --- Update Ex ---
    for(int i = 0; i < Nx; i++) {
      for(int j = 1; j < Ny; j++) {
        Ex[i][j] += (dt/EPSILON_0) * ((Hz[i][j] - Hz[i][j-1]) / dy);
      }
    }

    // --- Update Ey ---
    for(int i = 1; i < Nx; i++) {
      for(int j = 0; j < Ny; j++) {
        Ey[i][j] -= (dt/EPSILON_0) * ((Hz[i][j] - Hz[i-1][j]) / dx);
      }
    }

    // --- Inject Gaussian hard source into Ey ---
    double time = t * dt;
    double source = amplitude * exp(-pow((time - pulse_delay)/pulse_width, 2.0));
    Ey[src_i][int(0.25 * src_j)] += source;

    // (optional) monitor Ey at source
    cout<<"Iteration :"<<t<<endl;
    if(t % data_capture_interval == 0) {
      // cout<<"Iteration :"<<t<<endl;
      std::ofstream Hz_field("data/Hz/Hz" + to_string(t) + ".txt"); 
      for(int i = 0;i<N_space_cells;i++) {
        for(int j = 0;j<N_space_cells;j++) {
          Hz_field<<Hz[i][j]<<" ";
        }  
        Hz_field<<endl;
      }
      Hz_field.close();
    }
      
  }

  return 0;
}
