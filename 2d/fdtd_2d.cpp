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

double pulse_width = 1e-14;
double pulse_delay = 4e-14;
double amplitude   = 1.0;

double step_size = center_wavelength / cells_per_wavelength;
double dx = step_size;
double dy = step_size;
int    N_space_cells = static_cast<int>(simulation_size / step_size);
double dt =  0.75 / (c0 * sqrt( (1.0/(dx*dx)) + (1.0/(dy*dy)) ));

int    N_time_steps = 2000;
int data_capture_interval = 50;

int main() {    
  int Nx = N_space_cells;
  int Ny = N_space_cells;

  // Correctly initialize single 1D vectors for contiguous memory
  vector<double> Ex(Nx * Ny, 0.0);
  vector<double> Ey(Nx * Ny, 0.0);
  vector<double> Hz(Nx * Ny, 0.0);

  // source position (center of domain)
  int src_i = Nx/2;
  int src_j = 3*Ny/4;

  double h_coeff = (dt/(MU_0*dx));
  double e_coeff = (dt/(EPSILON_0*dx));

  cout<<"*************************************************************"<<endl;
  cout << "N_space_cells = " << N_space_cells << endl;
  cout << "N_time_steps = " << N_time_steps << endl;
  cout << "Step_size = " << dx << endl;
  cout << "dt = " << dt << endl;
  cout<<"H_coeff = "<<h_coeff<<endl;
  cout<<"E_coeff = "<<e_coeff<<endl;
  cout<<"*************************************************************"<<endl;

  for(int t = 0; t < N_time_steps; t++) {
    // --- Update Hz (using E from previous time step) ---
    for(int i = 0; i < Nx-1; ++i) {
      for(int j = 0; j < Ny-1; ++j) {
        int idx = i*Ny + j;
        Hz[idx] -= h_coeff * ((Ey[(i+1)*Ny + j] - Ey[idx]) - (Ex[i*Ny + j+1] - Ex[idx]));
      }
    }

    // --- Update Ex (using H from current time step) ---
    for(int i = 0; i < Nx; ++i) {
      for(int j = 1; j < Ny; ++j) {
        int idx = i*Ny + j;
        Ex[idx] += e_coeff * (Hz[idx] - Hz[i*Ny + j-1]);
      }
    }

    // --- Update Ey (using H from current time step) ---
    for(int i = 1; i < Nx; ++i) {
      for(int j = 0; j < Ny; ++j) {
        int idx = i*Ny + j;
        Ey[idx] -= e_coeff * (Hz[idx] - Hz[(i-1)*Ny + j]);
      }
    }

    // --- Inject Gaussian hard source into Ey ---
    double time = t * dt;
    double source = amplitude * exp(-pow((time - pulse_delay)/pulse_width, 2.0));
    Ey[src_i*Ny + src_j] += source; // Corrected index

    // (optional) monitor Ey at source
    cout << "Iteration :" << t << endl;
    if(t % data_capture_interval == 0) {
      std::ofstream Hz_field("data/Hz/Hz" + to_string(t) + ".txt"); 
      for(int i = 0; i < N_space_cells; ++i) {
        for(int j = 0; j < N_space_cells; ++j) {
          Hz_field << Hz[i*Ny + j] << " "; // Corrected index
        }
        Hz_field << endl;
      }
      Hz_field.close();
    }
  }

  return 0;
}