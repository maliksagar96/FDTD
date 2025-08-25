#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <json/json.h>

using namespace std;

// --- Physical constants ---
double EPSILON_0 = 8.8541878128e-12;
double MU_0      = 1.256637062e-6;
double c0        = 2.99792458e8;

// --- Simulation parameters ---
double frequency, center_wavelength, simulation_size;
int cells_per_wavelength;
double pulse_width, pulse_delay, amplitude;
double step_size;
double dx, dy, dt;
int vacuum_cells, N_time_steps, data_capture_interval;
double pml_cond_e;

vector<double> sigma_x,sigma_y, sigma_h_x, sigma_h_y;
vector<double> Ca_x,Cb_x, Ca_y,Cb_y;

int Nx, Ny, pml_size;

void init(Json::Value root) {
  frequency = root["frequency"].asDouble();
  N_time_steps = root["Iterations"].asInt();
  data_capture_interval = root["data_capture_interval"].asInt();
  pml_size = root["pml_size"].asInt();
  
  
  
  center_wavelength = c0/frequency;
  simulation_size = 30 * center_wavelength;
  cells_per_wavelength = 50;
  pulse_width = dt*100;
  pulse_delay = 4*pulse_width;
  amplitude = 1.0;
  step_size = center_wavelength / cells_per_wavelength;
  dx = step_size;
  dy = step_size;
  vacuum_cells = static_cast<int>(simulation_size / step_size);
  dt = 0.75 / (c0 * sqrt((1.0/(dx*dx)) + (1.0/(dy*dy))));
  
  Nx = vacuum_cells + pml_size;
  Ny = vacuum_cells + pml_size;

  sigma_x.resize(Nx, 0);
  sigma_y.resize(Ny, 0);
  sigma_h_x.resize(Nx, 0);
  sigma_h_y.resize(Ny, 0);
  Ca_x.resize(Nx, 0);
  Cb_x.resize(Nx, 0);
  Ca_y.resize(Ny, 0);
  Cb_y.resize(Ny, 0);

  pml_cond_e = -log(1e-3) * c0 * EPSILON_0 / (2 * pml_size * dx);
  
  for(int i = vacuum_cells;i<Nx;i++) {
    sigma_x[i] = pml_cond_e;
    sigma_y[i] = pml_cond_e;
  }
  

}

int main() {
 std::ifstream file("../control_file.scf", std::ifstream::binary);
  if (!file.is_open()) {
    std::cerr << "Error: could not open file\n";
    return 1;
  }

  Json::Value root;
  file >> root; 

  init(root);

  // Initialize fields (1D contiguous storage)
  vector<double> Ex(Nx * Ny, 0.0);
  vector<double> Ey(Nx * Ny, 0.0);
  vector<double> Hz(Nx * Ny, 0.0);

  // Source position
  int src_i = Nx / 2;
  int src_j = 3 * Ny / 4;

  // double h_coeff = dt / (MU_0 * dx);
  // double e_coeff = (dt / (EPSILON_0 * dx));

  double h_coeff = dt / (sqrt(MU_0*EPSILON_0) * dx);
  double e_coeff = (dt / (sqrt(MU_0*EPSILON_0) * dx));

  ofstream sig_x("sigma_x.txt");
  ofstream sig_y("sigma_y.txt");
  ofstream cax("cax.txt");
  ofstream cay("cay.txt");
  ofstream cbx("cbx.txt");
  ofstream cby("cby.txt");

  for (int i = 0; i < Nx; ++i) {
    double s = sigma_x[i];
    double denom = 1.0 + (s * dt) / (2.0 * EPSILON_0);
    Ca_x[i] = (1.0 - (s * dt) / (2.0 * EPSILON_0)) / denom;
    Cb_x[i] = (dt / (sqrt(MU_0*EPSILON_0)*dx)) / denom;
   
    sig_x<<sigma_x[i]<<endl;
    cax<<Ca_x[i]<<endl;
    cbx<<Cb_x[i]<<endl;    
  }

  for (int j = 0; j < Ny; ++j) {
    double s = sigma_y[j];
    double denom = 1.0 + (s * dt) / (2.0 * EPSILON_0);
    Ca_y[j] = (1.0 - (s * dt) / (2.0 * EPSILON_0)) / denom;
    Cb_y[j] = (dt / (sqrt(MU_0*EPSILON_0)*dy)) / denom;        

    sig_y<<sigma_y[j]<<endl;
    cay<<Ca_y[j]<<endl;
    cby<<Cb_y[j]<<endl;    
  }

  sig_x.close();
  sig_y.close();
  cax.close();
  cay.close();
  cbx.close();
  cby.close();

  double omega = 2*M_PI*frequency;  

  cout<<"*************************************************************"<<endl;
  
  cout << "N_time_steps = " << N_time_steps << endl;
  cout << "dx = " << dx << endl;
  cout << "dt = " << dt << endl;
  cout<<"H_coeff = "<<h_coeff<<endl;
  cout<<"E_coeff = "<<e_coeff<<endl;
  cout<<"*************************************************************"<<endl;
  double t0 = 100 * dt;
  double tau = 20*dt;

// --- FDTD main loop ---
  for(int t = 0; t < N_time_steps; t++) {

    // --- Update Hz ---
    for(int i = 0; i < Nx-1; ++i) {
      for(int j = 0; j < Ny-1; ++j) {
        int idx = i*Ny + j;
        Hz[idx] -= h_coeff * ((Ey[(i+1)*Ny + j] - Ey[idx]) - (Ex[i*Ny + j+1] - Ex[idx]));
      }
    }

  // --- Update Ex ---
    for(int i = 0; i < Nx; ++i) {
      for(int j = 1; j < Ny; ++j) {
        int idx = i*Ny + j;
        Ex[idx] = Ca_x[i]*Ex[idx] + Cb_x[i] * (Hz[idx] - Hz[i*Ny + j-1]);
        // Ex[idx] += e_coeff * (Hz[idx] - Hz[i*Ny + j-1]);
      }
    }

    // --- Update Ey ---
    for (int i = 1; i < Nx; ++i) {
      for (int j = 0; j < Ny; ++j) {
        int idx = i * Ny + j;
        Ey[idx] = Ca_y[j] * Ey[idx] -
                  Cb_y[j] * (Hz[idx] - Hz[(i - 1) * Ny + j]);
      }
    }

    // --- Inject Gaussian hard source into Ey ---
    double time = t * dt;
    // double source = amplitude * exp(-pow((time - pulse_delay)/pulse_width, 2.0));
    // double source = sin(omega * time);
    double source = sin(omega * time) * exp(-pow((time - t0)/tau, 2));
    Ey[src_i*Ny + src_j] += source;

    cout << "Iteration: " << t << endl;

    // --- Capture data asynchronously ---
    if(t % data_capture_interval == 0) {
      std::ofstream Hz_field("../data/Hz/Hz" + to_string(t) + ".txt"); 
      std::ofstream Ex_field("../data/Ex/Ex" + to_string(t) + ".txt"); 
      std::ofstream Ey_field("../data/Ey/Ey" + to_string(t) + ".txt"); 

      for(int i = 0; i < Nx; ++i) {
        for(int j = 0; j < Ny; ++j) {
          Hz_field << Hz[i*Nx + j] << " "; 
          Ex_field << Ex[i*Nx + j] << " "; 
          Ey_field << Ey[i*Nx + j] << " "; 
        }
        Hz_field << endl;
        Ex_field << endl;
        Ey_field << endl;
      }
      Hz_field.close();
      Ex_field.close();
      Ey_field.close();
    }
  }

  return 0;
}
