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
double frequency = 1e9;
double center_wavelength = c0/frequency;
double simulation_size = 30 * center_wavelength;
int cells_per_wavelength = 50;

double pulse_width = 1e-14;
double pulse_delay = 4e-14;
double amplitude   = 1.0;

double step_size = center_wavelength / cells_per_wavelength;
double dx = step_size;
double dy = step_size;
int N_space_cells = static_cast<int>(simulation_size / step_size);
double dt = 0.75 / (c0 * sqrt((1.0/(dx*dx)) + (1.0/(dy*dy))));

int N_time_steps = 2000;
int data_capture_interval = 50;

// --- Threading variables ---
queue<pair<int, vector<double>>> write_queue;
mutex queue_mutex;
condition_variable cv;
bool finished = false;

// --- Writer thread function ---
void writer_thread_func() {
  while (!finished || !write_queue.empty()) {
    unique_lock<mutex> lock(queue_mutex);
    cv.wait(lock, []{ return finished || !write_queue.empty(); });

    while (!write_queue.empty()) {
      auto [timestep, Hz_snapshot] = write_queue.front();
      write_queue.pop();
      lock.unlock(); // unlock while writing

      ofstream Hz_field("../data/Hz/Hz" + to_string(timestep) + ".txt");
      int N = sqrt(Hz_snapshot.size());
      for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
          Hz_field << Hz_snapshot[i*N + j] << " ";
          }
          Hz_field << "\n";
      }
      Hz_field.close();

      lock.lock(); // lock again for queue
    }
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

  N_time_steps = root["Iterations"].asInt();
  data_capture_interval = root["data_capture_interval"].asInt();

  int Nx = N_space_cells;
  int Ny = N_space_cells;

  // Initialize fields (1D contiguous storage)
  vector<double> Ex(Nx * Ny, 0.0);
  vector<double> Ey(Nx * Ny, 0.0);
  vector<double> Hz(Nx * Ny, 0.0);

  // Source position
  int src_i = Nx / 2;
  int src_j = 3 * Ny / 4;

  double h_coeff = dt / (sqrt(MU_0*EPSILON_0) * dx);
  double e_coeff = (dt / (sqrt(MU_0*EPSILON_0) * dx));

  // double h_coeff = dt / (MU_0 * dx);
  // double e_coeff = (dt / (EPSILON_0 * dx));

  double omega = 2*M_PI*frequency;  

  cout<<"*************************************************************"<<endl;
  cout << "Cells = " << N_space_cells * N_space_cells << endl;
  cout << "N_time_steps = " << N_time_steps << endl;
  cout << "dx = " << dx << endl;
  cout << "dt = " << dt << endl;
  cout<<"H_coeff = "<<h_coeff<<endl;
  cout<<"E_coeff = "<<e_coeff<<endl;
  cout<<"*************************************************************"<<endl;

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
        Ex[idx] += e_coeff * (Hz[idx] - Hz[i*Ny + j-1]);
      }
    }

    // --- Update Ey ---
    for(int i = 1; i < Nx; ++i) {
      for(int j = 0; j < Ny; ++j) {
        int idx = i*Ny + j;
        Ey[idx] -= e_coeff * (Hz[idx] - Hz[(i-1)*Ny + j]);
      }
    }

    // --- Inject Gaussian hard source into Ey ---
    double time = t * dt;
    // double source = amplitude * exp(-pow((time - pulse_delay)/pulse_width, 2.0)) * sin(omega * time);
    double source = sin(omega * time);
    Ey[src_i*Ny + src_j] += source;

    cout << "Iteration: " << t << endl;

    // --- Capture data asynchronously ---
    if(t % data_capture_interval == 0) {
      std::ofstream Hz_field("../data/Hz/Hz" + to_string(t) + ".txt"); 
      std::ofstream Ex_field("../data/Ex/Ex" + to_string(t) + ".txt"); 
      std::ofstream Ey_field("../data/Ey/Ey" + to_string(t) + ".txt"); 

      for(int i = 0; i < N_space_cells; ++i) {
        for(int j = 0; j < N_space_cells; ++j) {
          Hz_field << Hz[i*Ny + j] << " "; 
          Ex_field << Ex[i*Ny + j] << " "; 
          Ey_field << Ey[i*Ny + j] << " "; 
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
