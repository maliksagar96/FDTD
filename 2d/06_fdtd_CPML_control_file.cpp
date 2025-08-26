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

// --- Domain (simulation box) ---
double simulation_size[4];   // [x_min, y_min, x_max, y_max]
double source_point[2];      // source location (x,y)

// --- Grid & resolution ---
int Nx, Ny;                  // number of cells in x,y
int cells_per_wavelength;
double step_size;
double dx, dy;

// --- PML ---
int pml_size;
double pml_cond_e;
vector<double> sigma_x_r_pml, sigma_x_l_pml, sigma_y_t_pml, sigma_y_b_pml;
vector<double> kappa_x_l, kappa_x_r, kappa_y_b, kappa_y_t;
vector<double> a_x_r, a_x_l, a_y_b, a_y_t;
vector<double> be_x_r, be_x_l, be_y_t, be_y_b;
vector<double> ce_x_r, ce_x_l, ce_y_t, ce_y_b;
vector<double> inv_kappa_x, inv_kappa_y;

// --- Material / conductivity arrays (per-axis) ---
vector<double> sigma_x, sigma_y;
vector<double> sigma_h_x, sigma_h_y;

// --- Time / CFL ---
double CFL;
double dt;
int N_time_steps;
int data_capture_interval;

// --- Source / waveform parameters ---
double frequency;
double center_wavelength;
double pulse_width, pulse_delay, amplitude;

// --- CPML helper fields ---
vector<double> Psi_Ey_x, Psi_Ex_y, Psi_Hz_y, Psi_Hz_x;

// --- Update coefficients per cell ---
vector<double> Ca_x, Cb_x, Ca_y, Cb_y;

// --- Derived / bookkeeping grid counts ---
int vacuum_cells_x, vacuum_cells_y;

// --- Async I/O (writer queue) ---
std::queue<std::tuple<int, int, int,
                      std::vector<double>,
                      std::vector<double>,
                      std::vector<double>>> write_queue;
std::mutex queue_mutex;
std::condition_variable cv;
bool finished = false;


void initPML() {

  double m = 4;
  double kappa_max = 5.0;  
  double a_max     = 0.5;  

  pml_cond_e = -(m + 1) * log(1e-9) * c0 * EPSILON_0 / (2 * pml_size * dx);

  // ---------- Right (+x) PML ----------
  for (int i = 0; i < pml_size; i++) {
    double x = static_cast<double>(i) / pml_size;

    sigma_x_r_pml[i] = pml_cond_e * pow(x, m);
    kappa_x_r[i]     = 1.0 + (kappa_max - 1.0) * pow(x, m);
    a_x_r[i]         = a_max * ((pml_size - i) / static_cast<double>(pml_size));

    be_x_r[i] = exp(-(sigma_x_r_pml[i] / kappa_x_r[i] + a_x_r[i]) * dt / EPSILON_0);
    ce_x_r[i] = sigma_x_r_pml[i] * (be_x_r[i] - 1.0) /
                ((sigma_x_r_pml[i] + kappa_x_r[i] * a_x_r[i]) * kappa_x_r[i]);

    inv_kappa_x[Nx - pml_size + i] = 1.0 / (kappa_x_r[i] * dx);
  }

  // ---------- Left (−x) PML ----------
  for (int i = 0; i < pml_size; i++) {
    double x = static_cast<double>(pml_size - i) / pml_size;

    sigma_x_l_pml[i] = pml_cond_e * pow(x, m);
    kappa_x_l[i]     = 1.0 + (kappa_max - 1.0) * pow(x, m);
    a_x_l[i]         = a_max * (static_cast<double>(i) / pml_size);

    be_x_l[i] = exp(-(sigma_x_l_pml[i] / kappa_x_l[i] + a_x_l[i]) * dt / EPSILON_0);
    ce_x_l[i] = sigma_x_l_pml[i] * (be_x_l[i] - 1.0) /
                ((sigma_x_l_pml[i] + kappa_x_l[i] * a_x_l[i]) * kappa_x_l[i]);

    inv_kappa_x[i] = 1.0 / (kappa_x_l[i] * dx);
  }

  // ---------- Top (+y) PML ----------
  for (int j = 0; j < pml_size; j++) {
    double y = static_cast<double>(j) / pml_size;

    sigma_y_t_pml[j] = pml_cond_e * pow(y, m);
    kappa_y_t[j]     = 1.0 + (kappa_max - 1.0) * pow(y, m);
    a_y_t[j]         = a_max * ((pml_size - j) / static_cast<double>(pml_size));

    be_y_t[j] = exp(-(sigma_y_t_pml[j] / kappa_y_t[j] + a_y_t[j]) * dt / EPSILON_0);
    ce_y_t[j] = sigma_y_t_pml[j] * (be_y_t[j] - 1.0) /
                ((sigma_y_t_pml[j] + kappa_y_t[j] * a_y_t[j]) * kappa_y_t[j]);

    inv_kappa_y[Ny - pml_size + j] = 1.0 / (kappa_y_t[j] * dy);
  }

  // ---------- Bottom (−y) PML ----------
  for (int j = 0; j < pml_size; j++) {
    double y = static_cast<double>(pml_size - j) / pml_size;

    sigma_y_b_pml[j] = pml_cond_e * pow(y, m);
    kappa_y_b[j]     = 1.0 + (kappa_max - 1.0) * pow(y, m);
    a_y_b[j]         = a_max * (static_cast<double>(j) / pml_size);

    be_y_b[j] = exp(-(sigma_y_b_pml[j] / kappa_y_b[j] + a_y_b[j]) * dt / EPSILON_0);
    ce_y_b[j] = sigma_y_b_pml[j] * (be_y_b[j] - 1.0) /
                  ((sigma_y_b_pml[j] + kappa_y_b[j] * a_y_b[j]) * kappa_y_b[j]);

    inv_kappa_y[j] = 1.0 / (kappa_y_b[j] * dy);
  }
}


void init(Json::Value root) {
  
  frequency = root["frequency"].asDouble();
  N_time_steps = root["Iterations"].asInt();
  data_capture_interval = root["data_capture_interval"].asInt();
  pml_size = root["pml_size"].asInt();  

  simulation_size[0] = root["simulation_size"][0].asDouble();
  simulation_size[1] = root["simulation_size"][1].asDouble();
  simulation_size[2] = root["simulation_size"][2].asDouble();
  simulation_size[3] = root["simulation_size"][3].asDouble();
  cells_per_wavelength = root["cells_per_wavelength"].asInt();
  amplitude = root["amplitude"].asDouble();
  CFL = root["CFL"].asDouble();
  source_point[0] = root["source_point"][0].asDouble();
  source_point[1] = root["source_point"][1].asDouble();

  if((source_point[0] < simulation_size[0] || source_point[0] > simulation_size[2] || 
      source_point[1] < simulation_size[1] || source_point[1] > simulation_size[3])) {
    cerr<<"Source point not inside the simualtion domain."<<endl;
    exit(1);
  }

  center_wavelength = c0/frequency;      
  step_size = center_wavelength / cells_per_wavelength;
  dx = step_size;
  dy = step_size;
  dt = CFL / (c0 * sqrt((1.0/(dx*dx)) + (1.0/(dy*dy))));
  vacuum_cells_x = static_cast<int>(abs(simulation_size[2] - simulation_size[0]) / dx);
  vacuum_cells_y = static_cast<int>(abs(simulation_size[1] - simulation_size[3]) / dy);
  pulse_width = dt*100;
  pulse_delay = 4*pulse_width;
  
  Nx = vacuum_cells_x + pml_size;
  Ny = vacuum_cells_y + pml_size;

  sigma_x.resize(Nx, 0);
  sigma_y.resize(Ny, 0);
  sigma_h_x.resize(Nx, 0);
  sigma_h_y.resize(Ny, 0);
  
  Ca_x.resize(Nx, 1);
  Cb_x.resize(Nx, dt/EPSILON_0);
  Ca_y.resize(Ny, 1);
  Cb_y.resize(Ny, dt/EPSILON_0);

  inv_kappa_x.resize(Nx, 1/dx);
  inv_kappa_y.resize(Ny, 1/dy);

  Psi_Ex_y.resize(Nx * Ny, 0.0);
  Psi_Ey_x.resize(Nx * Ny, 0.0);
  Psi_Hz_y.resize(Nx * Ny, 0.0);
  Psi_Hz_x.resize(Nx * Ny, 0.0);

  sigma_x_r_pml.resize(pml_size, 0);
  sigma_x_l_pml.resize(pml_size, 0);
  sigma_y_t_pml.resize(pml_size, 0);
  sigma_y_b_pml.resize(pml_size, 0);

  kappa_x_r.resize(pml_size, 0);
  kappa_x_l.resize(pml_size, 0);
  kappa_y_t.resize(pml_size, 0);
  kappa_y_b.resize(pml_size, 0);

  a_x_l.resize(pml_size, 0);
  a_x_r.resize(pml_size, 0);
  a_y_b.resize(pml_size, 0);
  a_y_t.resize(pml_size, 0);
  
  be_x_r.resize(pml_size, 0);
  be_x_l.resize(pml_size, 0);
  be_y_t.resize(pml_size, 0);
  be_y_b.resize(pml_size, 0);
  

  ce_x_r.resize(pml_size, 0);
  ce_x_l.resize(pml_size, 0);
  ce_y_t.resize(pml_size, 0);
  ce_y_b.resize(pml_size, 0);

  initPML();

}


// Worker function
void writer_thread_func() {
  while (!finished || !write_queue.empty()) {
    std::unique_lock<std::mutex> lock(queue_mutex);
    cv.wait(lock, [] { return !write_queue.empty() || finished; });

    while (!write_queue.empty()) {
      auto [t, Nx, Ny, Hz_frame, Ex_frame, Ey_frame] = write_queue.front();
      write_queue.pop();
      lock.unlock();

      // Open files in text mode
      std::ofstream Hz_field("../data/Hz/Hz" + std::to_string(t) + ".txt");
      // std::ofstream Ex_field("../data/Ex/Ex" + std::to_string(t) + ".txt");
      std::ofstream Ey_field("../data/Ey/Ey" + std::to_string(t) + ".txt");

      // Write Hz, Ex, Ey in text
      for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
          int idx = i * Ny + j;
          Hz_field << Hz_frame[idx] << " ";
          // Ex_field << Ex_frame[idx] << " ";
          Ey_field << Ey_frame[idx] << " ";
        }
        Hz_field << "\n";
        // Ex_field << "\n";
        Ey_field << "\n";
      }

      Hz_field.close();
      // Ex_field.close();
      Ey_field.close();

      lock.lock();
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

  init(root);

  // Initialize fields (1D contiguous storage)
  vector<double> Ex(Nx * Ny, 0.0);
  vector<double> Ey(Nx * Ny, 0.0);
  vector<double> Hz(Nx * Ny, 0.0);

  // Source position
  int src_i = pml_size + round((source_point[0] - simulation_size[0])/dx);
  int src_j = pml_size + round((source_point[1] - simulation_size[1]) / dy);

  // double h_coeff = dt / (MU_0 * dx);
  // double e_coeff = (dt / (EPSILON_0 * dx));

  double h_coeff = dt / (sqrt(MU_0*EPSILON_0));
  double e_coeff = (dt / (sqrt(MU_0*EPSILON_0)));

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
    Cb_x[i] = (dt / (sqrt(MU_0*EPSILON_0))) / denom;
   
    sig_x<<sigma_x[i]<<endl;
    cax<<Ca_x[i]<<endl;
    cbx<<Cb_x[i]<<endl;    
  }

  for (int j = 0; j < Ny; ++j) {
    double s = sigma_y[j];
    double denom = 1.0 + (s * dt) / (2.0 * EPSILON_0);
    Ca_y[j] = (1.0 - (s * dt) / (2.0 * EPSILON_0)) / denom;
    Cb_y[j] = (dt / (sqrt(MU_0*EPSILON_0))) / denom;        

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

  std::thread writer_thread(writer_thread_func);

// --- FDTD main loop ---
  for(int t = 0; t < N_time_steps; t++) {

    for (int i = 0; i < Nx - 1; i++) {
      for (int j = 0; j < Ny - 1; j++) {
        int idx = i * Ny + j;

        double dEy_dx = inv_kappa_x[i] * (Ey[(i + 1) * Ny + j] - Ey[idx]);
        double dEx_dy = inv_kappa_y[j] * (Ex[i * Ny + j + 1] - Ex[idx]);

        Hz[idx] -= h_coeff * (dEy_dx - dEx_dy);
      }
    }

    // Bottom/top for Ex (y-PML)
    for (int j = 0; j < pml_size; j++) {
      for (int i = 0; i < Nx; i++) {
        int idx = i * Ny + j;
        Psi_Ex_y[idx] = be_y_b[j] * Psi_Ex_y[idx] +
                        ce_y_b[j] * (Hz[idx] - Hz[i * Ny + j - 1]) / dy;
      }
      for (int i = 0; i < Nx; i++) {
        int idx = i * Ny + (Ny - pml_size + j);
        Psi_Ex_y[idx] = be_y_t[j] * Psi_Ex_y[idx] +
                        ce_y_t[j] * (Hz[idx] - Hz[i * Ny + (Ny - pml_size + j - 1)]) / dy;
      }
    }

    // --- Update Ex with CPML ---
    for (int i = 0; i < Nx; i++) {
      for (int j = 1; j < Ny; j++) {
        int idx = i * Ny + j;

        double curlHz = inv_kappa_y[j] * (Hz[idx] - Hz[i * Ny + j - 1]);

        if (j < pml_size)        curlHz += Psi_Ex_y[idx];        // bottom
        if (j >= Ny - pml_size)  curlHz += Psi_Ex_y[idx];        // top

        Ex[idx] = Ca_x[i] * Ex[idx] + Cb_x[i] * curlHz;
      }
    }

    // Left/right for Ey (x-PML)
    for (int i = 0; i < pml_size; i++) {
      for (int j = 0; j < Ny; j++) {
        int idx = i * Ny + j;
        Psi_Ey_x[idx] = be_x_l[i] * Psi_Ey_x[idx] +
                        ce_x_l[i] * (Hz[idx] - Hz[(i - 1) * Ny + j]) / dx;
      }
      for (int j = 0; j < Ny; j++) {
        int idx = (Nx - pml_size + i) * Ny + j;
        Psi_Ey_x[idx] = be_x_r[i] * Psi_Ey_x[idx] +
                        ce_x_r[i] * (Hz[idx] - Hz[(Nx - pml_size + i - 1) * Ny + j]) / dx;
      }
    }

    // --- Update Ey with CPML ---
    for (int i = 1; i < Nx; i++) {
      for (int j = 0; j < Ny; j++) {
        int idx = i * Ny + j;

        double curlHz = inv_kappa_x[i] * (Hz[idx] - Hz[(i - 1) * Ny + j]);

        if (i < pml_size)        curlHz += Psi_Ey_x[idx];        // left
        if (i >= Nx - pml_size)  curlHz += Psi_Ey_x[idx];        // right

        Ey[idx] = Ca_y[j] * Ey[idx] - Cb_y[j] * curlHz;
      }
    }  

    // --- Inject Gaussian hard source into Ey ---
    double time = t * dt;
    // double source = amplitude * exp(-pow((time - pulse_delay)/pulse_width, 2.0));
    double source = sin(omega * time) * exp(-pow((time - t0)/tau, 2));
    // double source = sin(omega * time);
    Ey[src_i*Ny + src_j] += source;

    cout << "Iteration: " << t << endl;

    // --- Capture data asynchronously ---
    if(t % data_capture_interval == 0) {

      while (true) {
        {
          std::lock_guard<std::mutex> lock(queue_mutex);
          if (write_queue.size() <= 10) break;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
      }

      std::lock_guard<std::mutex> lock(queue_mutex);
      std::vector<double> Hz_copy = Hz;
      std::vector<double> Ex_copy = Ex;
      std::vector<double> Ey_copy = Ey;
      write_queue.push({t, Nx, Ny, std::move(Hz_copy), std::move(Ex_copy), std::move(Ey_copy)});
      cv.notify_one();
    }
  }

  {
    std::lock_guard<std::mutex> lock(queue_mutex);
    finished = true;
  }

  cv.notify_one();
  writer_thread.join();

  return 0;
}
