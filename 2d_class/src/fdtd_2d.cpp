#include <fdtd_2d.h>
#include <publish.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <json/json.h>
#include <sstream>     


FDTD_2D::FDTD_2D(const Json::Value& root) {
  read_json(root);
  compute_grid_parameters();
  initialize_fields();
  init_PML();  
  setup_source(); 
}

FDTD_2D::~FDTD_2D() = default;


void FDTD_2D::read_json(const Json::Value& root) {

  EPSILON_0 = 8.8541878128e-12;
  MU_0      = 1.256637062e-6;
  c0        = 2.99792458e8;

  frequency = root["frequency"].asDouble();
  N_time_steps = root["Iterations"].asInt();
  data_capture_interval = root["data_capture_interval"].asInt();
  pml_size = root["pml_size"].asInt();  
  source_type = root["source_type"].asString();

  simulation_size[0] = root["simulation_size"][0].asDouble();
  simulation_size[1] = root["simulation_size"][1].asDouble();
  simulation_size[2] = root["simulation_size"][2].asDouble();
  simulation_size[3] = root["simulation_size"][3].asDouble();
  cells_per_wavelength = root["cells_per_wavelength"].asInt();
  amplitude = root["amplitude"].asDouble();
  CFL = root["CFL"].asDouble();
  source_point[0] = root["source_point"][0].asDouble();
  source_point[1] = root["source_point"][1].asDouble();
  std::cout<<"Source_point = "<<source_point[0]<<", source_point[1] = "<<source_point[1]<<std::endl;
  std::cout<<"PML_size = "<<pml_size<<std::endl;
  

  if((source_point[0] < simulation_size[0] || source_point[0] > simulation_size[2] || 
      source_point[1] < simulation_size[1] || source_point[1] > simulation_size[3])) {
    std::cerr<<"Source point not inside the simualtion domain."<<std::endl;
    exit(1);
  }

  std::cout<<"************JSON Parameters***********************"<<std::endl;
  std::cout<<"Frequency = "<<frequency<<std::endl;
  std::cout<<"Iterations = "<<N_time_steps<<std::endl;
  std::cout<<"data_capture_interval = "<<data_capture_interval<<std::endl;
  std::cout<<"simulation_size = "<<simulation_size[0]<<", "<<simulation_size[1]<<", "<<simulation_size[2]<<", "<<simulation_size[3]<<std::endl;
  std::cout<<"source_point = "<<source_point[0]<<", "<<source_point[1]<<std::endl;
  std::cout<<"cells_per_wavelength = "<<cells_per_wavelength<<std::endl;
  std::cout<<"amplitude = "<<amplitude<<std::endl;
  std::cout<<"CFL = "<<CFL<<std::endl;
  std::cout<<"pml_size = "<<pml_size<<std::endl;
  std::cout<<"Reading JSON file completed."<<std::endl;
}

void FDTD_2D::compute_grid_parameters() {

  center_wavelength = c0/frequency;      
  omega = 2 * M_PI * frequency;
  step_size = center_wavelength / cells_per_wavelength;
  dx = step_size;
  dy = step_size;
  dt = CFL / (c0 * sqrt((1.0/(dx*dx)) + (1.0/(dy*dy))));
  vacuum_cells_x = static_cast<int>(std::fabs(simulation_size[2] - simulation_size[0]) / dx);
  vacuum_cells_y = static_cast<int>(std::fabs(simulation_size[3] - simulation_size[1]) / dy);

  pulse_width = dt/(2*M_PI*0.02);
  pulse_delay = 4*pulse_width;
  
  Nx = vacuum_cells_x + 2*pml_size;
  Ny = vacuum_cells_y + 2*pml_size;

  src_i = pml_size + int((source_point[0] - simulation_size[0])/dx);
  src_j = pml_size + int((source_point[1] - simulation_size[1])/dy);

  h_coeff = dt / (sqrt(MU_0*EPSILON_0));
  e_coeff = dt / (sqrt(MU_0*EPSILON_0));
  std::cout<<"Completed grid computation."<<std::endl;
}

void FDTD_2D::initialize_fields() {
  // --- Main fields ---
  Ex.resize(Nx * Ny, 0.0);
  Ey.resize(Nx * Ny, 0.0);
  Ez.resize(Nx * Ny, 0.0);

  Hx.resize(Nx * Ny, 0.0);
  Hy.resize(Nx * Ny, 0.0);
  Hz.resize(Nx * Ny, 0.0);

  // --- Conductivities (electric + magnetic) ---
  sigma_x.resize(Nx, 0.0); sigma_y.resize(Ny, 0.0);
  sigma_x_h.resize(Nx, 0.0); sigma_y_h.resize(Ny, 0.0);

  // --- Coefficients for E ---
  Ca_x.resize(Nx, 1.0); Cb_x.resize(Nx, dt / EPSILON_0);
  Ca_y.resize(Ny, 1.0); Cb_y.resize(Ny, dt / EPSILON_0);

  // --- Coefficients for H ---
  Da_x.resize(Nx, 1.0); Db_x.resize(Nx, dt / MU_0);
  Da_y.resize(Ny, 1.0); Db_y.resize(Ny, dt / MU_0);

  // --- Inverse kappa ---
  inv_kappa_x.resize(Nx, 1.0 / dx);
  inv_kappa_y.resize(Ny, 1.0 / dy);
  inv_kappa_h_x.resize(Nx, 1.0 / dx);
  inv_kappa_h_y.resize(Ny, 1.0 / dy);

  // --- Psi auxiliary fields ---
  Psi_Ex_y.resize(Nx * Ny, 0.0);
  Psi_Ey_x.resize(Nx * Ny, 0.0);
  Psi_Hx_y.resize(Nx * Ny, 0.0);
  Psi_Hy_x.resize(Nx * Ny, 0.0);

  // --- PML arrays (electric) ---
  sigma_x_r_pml.resize(pml_size, 0.0); sigma_x_l_pml.resize(pml_size, 0.0);
  sigma_y_t_pml.resize(pml_size, 0.0); sigma_y_b_pml.resize(pml_size, 0.0);

  kappa_x_r.resize(pml_size, 0.0); kappa_x_l.resize(pml_size, 0.0);
  kappa_y_t.resize(pml_size, 0.0); kappa_y_b.resize(pml_size, 0.0);

  a_x_l.resize(pml_size, 0.0); a_x_r.resize(pml_size, 0.0);
  a_y_b.resize(pml_size, 0.0); a_y_t.resize(pml_size, 0.0);

  be_x_r.resize(pml_size, 0.0); be_x_l.resize(pml_size, 0.0);
  be_y_t.resize(pml_size, 0.0); be_y_b.resize(pml_size, 0.0);

  ce_x_r.resize(pml_size, 0.0); ce_x_l.resize(pml_size, 0.0);
  ce_y_t.resize(pml_size, 0.0); ce_y_b.resize(pml_size, 0.0);

  // --- PML arrays (magnetic) ---
  sigma_x_r_h.resize(pml_size, 0.0); sigma_x_l_h.resize(pml_size, 0.0);
  sigma_y_t_h.resize(pml_size, 0.0); sigma_y_b_h.resize(pml_size, 0.0);

  kappa_x_r_h.resize(pml_size, 0.0); kappa_x_l_h.resize(pml_size, 0.0);
  kappa_y_t_h.resize(pml_size, 0.0); kappa_y_b_h.resize(pml_size, 0.0);

  a_x_r_h.resize(pml_size, 0.0); a_x_l_h.resize(pml_size, 0.0);
  a_y_t_h.resize(pml_size, 0.0); a_y_b_h.resize(pml_size, 0.0);

  bh_x_r.resize(pml_size, 0.0); bh_x_l.resize(pml_size, 0.0);
  bh_y_t.resize(pml_size, 0.0); bh_y_b.resize(pml_size, 0.0);

  ch_x_r.resize(pml_size, 0.0); ch_x_l.resize(pml_size, 0.0);
  ch_y_t.resize(pml_size, 0.0); ch_y_b.resize(pml_size, 0.0);

  std::cout << "Initialised Fields" << std::endl;
}


void FDTD_2D::init_PML() {

  double m = 3;
  double kappa_max = 5.0;  
  double a_max     = 0.5;  

  pml_cond_e = -(m + 1) * log(1e-9) * c0 * EPSILON_0 / (2 * pml_size * dx);

  // ---------- Right (+x) PML ----------
  for (int i = 0; i < pml_size; i++) {
    double x  = static_cast<double>(i) / pml_size;
    double xh = (i + 0.5) / pml_size; // staggered for H

    // Electric
    sigma_x_r_pml[i] = pml_cond_e * pow(x, m);    
    kappa_x_r[i]     = 1.0 + (kappa_max - 1.0) * pow(x, m);
    a_x_r[i]         = a_max * ((pml_size - i) / static_cast<double>(pml_size));
    be_x_r[i]        = exp(-(sigma_x_r_pml[i] / kappa_x_r[i] + a_x_r[i]) * dt / EPSILON_0);
    ce_x_r[i]        = sigma_x_r_pml[i] * (be_x_r[i] - 1.0) /
                      ((sigma_x_r_pml[i] + kappa_x_r[i] * a_x_r[i]) * kappa_x_r[i]);
    inv_kappa_x[Nx - pml_size + i] = 1.0 / (kappa_x_r[i] * dx);

    // Magnetic
    sigma_x_r_h[i] = pml_cond_e * pow(xh, m);    
    kappa_x_r_h[i] = 1.0 + (kappa_max - 1.0) * pow(xh, m);
    a_x_r_h[i]     = a_max * ((pml_size - i + 0.5) / static_cast<double>(pml_size));
    bh_x_r[i]      = exp(-(sigma_x_r_h[i] / kappa_x_r_h[i] + a_x_r_h[i]) * dt / EPSILON_0);
    ch_x_r[i]      = sigma_x_r_h[i] * (bh_x_r[i] - 1.0) /
                    ((sigma_x_r_h[i] + kappa_x_r_h[i] * a_x_r_h[i]) * kappa_x_r_h[i]);
    inv_kappa_h_x[Nx - pml_size + i] = 1.0 / (kappa_x_r_h[i] * dx);
  }

  // ---------- Left (−x) PML ----------
  for (int i = 0; i < pml_size; i++) {
    double x  = static_cast<double>(pml_size - i) / pml_size;
    double xh = (pml_size - i - 0.5) / pml_size;

    // Electric
    sigma_x_l_pml[i] = pml_cond_e * pow(x, m);    
    kappa_x_l[i]     = 1.0 + (kappa_max - 1.0) * pow(x, m);
    a_x_l[i]         = a_max * (static_cast<double>(i) / pml_size);
    be_x_l[i]        = exp(-(sigma_x_l_pml[i] / kappa_x_l[i] + a_x_l[i]) * dt / EPSILON_0);
    ce_x_l[i]        = sigma_x_l_pml[i] * (be_x_l[i] - 1.0) /
                      ((sigma_x_l_pml[i] + kappa_x_l[i] * a_x_l[i]) * kappa_x_l[i]);
    inv_kappa_x[i]   = 1.0 / (kappa_x_l[i] * dx);

    // Magnetic
    sigma_x_l_h[i] = pml_cond_e *  pow(xh, m);    
    kappa_x_l_h[i] = 1.0 + (kappa_max - 1.0) * pow(xh, m);
    a_x_l_h[i]     = a_max * ((i + 0.5) / pml_size);
    bh_x_l[i]      = exp(-(sigma_x_l_h[i] / kappa_x_l_h[i] + a_x_l_h[i]) * dt / EPSILON_0);
    ch_x_l[i]      = sigma_x_l_h[i] * (bh_x_l[i] - 1.0) /
                    ((sigma_x_l_h[i] + kappa_x_l_h[i] * a_x_l_h[i]) * kappa_x_l_h[i]);
    inv_kappa_h_x[i] = 1.0 / (kappa_x_l_h[i] * dx);
  }

  // ---------- Top (+y) PML ----------
  for (int j = 0; j < pml_size; j++) {
    double y  = static_cast<double>(j) / pml_size;
    double yh = (j + 0.5) / pml_size;

    // Electric
    sigma_y_t_pml[j] = pml_cond_e * pow(y, m);    
    kappa_y_t[j]     = 1.0 + (kappa_max - 1.0) * pow(y, m);
    a_y_t[j]         = a_max * ((pml_size - j) / static_cast<double>(pml_size));
    be_y_t[j]        = exp(-(sigma_y_t_pml[j] / kappa_y_t[j] + a_y_t[j]) * dt / EPSILON_0);
    ce_y_t[j]        = sigma_y_t_pml[j] * (be_y_t[j] - 1.0) /
                      ((sigma_y_t_pml[j] + kappa_y_t[j] * a_y_t[j]) * kappa_y_t[j]);
    inv_kappa_y[Ny - pml_size + j] = 1.0 / (kappa_y_t[j] * dy);

    // Magnetic
    sigma_y_t_h[j] = pml_cond_e * pow(yh, m);    
    kappa_y_t_h[j] = 1.0 + (kappa_max - 1.0) * pow(yh, m);
    a_y_t_h[j]     = a_max * ((pml_size - j + 0.5) / static_cast<double>(pml_size));
    bh_y_t[j]      = exp(-(sigma_y_t_h[j] / kappa_y_t_h[j] + a_y_t_h[j]) * dt / EPSILON_0);
    ch_y_t[j]      = sigma_y_t_h[j] * (bh_y_t[j] - 1.0) /
                    ((sigma_y_t_h[j] + kappa_y_t_h[j] * a_y_t_h[j]) * kappa_y_t_h[j]);
    inv_kappa_h_y[Ny - pml_size + j] = 1.0 / (kappa_y_t_h[j] * dy);
  }

  // ---------- Bottom (−y) PML ----------
  for (int j = 0; j < pml_size; j++) {
    double y  = static_cast<double>(pml_size - j) / pml_size;
    double yh = (pml_size - j - 0.5) / pml_size;

    // Electric
    sigma_y_b_pml[j] = pml_cond_e * pow(y, m);    
    kappa_y_b[j]     = 1.0 + (kappa_max - 1.0) * pow(y, m);
    a_y_b[j]         = a_max * (static_cast<double>(j) / pml_size);
    be_y_b[j]        = exp(-(sigma_y_b_pml[j] / kappa_y_b[j] + a_y_b[j]) * dt / EPSILON_0);
    ce_y_b[j]        = sigma_y_b_pml[j] * (be_y_b[j] - 1.0) /
                      ((sigma_y_b_pml[j] + kappa_y_b[j] * a_y_b[j]) * kappa_y_b[j]);
    inv_kappa_y[j]   = 1.0 / (kappa_y_b[j] * dy);

    // Magnetic
    sigma_y_b_h[j] = pml_cond_e * (MU_0/EPSILON_0) * pow(yh, m);    
    kappa_y_b_h[j] = 1.0 + (kappa_max - 1.0) * pow(yh, m);
    a_y_b_h[j]     = a_max * ((j + 0.5) / pml_size);
    bh_y_b[j]      = exp(-(sigma_y_b_h[j] / kappa_y_b_h[j] + a_y_b_h[j]) * dt / EPSILON_0);
    ch_y_b[j]      = sigma_y_b_h[j] * (bh_y_b[j] - 1.0) /
                    ((sigma_y_b_h[j] + kappa_y_b_h[j] * a_y_b_h[j]) * kappa_y_b_h[j]);
    inv_kappa_h_y[j] = 1.0 / (kappa_y_b_h[j] * dy);
  }

  for (int i = 0; i < Nx; ++i) {
    double s = sigma_x[i];
    double denom = 1.0 + (s * dt) / (2.0 * EPSILON_0);
    Ca_x[i] = (1.0 - (s * dt) / (2.0 * EPSILON_0)) / denom;
    Cb_x[i] = (dt / (sqrt(MU_0*EPSILON_0))) / denom;   
  }

  for (int j = 0; j < Ny; ++j) {
    double s = sigma_y[j];
    double denom = 1.0 + (s * dt) / (2.0 * EPSILON_0);
    Ca_y[j] = (1.0 - (s * dt) / (2.0 * EPSILON_0)) / denom;
    Cb_y[j] = (dt / (sqrt(MU_0*EPSILON_0))) / denom;        
  }

  for (int i = 0; i < Nx; ++i) {
    double s = sigma_x_h[i];  // magnetic conductivity in x-PML
    double denom = 1.0 + (s * dt) / (2.0 * MU_0);
    Da_x[i] = (1.0 - (s * dt) / (2.0 * MU_0)) / denom;
    Db_x[i] = (dt / (sqrt(MU_0*EPSILON_0))) / denom;
  }

  for (int j = 0; j < Ny; ++j) {
    double s = sigma_y_h[j];  // magnetic conductivity in y-PML
    double denom = 1.0 + (s * dt) / (2.0 * MU_0);
    Da_y[j] = (1.0 - (s * dt) / (2.0 * MU_0)) / denom;
    Db_y[j] = (dt / (sqrt(MU_0*EPSILON_0))) / denom;
  }

  std::cout<<"Initialising PML"<<std::endl;

}

void FDTD_2D::writer_thread_func() {
  while (!finished || !write_queue.empty()) {
    std::unique_lock<std::mutex> lock(queue_mutex);
    cv.wait(lock, [this] { return !write_queue.empty() || finished; });

    while (!write_queue.empty()) {
      auto [t, Nx, Ny, Hz_frame, Ex_frame, Ey_frame] = write_queue.front();
      write_queue.pop();
      lock.unlock();

      std::ofstream Hz_field("data/Hz/Hz" + std::to_string(t) + ".txt");
      std::ofstream Ey_field("data/Ey/Ey" + std::to_string(t) + ".txt");

      for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
          int idx = i * Ny + j;
          Hz_field << Hz_frame[idx] << " ";
          Ey_field << Ey_frame[idx] << " ";
        }
        Hz_field << "\n";
        Ey_field << "\n";
      }

      Hz_field.close();
      Ey_field.close();

      lock.lock();
    }
  }
}

void FDTD_2D::setup_source() {
    if (source_type == "sin") {
        source_fn = [=](double t){ return amplitude * sin(omega * t); };
    }
    else if (source_type == "gaussian") {
        source_fn = [=](double t){ return amplitude * exp(-pow((t - pulse_delay)/pulse_width, 2.0)); };
    }
    else if (source_type == "gaussian_sin") {
        source_fn = [=](double t){ return amplitude * sin(omega * t) * exp(-pow((t - pulse_delay)/pulse_width, 2.0)); };
    }
}

double FDTD_2D::source(double time) {
    return source_fn(time);
}

void FDTD_2D::run_TEz() {

  std::thread writer_thread(&FDTD_2D::writer_thread_func, this);

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
    for (int j = 1; j < pml_size; j++) {
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
        Ex[idx] = Ca_x[i] * Ex[idx] + Cb_x[i] * (curlHz + Psi_Ex_y[idx]);
      }
    }

    // Left/right for Ey (x-PML)
    for (int i = 1; i < pml_size; i++) {
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
        Ey[idx] = Ca_y[j] * Ey[idx] - Cb_y[j] * (curlHz + Psi_Ey_x[idx]);
      }
    }  

    // --- Inject Gaussian hard source into Ey ---
    double time = t * dt;  
    Hz[src_i*Ny + src_j] += source(time);
    
    std::cout << "Iteration: " << t << std::endl;

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

}

