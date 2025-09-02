#ifndef FDTD_2D_H
#define FDTD_2D_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <json/json.h>
#include <functional>

// =========================
// Class Declaration
// =========================
class FDTD_2D {
  public:
    // Constructor / Destructor
    FDTD_2D(const Json::Value&);
    ~FDTD_2D();

    // Initialization
    void read_json(const Json::Value&);
    void compute_grid_parameters();
    void initialize_fields();
    void init_PML();
    void init(const Json::Value&);
  
    // Main solver loop
    void run_TEz();
    void run_TMz();

  private:
    // --- Physical constants ---
    double EPSILON_0;
    double MU_0;
    double c0;

    // --- Domain ---
    double simulation_size[4];   // [x_min, y_min, x_max, y_max]
    double source_point[2];      // source location (x,y)

    // --- Grid ---
    int Nx, Ny;
    int cells_per_wavelength;
    double step_size, dx, dy;

    // --- PML ---
    int pml_size;
    double pml_cond_e;
    // Electric PML parameters
    std::vector<double> sigma_x_r_pml, sigma_x_l_pml, sigma_y_t_pml, sigma_y_b_pml;
    std::vector<double> kappa_x_l, kappa_x_r, kappa_y_b, kappa_y_t;
    std::vector<double> a_x_r, a_x_l, a_y_b, a_y_t;
    std::vector<double> be_x_r, be_x_l, be_y_t, be_y_b;
    std::vector<double> ce_x_r, ce_x_l, ce_y_t, ce_y_b;
    std::vector<double> inv_kappa_x, inv_kappa_y;

    // Magnetic PML parameters
    std::vector<double> sigma_x_r_h, sigma_x_l_h, sigma_y_t_h, sigma_y_b_h;
    std::vector<double> kappa_x_r_h, kappa_x_l_h, kappa_y_t_h, kappa_y_b_h;
    std::vector<double> a_x_r_h, a_x_l_h, a_y_t_h, a_y_b_h;
    std::vector<double> bh_x_r, bh_x_l, bh_y_t, bh_y_b;
    std::vector<double> ch_x_r, ch_x_l, ch_y_t, ch_y_b;
    std::vector<double> inv_kappa_h_x, inv_kappa_h_y;

		double h_coeff, e_coeff;		

    // --- Material arrays ---
    std::vector<double> sigma_x, sigma_y;
    std::vector<double> sigma_x_h, sigma_y_h;

    // --- Time / CFL ---
    double CFL, dt;
    int N_time_steps;
    int data_capture_interval;

    // --- Source ---
    double frequency, center_wavelength, omega;
    double pulse_width, pulse_delay, amplitude;

    // --- CPML fields ---
    // For Ex updates (needs ∂Hz/∂y and ∂Hy/∂z)
    std::vector<double> Psi_Ex_y, Psi_Ex_z;

    // For Ey updates (needs ∂Hx/∂z and ∂Hz/∂x)
    std::vector<double> Psi_Ey_z, Psi_Ey_x;

    // For Ez updates (needs ∂Hy/∂x and ∂Hx/∂y)
    std::vector<double> Psi_Ez_x, Psi_Ez_y;

    // For Hx updates (needs ∂Ey/∂z and ∂Ez/∂y)
    std::vector<double> Psi_Hx_z, Psi_Hx_y;

    // For Hy updates (needs ∂Ez/∂x and ∂Ex/∂z)
    std::vector<double> Psi_Hy_x, Psi_Hy_z;

    // For Hz updates (needs ∂Ex/∂y and ∂Ey/∂x)
    std::vector<double> Psi_Hz_y, Psi_Hz_x;

    // --- Update coefficients ---
    std::vector<double> Ca_x, Ca_y, Ca_z, Cb_x, Cb_y, Cb_z;
    std::vector<double> Da_x, Da_y, Da_z, Db_x, Db_y, Db_z;

    // --- Grid bookkeeping ---
    int vacuum_cells_x, vacuum_cells_y;

    // --- Async I/O ---
		std::queue<std::tuple<int,int,int,
											std::vector<double>,
											std::vector<double>,
											std::vector<double>>> write_queue;

    std::mutex queue_mutex;
    std::condition_variable cv;
    bool finished = false;

    void writer_thread_func();
    
    // --- Fields ---
    std::vector<double> Ex, Ey, Ez, Hx, Hy, Hz;

    // --- Source ---
    int src_i, src_j;
    double source(double);
    std::string source_type;
    std::function<double(double)> source_fn; 
    void setup_source();  
};

#endif // FDTD_2D_H
