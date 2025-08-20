#pragma once
#include <string>
#include <vector>

class FDTD_1D {
public:
    explicit FDTD_1D(const std::string &control_file);
    ~FDTD_1D() = default;

    void initialize();
    void run();

private:
    // constants
    static constexpr double EPSILON_0 = 8.8541878128e-12;
    static constexpr double MU_0      = 1.256637062e-6;
    static constexpr double c0        = 2.99792458e8;
    static constexpr double imp0      = 376.730313668;

    // parameters (from JSON)
    double simulation_size = 20e-6;
    double frequency = 1.93414489032258e14;
    int    cells_per_wavelength = 100;
    double CFL = 1.0;
    double simulation_time = 1e-12;	

    int    pml_size = 10;
    double m = 3.0;
    double a_max = 0.05;
    double kappa_max = 5.0;

    int    source_position = 10;
    double pulse_width = 1e-14;
    double pulse_delay = 4e-14;
    double amplitude = 1.0;
		double omega0;

    // derived + arrays (same as before) …
    double center_wavelength = 0.0;
    double step_size = 0.0;
    int    N_space_cells = 0;
    double dt = 0.0;
    int    N_time_steps = 0;
    int data_capture_interval = 100;

    std::vector<double> Ex, Ex_prev, Hz, Hz_prev;
    std::vector<double> P_Ex_y, P_Hz_y;
    std::vector<double> eps, sigma, sigma_h, refractive_index;
    std::vector<double> inv_kappa_dy, inv_kappa_h_dy;

    std::vector<double> sigma_pml, sigma_h_pml, kappa_far, kappa_h_far;
    std::vector<double> a_pml, a_h_pml, kappa_left, kappa_h_left;
    std::vector<double> be_y_f, ce_y_f, bh_y_f, ch_y_f;
    std::vector<double> be_y, ce_y, bh_y, ch_y;

    std::vector<double> e_coeff1, e_coeff2, h_coeff1, h_coeff2;

    void load_parameters(const std::string &control_file);
    void setup_grid_and_arrays();
    void setup_pml();
    void setup_coefficients();
    double signal(double time);
    void output_fields(int step);
};
