#ifndef FDTD_1D_H
#define FDTD_1D_H

#include <iostream>
#include <vector> 
#include <cmath>
#include <string>
#include <functional>

using namespace std;

class FDTD_1D {
public:
  FDTD_1D(string);
  ~FDTD_1D();

  void init();
  void set_spacing();            // keep the one with int
  void set_dt();
  void set_fields();
  void setup_source_function();
  void set_PML_layers();               // added this
  void set_field_coeffs();
  void set_sigma_graded();       //For graded PML
  void print_initial_state(); 


  void fdtd_1d_basic();
  void fdtd_1d_pml();
  void fdtd_1d_graded_pml();
  void fdtd_1d_graded_Cpml();

  void writeArrayToFile(const string, const string, int);
  void write_Vector_To_File(vector<double>&, string);

private:
  const double epsilon_0 = 8.8541878128e-12;
  const double mu_0 = 1.256637062e-6;
  const double c0 = 2.99792458e8;
  const double eta0 = sqrt(mu_0 / epsilon_0);

  double delta_z, CFL, bandwidh_percent, frequency, spread, tau, dt;
  double bandwidth_percent, bandwidth, source_position;
  int pml_layers, source_x, data_capture_interval;
  vector<double> domain = {0,0};
  
  string source_type;

  vector<double> Ex, Hy, psiEx_z, psiHy_z, e1_coeff, e2_coeff, h1_coeff, h2_coeff;
// Electric field PML parameters
vector<double> sigma_z_e, kappa_z_e, pml_a_z_e, pml_b_z_e, pml_c_z_e;

// Magnetic field PML parameters
vector<double> sigma_z_h, kappa_z_h, pml_a_z_h, pml_b_z_h, pml_c_z_h;

  vector<double> mu_r, epsilon_r;

  int nt, nz;
  function<double(int)> source_function;
};

#endif
