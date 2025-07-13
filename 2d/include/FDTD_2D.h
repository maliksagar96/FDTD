
#ifndef FDTD_2D_H
#define FDTD_2D_H

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <functional>
#include <jsoncpp/json/value.h>
#include <json/json.h>

using namespace std;

class FDTD_2D {
  public:
    FDTD_2D(string filename);
    ~FDTD_2D();

    void init();

    void set_spacing();

    void set_dt();

    void set_fields();

    void setup_source_function(); 

    void point_source(vector<vector<double>>&, int);
    
    void print_initial_state();

    void fdtd_2d_basic_TMz();

    void fdtd_2d_basic_TMz_PML();
    
    void fdtd_2d_basic_TEz();

    void fdtd_2d_basic_TEz_PML();

    void set_PML_kappa();

    void set_PML_sigma();

    void set_PML_a();

    void set_PML_b();

    void set_PML_c();

    void set_PML_layers();

    void fdtd_2d_PML();

    void write_Matrix_To_File(const vector<vector<double>>&, const string&);

  private:
    const double epsilon_0 = 8.854187817e-12;
    const double mu_0 = 4*M_PI*1e-7;
    const double eta_0 = sqrt(mu_0/epsilon_0);//impedance of free space
    int globalIteration;
    int nx, ny, nt, nc, pml_layers;
    double delta_x, delta_y;
    double frequency_tilde;
    double bandwidth_percent, bandwidth_tilde;
    double tau;
    double dt;
    double epsilon;
    double mu;
    double omega;
    double t0, spread_tilde;
    double source_position[2], domain[4];
    double CFL;
    int data_capture_interval;
    int source_x, source_y;
    string source_type;

    vector<double> kappa_x, kappa_y, pml_a_x, pml_a_y, pml_b_x, pml_b_y, pml_c_x, pml_c_y;
    

    vector<vector<double>> Ca, Cb, ca_ex, cb_ex, ca_ey, cb_ey, da_Hzx, da_Hzy, db_Hzx, db_Hzy;
    vector<vector<double>> Ex, Ey, Ez,Dx, Dy, Dz, Hx, Hy, Hz, epsilon_r, mu_r, sigma;
    vector<vector<double>> psiEz_x, psiEz_y, psiHy_x, psiHy_y, psiEx_x, psiEx_y, psiHx_x, psiHx_y;

    function<double(int)> source_function;

};


#endif