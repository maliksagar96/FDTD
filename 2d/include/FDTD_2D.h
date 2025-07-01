#ifndef FDTD_2D_H
#define FDTD_2D_H

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

class FDTD_2D {
  public:
    FDTD_2D();
    ~FDTD_2D();
    
    /**
     * @brief Set the number of time steps for the FDTD simulation.  
     * @param steps The number of time steps to run the simulation.
     */
    void set_time_steps(int steps);

    void set_spacing(int, int);

    void set_excitation(double t0, double spread);
    
    void write_Matrix_To_File(const vector<vector<double>>&, const string&);

    void fdtd_2d_basic();

    void set_PML_layers(int, int);

    void fdtd_2d_PML();

  private:
    int globalIteration;
    int nx, ny, nt, nc;
    double dx, dy;
    double dt;
    double epsilon;
    double mu;
    double omega;
    double t0, spread;

    vector<vector<double>> ca_ex, cb_ex, ca_ey, cb_ey, da_Hzx, da_Hzy, db_Hzx, db_Hzy;
    vector<vector<double>> ez, dz, hx, hy, epsilon_r, mu_r, sigma;
};


#endif