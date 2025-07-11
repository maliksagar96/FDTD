#ifndef FDTD_2D_H
#define FDTD_2D_H

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <functional>

using namespace std;

class FDTD_2D {
  public:
    FDTD_2D(string filename);
    ~FDTD_2D();
    
    /**
     * @brief Set the number of time steps for the FDTD simulation.  
     * @param steps The number of time steps to run the simulation.
     */
    void set_time_steps(int steps);

    void set_spacing(int, int);

    void set_excitation(double t0, double spread);
    
    void write_Matrix_To_File(const vector<vector<double>>&, const string&);

    void fdtd_2d_basic_TMz();
    
    void fdtd_2d_basic_TEz();

    void set_PML_layers(int, int);

    void fdtd_2d_PML();

    void setup_source_function(); 

  private:
    int globalIteration;
    int nx, ny, nt, nc;
    double delta_x, delta_y;
    double frequency;
    double bandwidth;
    double tau;
    double dt;
    double epsilon;
    double mu;
    double omega;
    double t0, spread;
    double source_position[2], domain[4];
    string source_type;

    vector<vector<double>> ca_ex, cb_ex, ca_ey, cb_ey, da_Hzx, da_Hzy, db_Hzx, db_Hzy;
    vector<vector<double>> ex, ey, ez, dx, dy, dz, hx, hy, hz, epsilon_r, mu_r, sigma;

    function<double(int)> source_function;

};


#endif