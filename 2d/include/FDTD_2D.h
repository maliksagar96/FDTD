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
    

  ///Set number of iterations.
  void set_time_steps(int);

  ///Set spacing in 1D.
  void set_spacing(int, int);

    /**
 * @brief Set Gaussian‑pulse excitation parameters.
 * @param t0     Time‑center of the pulse (in time steps)
 * @param spread Standard deviation controlling pulse width
 */
  void set_excitation(double t0, double spread);
    void fdtd_2d_basic();
    void write_Matrix_To_File(const vector<vector<double>>&, const string&);

    private:
      int globalIteration;
      int nx, ny, nt, nc;
      double dx, dy;
      double dt;
      double epsilon;
      double mu;
      double omega;
      double t0, spread;

      vector<vector<double>> ez, dz, hx, hy, epsilon_r;
};



#endif