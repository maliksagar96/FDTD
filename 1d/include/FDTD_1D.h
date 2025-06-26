/*
This simulation is done using the first chapter of the book Electromagnetic simulation by Dennis M. Sullivan.
*/


/*



*/

#ifndef FDTD_1D_H
#define FDTD_1D_H

#include <iostream>
#include <vector> 
#include <cmath>
#include <string>

using namespace std;

// Class-based FDTD
class FDTD_1D {
public:
  FDTD_1D();
  ~FDTD_1D();

  ///Set number of iterations.
  void set_time_steps(int nt);

  ///Set spacing in 1D.
  void set_spacing(int nx);

 /**
 * @brief Set Gaussian‑pulse excitation parameters.
 * @param t0     Time‑center of the pulse (in time steps)
 * @param spread Standard deviation controlling pulse width
 */
  void set_excitation(double t0, double spread);
  
  
  /**
  * @brief This is the first function anyone would write to see FDTD in action.
    There are no boundary conditons yet and the waves reflect off the boundaries. This gives good insight to start the FDTD code and is the easiest code one can write.
  */
  void fdtd_1d_basic();

  /**
  *@brief Mur's absorbing boundary condition added. 
  */
  void fdtd_1d_abc();

  /**
  * @brief The simplemost implementation of dielectric material. I will mark some particular area as dielectric material. 
    Right now the dielectric area is hard coded. 
    The reflection from the first dielectric face agress with the interaction with real dielectric medium. 
    The problem lies in the ABC implementation in the dielectric region. 
  */
  void fdtd_1d_dielectric();

  /// Writes 1D array to a text file which can then be visualised.
  void writeArrayToFile(const string, const string, int);

private:
  const double epsilon = 8.854e-12;
  const double mu = 4 * M_PI * 1e-7;
  
  vector<double> ex, hy;
  double t0, spread;
  int nt, globalIteration, nx, nc;
};



#endif
