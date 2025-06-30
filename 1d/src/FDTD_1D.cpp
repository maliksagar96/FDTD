#include <iostream> 
#include <vector> 
#include <cmath>
#include <string>
#include <fstream>
#include "FDTD_1D.h"

using namespace std;

// Constructor
FDTD_1D::FDTD_1D() : globalIteration(0), nt(0), nx(0), nc(0), t0(0.0), spread(1.0) {}

// Destructor
FDTD_1D::~FDTD_1D() {}

// Set number of iterations.
void FDTD_1D::set_time_steps(int nt) {
  this->nt = nt;
}

// Set spatial size and initialize fields
void FDTD_1D::set_spacing(int nx) {
  this->nx = nx;
  this->nc = nx / 2;
  ex.resize(nx, 0.0);
  hy.resize(nx, 0.0);
}

// Set excitation parameters. The center of the gaussian in time and the std dev.
void FDTD_1D::set_excitation(double t0, double spread) {
  this->t0 = t0;
  this->spread = spread;
}

// Write fields to text file
void FDTD_1D::writeArrayToFile(const string field, const string data_directory, int t) {
  string filename = data_directory + field + "_" + to_string(t) + ".txt";
  ofstream outFile(filename);
  if (!outFile) {
    cerr << "Error opening file: " << filename << "\n";
    return;
  }

  if (field == "ElectricField") {
    for (double val : ex) outFile << val << "\n";
  } else if (field == "MagneticField") {
    for (double val : hy) outFile << val << "\n";
  }

  outFile.close();
}

//No boundary condition applied. 
void FDTD_1D::fdtd_1d_basic() {
  for (int t = 0; t < nt; t++) {    

    /// If we keep the spacing as .5 delta x/c then we reduce the discretised maxwell equations in the following form. 
    /// read the first chapter of the Sullivan book. 
    for (int i = 1; i < nx; i++) {
      ex[i] += 0.5 * (hy[i - 1] - hy[i]);
    }

    ex[nc] = exp(-0.5 * pow((t0 - t) / spread, 2));

    for (int i = 0; i < nx - 1; i++) {
      hy[i] += 0.5 * (ex[i] - ex[i + 1]);
    }
    
    writeArrayToFile("ElectricField","../data/Ex/", t);
    writeArrayToFile("MagneticField","../data/Hy/", t);
  }
}

/**
 @brief : 1d FDTD code with boundary condition applied. In plain and simple language the boundary conditons is such that it takes 2 time steps for ligth to cross one simulation cell. 
 The 0th and the nth cell are updated every 2nd time step.  
*/
void FDTD_1D::fdtd_1d_abc() {

  // Taking 2 elements one is the 0th element and the other is nth element.
  vector<double> ex_boundary_even(2, 0);
  vector<double> ex_boundary_odd(2, 0);

  double ex_low_m1, ex_low_m2, ex_high_m1, ex_high_m2;
  
  for (int t = 0; t < nt; t++) {
    
    for (int i = 1; i < nx-1; i++) {
      ex[i] += 0.5 * (hy[i - 1] - hy[i]);
    }

    ex[nc] = exp(-0.5 * pow((t0 - t) / spread, 2));

    if((t > 0) && (t%2 == 0)) {
      ex[0] =    ex_boundary_even[0];
      ex[nx-1] = ex_boundary_even[1];
    }

    if((t > 1) && (t%2 == 1)) {
      ex[0] =    ex_boundary_odd[0];
      ex[nx-1] = ex_boundary_odd[1];
    }

    // // Storing the value of boundary points every 2nd iteration. This is the condition in the Sullivan book. I have written another code which is more readable and gives the same affect. 
    // ex[0] = ex_low_m2;
    // ex_low_m2 = ex_low_m1;
    // ex_low_m1 = ex[1];
    // ex[nx-1] = ex_high_m2;
    // ex_high_m2 = ex_high_m1;
    // ex_high_m1 = ex[nx-2];

    if((t >=0) && (t%2 == 0)) {
      ex_boundary_even[0] = ex[1];
      ex_boundary_even[1] = ex[nx-2];
    }
    
    if((t >=1) && (t%2 == 1)) {
      ex_boundary_odd[0] = ex[1];
      ex_boundary_odd[1] = ex[nx-2];
    }

    for (int i = 0; i < nx - 1; i++) {
      hy[i] += 0.5 * (ex[i] - ex[i + 1]);
    }

    //wrting to text file. 
    writeArrayToFile("ElectricField","../data/Ex/", t);
    writeArrayToFile("MagneticField","../data/Hy/", t);
  }

}

void FDTD_1D::fdtd_1d_dielectric() {
 
  // Taking 2 elements one is the 0th element and the other is nth element.
  vector<double> ex_boundary_even(2, 0);
  vector<double> ex_boundary_odd(2, 0);

  vector<double> dielectric_constant(nx, 1);

  for(int i = 175;i<200;i++) {
    dielectric_constant[i] = 0.2;
  }
  
  for (int t = 0; t < nt; t++) {
    
    for (int i = 1; i < nx; i++) {
      ex[i] += 0.5 * dielectric_constant[i] * (hy[i - 1] - hy[i]);
    }

    ex[nc] = exp(-0.5 * pow((t0 - t) / spread, 2));

    if((t > 0) && (t%2 == 0)) {
      ex[0] =    ex_boundary_even[0];
      ex[nx-1] = ex_boundary_even[1];
    }

    if((t > 1) && (t%2 == 1)) {
      ex[0] =    ex_boundary_odd[0];
      ex[nx-1] = ex_boundary_odd[1];
    }

    if((t >=0) && (t%2 == 0)) {
      ex_boundary_even[0] = ex[1];
      ex_boundary_even[1] = ex[nx-2];
    }
    
    if((t >=1) && (t%2 == 1)) {
      ex_boundary_odd[0] = ex[1];
      ex_boundary_odd[1] = ex[nx-2];
    }

    for (int i = 0; i < nx - 1; i++) {
      hy[i] += 0.5 * (ex[i] - ex[i + 1]);
    }

    //wrting to text file. 
    writeArrayToFile("ElectricField","../data/Ex/", t);
    writeArrayToFile("MagneticField","../data/Hy/", t);
  }

  
}