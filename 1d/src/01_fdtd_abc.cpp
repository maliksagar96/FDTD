#include <iostream> 
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace std;

double EPSILON_0 = 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;
double MU_0 = 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6;
double c0 = 1/sqrt(EPSILON_0 * MU_0); 
double imp0 = sqrt(MU_0/EPSILON_0);

double simulation_size = 10e-6;
double step_size = 5e-9;
double N_space_cells = int(simulation_size/step_size);

double lambda = 350e-9;
double dt = step_size/c0;
double simulation_time = 1e-12;
double N_time_steps = int(simulation_time/dt);


vector<double> Ex(N_space_cells, 0), Ex_prev(N_space_cells, 0);
vector<double> Hz(N_space_cells, 0), Hz_prev(N_space_cells, 0);
vector<double> eps(N_space_cells, 1);


double signal(double time, double pulse_width, double pulse_delay, double omega0, double amplitude) {
  double envelope = exp(-pow(((time - pulse_delay)/pulse_width), 2));
  double source_signal = envelope * amplitude * sin(omega0*time);
  return source_signal;
}

int main() {

  double h_coeff = dt/(step_size*MU_0);
  vector<double> e_coeff(N_space_cells, 1);  
  vector<double> refractive_index(N_space_cells, 1);

  for(int i = 0;i<N_space_cells;i++) {
    refractive_index[i] = sqrt(eps[i]);
    e_coeff[i] = dt/(step_size*EPSILON_0*eps[i]);
  }

  double c = c0/refractive_index[0];
  double c_ = c0/refractive_index.back();

  double a = (c * dt - step_size)/(c * dt + step_size);
  double a_ = (c_ * dt - step_size)/(c_ * dt + step_size);

  //pulse setup
  double center_wavelenght = 1550e-9;
  double omega0 = 2*M_PI*c0/center_wavelenght;
  double pulse_width = 20e-15;
  double pulse_delay = 4 * pulse_width;

  int jsource = 10;
  double t_offset = refractive_index[jsource] * step_size/(2 * c0);
  double Z = imp0/refractive_index[jsource];
  double amplitude = 1;


  for(int n = 0;n<N_time_steps;n++) {    
    for(int j = 0;j<N_space_cells;j++) {
      Hz_prev[j] = Hz[j];
      Ex_prev[j] = Ex[j];
    }

    //Update Hz
    for(int j = 0;j<N_space_cells-1;j++) {
      Hz[j] = Hz_prev[j] + h_coeff * (Ex[j+1] - Ex[j]);      
    }

    //Magnetic Field Source
    Hz[jsource-1] -=  signal((n+0.5)*dt -t_offset, pulse_width, pulse_delay, omega0, amplitude)/Z;

    //Update Ex
    for(int j = 1;j<N_space_cells;j++) {
      Ex[j] = Ex_prev[j] + e_coeff[j] * (Hz[j] - Hz[j-1]);    
    }

    //Electric field source
    Ex[jsource] += signal((n+1)*dt, pulse_width, pulse_delay, omega0, amplitude);


    //Update Boundary Conditions
    Ex[0] = Ex_prev[1] + a * (Ex[1] - Ex_prev[0]);
    Ex.back() = Ex_prev[Ex_prev.size()-2] + a * (Ex[Ex.size() - 2] - Ex_prev.back());

    if(n%100 == 0) {
      double max_val = *std::max_element(Ex.begin(), Ex.end());
      double min_val = *std::min_element(Ex.begin(), Ex.end());
      cout<<"Time step : "<<n<<endl;
      cout<<"Emin = "<<min_val<<", Emax = "<<max_val<<endl;
      std::ofstream ex("data/Ex/ElectricField_" + to_string(n)  + ".txt");
      for(int i = 0;i<N_space_cells;i++) {
        ex<<Ex[i]<<endl;
      }
      ex.close();

    } 
  }

  return 0;
}