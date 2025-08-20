#include <iostream> 
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace std;

double EPSILON_0 = 8.8541878128e-12;
double MU_0 = 1.256637062e-6;
double c0 = 2.99792458e8; 
double imp0 = sqrt(MU_0/EPSILON_0);

double simulation_size = 20e-6;
double step_size = 5e-9;
double N_space_cells = int(simulation_size/step_size);

double dt = step_size/c0;
double simulation_time = 1e-13;
double N_time_steps = int(simulation_time/dt);


vector<double> Ex(N_space_cells, 0), Ex_prev(N_space_cells, 0);
vector<double> Hz(N_space_cells, 0), Hz_prev(N_space_cells, 0);
vector<double> eps(N_space_cells, 1);
vector<double> sigma(N_space_cells, 0);
vector<double> sigma_h(N_space_cells, 0);
vector<double> refractive_index(N_space_cells, 1);

double pml_size = 10;

double signal(double time, double pulse_width, double pulse_delay, double omega0, double amplitude) {
  double envelope = exp(-pow(((time - pulse_delay)/pulse_width), 2));
  double source_signal = envelope * amplitude * sin(omega0*time);
  return source_signal;
}

int main() {

  for(int i = 0;i<N_space_cells;i++) {
   refractive_index[i] = sqrt(eps[i]) ;
  }

  double pml_cond_e = -log(1e-3) * c0 * EPSILON_0 / (2 * refractive_index[N_space_cells-1 - pml_size] * pml_size * step_size);
  
  for(int i = 0;i<pml_size;i++) {
    sigma[N_space_cells-1-pml_size + i] = pml_cond_e;
    sigma_h[N_space_cells-1-pml_size + i] = pml_cond_e * MU_0 / (EPSILON_0 * eps[N_space_cells-1-pml_size + i]);
  }

  vector<double> denominator(N_space_cells, 1);
  vector<double> e_coeff1(N_space_cells, 1);
  vector<double> e_coeff2(N_space_cells, 1);

  vector<double> denominator_h(N_space_cells, 1);
  vector<double> h_coeff1(N_space_cells, 1);
  vector<double> h_coeff2(N_space_cells, 1);

  for(int i=0;i < N_space_cells;i++) {
    denominator[i] = EPSILON_0 * eps[i]/dt + sigma[i]/2;
    e_coeff1[i] = (EPSILON_0 * eps[i]/dt - sigma[i]/2)/denominator[i];
    e_coeff2[i] = 1 / (step_size * denominator[i]);

    denominator_h[i] = MU_0/dt + sigma_h[i]/2;
    h_coeff1[i] = (MU_0/dt - sigma_h[i]/2)/denominator_h[i];
    h_coeff2[i] = 1/(step_size*denominator_h[i]);
  }


  double c = c0/refractive_index[0];
  double c_ = c0/refractive_index.back();

  double a = (c * dt - step_size)/(c * dt + step_size);
  double a_ = (c_ * dt - step_size)/(c_ * dt + step_size);

  //pulse setup
  double center_wavelenght = 1550e-9;
  double omega0 = 2*M_PI*c0/center_wavelenght;
  double pulse_width = 10e-15;
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
      Hz[j] = h_coeff1[j] * Hz_prev[j] + h_coeff2[j] * (Ex[j+1] - Ex[j]);      
    }

    //Magnetic Field Source
    Hz[jsource-1] -=  signal((n+0.5)*dt -t_offset, pulse_width, pulse_delay, omega0, amplitude)/Z;

    //Update Ex
    for(int j = 1;j<N_space_cells;j++) {
      Ex[j] = e_coeff1[j] * Ex_prev[j] + e_coeff2[j] * (Hz[j] - Hz[j-1]);    
    }

    //Electric field source
    Ex[jsource] += signal((n+1)*dt, pulse_width, pulse_delay, omega0, amplitude);


    //Update Boundary Conditions
    Ex[0] = Ex_prev[1] + a * (Ex[1] - Ex_prev[0]);
    

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