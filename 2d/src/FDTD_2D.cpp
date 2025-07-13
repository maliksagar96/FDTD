/*
In the FDTD formulation what I am taking is the normalised Electric Field E~ = sqrt(epsilon_0/mu_0) * E. 
When we do this assumption every equation gets divided by the speed of light, which is 1/sqrt(epsilon_0*mu_0).
And then I am taking t~ = c*t. This assumption removes all the c's from the equations.

The equations are:
1. Hx(i,j) = Hx(i,j) + delta_t/delta_x*(Ey(i,j) - Ey(i+1,j))
2. Hy(i,j) = Hy(i,j) + delta_t/delta_y*(Ex(i,j) - Ex(i,j+1))
3. Ez(i,j) = Ez(i,j) + delta_t/delta_y*(Hy(i,j) - Hy(i,j+1))
4. Ex(i,j) = Ex(i,j) + delta_t/delta_x*(Hz(i,j) - Hz(i,j+1))
5. Ey(i,j) = Ey(i,j) + delta_t/delta_x*(Hz(i,j) - Hz(i+1,j))

With the above formulation the mesh size will still remain same. And we have to mesh our object according to the frequency of interest.
The mesh size is taken to be 1/10th of the wavelength. 

lambda = c/f
mesh size = lambda/10 = c/(10*f) = 1/10*f~

Since t~ = c*t, 
f~ = f/c

delta_t < c * sqrt(2) * delta_x
which leadas to 
delta_t_tilde < sqrt(2) * delta_x

So I will take delta_t_tilde =  delta_x/sqrt(2) if CFL number is 1.

I will be assuming the length of the PML layers to be equal to the wavelength.
So, number of PML layers = wavelength/delta_x.

*/

#include <vector> 
#include <cmath>
#include <string>
#include <fstream>
#include <jsoncpp/json/value.h>
#include <json/json.h>
#include <functional>

#include "FDTD_2D.h"


using namespace std;

FDTD_2D::FDTD_2D(string filename) {
  ifstream file(filename);
  Json::Value jroot;
  Json::CharReaderBuilder builder;
  string errs;

  Json::parseFromStream(builder, file, &jroot, &errs);

  //Setting Frequency
  frequency_tilde = jroot["frequency"].asDouble()/3e8; // See the comments in the header file to see the definition of frequency_tilde. 
  CFL = jroot["CFL"].asDouble();

  //Setting Domain size in x and y direction
  domain[0] = jroot["domain_x"][0].asDouble();
  domain[1] = jroot["domain_x"][1].asDouble();
  domain[2] = jroot["domain_y"][0].asDouble();
  domain[3] = jroot["domain_y"][1].asDouble();
  
  //Setting Bandwidth.
  bandwidth_percent = jroot["bandwidth"].asInt();

  //Setting Source Position
  source_position[0] = jroot["source_position"][0].asDouble();
  source_position[1] = jroot["source_position"][1].asDouble();
  

  if(source_position[0] > domain[1] || source_position[0] < domain[0] || source_position[1] > domain[3] || source_position[1] < domain[2]){
    cerr << "Source position is out of domain" << endl;
    exit(1);
  }

  //Setting Number of Iterations
  nt = jroot["iterations"].asInt();  
  data_capture_interval = jroot["data_capture_interval"].asInt();
  
  //Setting Source Type (gaussian or sin)
  source_type = jroot["source"].asString();
}

FDTD_2D::~FDTD_2D() = default;

void FDTD_2D::set_spacing(){
  delta_x = 1/(10*frequency_tilde);
  delta_y = delta_x;
  double wavelength = 1/frequency_tilde;
  pml_layers = int(wavelength/delta_x);

  nx = int((domain[1] - domain[0])/delta_x) + 2*pml_layers;
  ny = int((domain[3] - domain[2])/delta_y) + 2*pml_layers;
  cout<<"nx: "<<nx<<" ny: "<<ny<<endl;
}

void FDTD_2D::set_dt(){
  dt = CFL*delta_x/sqrt(2);  
  bandwidth_tilde = frequency_tilde*bandwidth_percent/100;
  spread_tilde = 1 / (bandwidth_tilde*2*M_PI);
  tau = 3*spread_tilde;
}

void FDTD_2D::setup_source_function() {
  source_x = pml_layers + int((domain[1] - source_position[0])/delta_x);
  source_y = pml_layers + int((domain[3] - source_position[1])/delta_y);
  if (source_type == "gaussian") {
      source_function = [this](int t) {
          return exp(-0.5 * pow((t * dt - tau) / spread_tilde, 2));
      };
  } else if (source_type == "sin") {
      source_function = [this](int t) {
          return sin(2 * M_PI * frequency_tilde * t * dt);
      };
  } else {
      cerr << "Unknown source type: " << source_type << endl;
      exit(1);
  }
}

void FDTD_2D::point_source(vector<vector<double>>& Ez, int t){  
  Ez[source_x][source_y] = source_function(t);  
}

void FDTD_2D::set_fields(){
  Ex.resize(nx, vector<double>(ny, 0.0));
  Ey.resize(nx, vector<double>(ny, 0.0));
  Ez.resize(nx, vector<double>(ny, 0.0));
  Dx.resize(nx, vector<double>(ny, 0.0));
  Dy.resize(nx, vector<double>(ny, 0.0));
  Dz.resize(nx, vector<double>(ny, 0.0));
  Hx.resize(nx, vector<double>(ny, 0.0));
  Hy.resize(nx, vector<double>(ny, 0.0));
  Hz.resize(nx, vector<double>(ny, 0.0));
  epsilon_r.resize(nx, vector<double>(ny, 1.0));
  mu_r.resize(nx, vector<double>(ny, 1.0));
  sigma.resize(nx, vector<double>(ny, 0.0));
  kappa_x.resize(nx, 1.0);
  kappa_y.resize(ny, 1.0);
  pml_a_x.resize(nx, 0.0);
  pml_a_y.resize(ny, 0.0);
  pml_b_x.resize(nx, 1.0);
  pml_b_y.resize(ny, 1.0);
  pml_c_x.resize(nx, 0.0);
  pml_c_y.resize(ny, 0.0);
  psiEz_x.resize(nx, vector<double>(ny, 0.0));
  psiEz_y.resize(nx, vector<double>(ny, 0.0));
  psiHy_x.resize(nx, vector<double>(ny, 0.0));
  psiHy_y.resize(nx, vector<double>(ny, 0.0));
  psiEx_x.resize(nx, vector<double>(ny, 0.0));
  psiEx_y.resize(nx, vector<double>(ny, 0.0));
  psiHx_x.resize(nx, vector<double>(ny, 0.0));
  psiHx_y.resize(nx, vector<double>(ny, 0.0));

}

void FDTD_2D::init() {
  set_spacing();
  set_dt();
  set_fields();
  setup_source_function();
  set_PML_layers();
}

void FDTD_2D::print_initial_state(){
  cout<<"delta_x: "<<delta_x<<" delta_y: "<<delta_y<<endl;
  cout<<"frequency_tilde: "<<frequency_tilde<<endl;
  cout<<"CFL: "<<CFL<<endl;
  cout<<"dt: "<<dt<<endl;
  cout<<"nx: "<<nx<<" ny: "<<ny<<endl;
  cout<<"source_type: "<<source_type<<endl;
  cout<<"source_position: "<<source_position[0]<<" "<<source_position[1]<<endl;
  cout<<"bandwidth_percent: "<<bandwidth_percent<<endl;
  cout<<"bandwidth_tilde: "<<bandwidth_tilde<<endl;
  cout<<"spread_tilde: "<<spread_tilde<<endl;
  cout<<"tau: "<<tau<<endl;
  cout<<"nt: "<<nt<<endl;
  cout<<"pml_layers: "<<pml_layers<<endl;
}

void FDTD_2D::set_PML_sigma() {
  double grading_order = 3.0;
  double sigma_max = 1.0;

  // Left boundary
  for (int i = 0; i < pml_layers; i++) {
    double sigma_x = sigma_max * pow((pml_layers - i - 1) / (double)pml_layers, grading_order);
    for (int j = 0; j < ny; j++) {
      sigma[i][j] += sigma_x;
    }
  }

  // Bottom boundary
  for (int j = 0; j < pml_layers; j++) {
    double sigma_y = sigma_max * pow((pml_layers - j - 1) / (double)pml_layers, grading_order);
    for (int i = 0; i < nx; i++) {
      sigma[i][j] += sigma_y;
    }
  }

  // Right boundary
  for(int i = nx-pml_layers;i<nx;i++){
    double sigma_x = sigma_max * pow((i - (nx - pml_layers - 1)) / (double)pml_layers, grading_order);  
    for(int j = 0;j<ny;j++){
      sigma[i][j] += sigma_x;
    }
  }

  // Top boundary
  for(int j = ny-pml_layers;j<ny;j++){
    double sigma_y = sigma_max * pow((j - (ny - pml_layers - 1)) / (double)pml_layers, grading_order);  
    for(int i = 0;i<nx;i++){
      sigma[i][j] += sigma_y;
    }
  }

  // write_Matrix_To_File(sigma, "../../data/sigma.txt");
  // exit(0);
}

void FDTD_2D::set_PML_kappa(){
  double kappa_max = 5.0; // or 10.0
  double grading_order = 3.0;
  // Left boundary
  for (int i = 0; i < pml_layers; i++) {
    kappa_x[i] = 1.0 + (kappa_max - 1.0) * pow((pml_layers - i - 1) / (double)pml_layers, grading_order);
  }
  
  // Right boundary
  for (int i = nx - pml_layers; i < nx; i++) {
    kappa_x[i] = 1.0 + (kappa_max - 1.0) * pow((i - (nx - pml_layers - 1)) / (double)pml_layers, grading_order);
  }
  
  // Bottom boundary
  for (int j = 0; j < pml_layers; j++) {
    kappa_y[j] = 1.0 + (kappa_max - 1.0) * pow((pml_layers - j - 1) / (double)pml_layers, grading_order);
  }
  
  // Top boundary
  for (int j = ny - pml_layers; j < ny; j++) {
    kappa_y[j] = 1.0 + (kappa_max - 1.0) * pow((j - (ny - pml_layers - 1)) / (double)pml_layers, grading_order);
  }
  
  // write_Matrix_To_File(sigma, "../../data/sigma.txt");
}

void FDTD_2D::set_PML_a(){
  double grading_order = 3.0;
  double a_max = 0.05;  

// Left boundary
for (int i = 0; i < pml_layers; i++) {
  double x = (pml_layers - i - 1) / (double)pml_layers;
  pml_a_x[i] = a_max * pow(x, grading_order);
}

// Right boundary
for (int i = nx - pml_layers; i < nx; i++) {
  double x = (i - (nx - pml_layers - 1)) / (double)pml_layers;
  pml_a_x[i] = a_max * pow(x, grading_order);
}

// Bottom boundary
for (int j = 0; j < pml_layers; j++) {
  double y = (pml_layers - j - 1) / (double)pml_layers;
  pml_a_y[j] = a_max * pow(y, grading_order);
}

// Top boundary
for (int j = ny - pml_layers; j < ny; j++) {
  double y = (j - (ny - pml_layers - 1)) / (double)pml_layers;
  pml_a_y[j] = a_max * pow(y, grading_order);
}

}

void FDTD_2D::set_PML_b() {
  // Compute pml_b_x
  for (int i = 0; i < nx; i++) {
      if (kappa_x[i] != 0.0) {
          double loss_term = (sigma[i][0] / kappa_x[i]) + pml_a_x[i];
          pml_b_x[i] = std::exp(-loss_term * dt / epsilon_0);
      } else {
          pml_b_x[i] = 1.0;
      }
  }

  // Compute pml_b_y
  for (int j = 0; j < ny; j++) {
      if (kappa_y[j] != 0.0) {
          double loss_term = (sigma[0][j] / kappa_y[j]) + pml_a_y[j];
          pml_b_y[j] = std::exp(-loss_term * dt * eta_0);
      } else {
          pml_b_y[j] = 1.0;
      }
  }
}

void FDTD_2D::set_PML_c() {
  for (int i = 0; i < nx; i++) {
      double sigma_x = sigma[i][0];
      double kappa = kappa_x[i];
      double a = pml_a_x[i];
      double b = pml_b_x[i];

      double denom = sigma_x * kappa + kappa * kappa * a;
      if (denom != 0.0)
          pml_c_x[i] = sigma_x * (b - 1.0) / denom;
      else
          pml_c_x[i] = 0.0;
  }

  for (int j = 0; j < ny; j++) {
      double sigma_y = sigma[0][j];
      double kappa = kappa_y[j];
      double a = pml_a_y[j];
      double b = pml_b_y[j];

      double denom = sigma_y * kappa + kappa * kappa * a;
      if (denom != 0.0)
          pml_c_y[j] = sigma_y * (b - 1.0) / denom;
      else
          pml_c_y[j] = 0.0;
  }
}

void FDTD_2D::set_PML_layers(){
  set_PML_sigma();
  set_PML_kappa();
  set_PML_a();
  set_PML_b();
  set_PML_c();
}

void FDTD_2D::write_Matrix_To_File(const vector<vector<double>>& field, const string& filename){
 ofstream fout(filename);
  if (!fout.is_open()) {
      cerr << "Error: Cannot open file " << filename << endl;
      return;
  }

  for (const auto& row : field) {
      for (const auto& val : row) {
          fout << val << " ";
      }
      fout << "\n";
  }

  fout.close();
  
}

void FDTD_2D::fdtd_2d_basic_TEz(){
   
  for(int t = 0;t<nt;t++) {
     
    for(int i = 1;i<nx;i++) { 
      for(int j = 1;j<ny;j++) {
        Dz[i][j] += (dt/delta_x)*(Hy[i][j] - Hy[i-1][j] - Hx[i][j] + Hx[i][j-1]);    
      }
    }
  
    for(int i = 1;i<nx;i++) { 
      for(int j = 1;j<ny;j++) {
        Ez[i][j] = epsilon_r[i][j] * Dz[i][j];    
      }
    }
  
    point_source(Ez, t);
  
    //hy and Hz
    for(int i = 0;i<nx-1;i++) {
      for(int j = 0;j<ny-1;j++) {
        Hx[i][j] += (dt/delta_y)*(Ez[i][j] - Ez[i][j+1]);
      }
    }

    
    for(int i = 0;i<nx-1;i++) {
      for(int j = 0;j<ny-1;j++) {
        Hy[i][j] += (dt/delta_x)*(Ez[i+1][j] - Ez[i][j]);
      }
    }
    if(t%data_capture_interval == 0){
      write_Matrix_To_File(Ez, "../../data/Ez/Ez_" + std::to_string(t) + ".txt");
    }
  }

}


void FDTD_2D::fdtd_2d_basic_TEz_PML() {
  for (int t = 0; t < nt; t++) {

    // --- Update Dz with Ez CPML ---
    for (int i = 1; i < nx; i++) {
      for (int j = 1; j < ny; j++) {
        double dHy_dx = (Hy[i][j] - Hy[i - 1][j]) / delta_x;
        double dHx_dy = (Hx[i][j] - Hx[i][j - 1]) / delta_y;

        psiEz_x[i][j] = pml_b_x[i] * psiEz_x[i][j] + pml_c_x[i] * dHy_dx;
        psiEz_y[i][j] = pml_b_y[j] * psiEz_y[i][j] + pml_c_y[j] * dHx_dy;

        Dz[i][j] += dt * ((1.0 / kappa_x[i]) * dHy_dx + psiEz_x[i][j] - (1.0 / kappa_y[j]) * dHx_dy - psiEz_y[i][j]);
      }
    }

    // --- Update Ez from Dz ---
    for (int i = 1; i < nx; i++) {
      for (int j = 1; j < ny; j++) {
        Ez[i][j] = Dz[i][j] / epsilon_r[i][j];
      }
    }

    point_source(Ez, t);

    // --- Update Hx with CPML ---
    for (int i = 0; i < nx - 1; i++) {
      for (int j = 0; j < ny - 1; j++) {
        double dEz_dy = (Ez[i][j] - Ez[i][j + 1]) / delta_y;

        psiHx_y[i][j] = pml_b_y[j] * psiHx_y[i][j] + pml_c_y[j] * dEz_dy;

        Hx[i][j] += dt * (
          (1.0 / (mu_0 * kappa_y[j])) * dEz_dy + psiHx_y[i][j]
        );
      }
    }

    // --- Update Hy with CPML ---
    for (int i = 0; i < nx - 1; i++) {
      for (int j = 0; j < ny - 1; j++) {
        double dEz_dx = (Ez[i + 1][j] - Ez[i][j]) / delta_x;
        psiHy_x[i][j] = pml_b_x[i] * psiHy_x[i][j] + pml_c_x[i] * dEz_dx;
        Hy[i][j] += dt * ((1.0 / (mu_0 * kappa_x[i])) * dEz_dx + psiHy_x[i][j]);
      }
    }

    if (t % data_capture_interval == 0) {
      write_Matrix_To_File(Ez, "../../data/Ez/Ez_" + std::to_string(t) + ".txt");
    }
  }
}


// void FDTD_2D::fdtd_2d_basic_TMz(){
   
//   for(int t = 0;t<nt;t++) {
     
//     for(int i = 1;i<nx;i++) { 
//       for(int j = 1;j<nx;j++) {
//         dz[i][j] += 0.5*(hy[i][j] - hy[i-1][j] - hx[i][j] + hx[i][j-1]);    
//       }
//     }
  
//     for(int i = 1;i<nx;i++) { 
//       for(int j = 1;j<nx;j++) {
//         ez[i][j] = epsilon_r[i][j] * dz[i][j];    
//       }
//     }
  
//     ez[nx/2][nx/2] = source_function(t);
  
//     //hy and Hz
//     for(int i = 0;i<nx-1;i++) {
//       for(int j = 0;j<nx-1;j++) {
//         hx[i][j] += 0.5*(ez[i][j] - ez[i][j+1]);
//       }
//     }

    
//     for(int i = 0;i<nx-1;i++) {
//       for(int j = 0;j<nx-1;j++) {
//         hy[i][j] += 0.5*(ez[i+1][j] - ez[i][j]);
//       }
//     }

//     write_Matrix_To_File(ez, "../data/Ez/Ez_" + std::to_string(t) + ".txt");

//   }

// }



// void FDTD_2D::fdtd_2d_PML(){
  
//   for(int t = 0;t<nt;t++) {
     
//     for(int i = 1;i<nx;i++) { 
//       for(int j = 1;j<nx;j++) {
//         dz[i][j] += 0.5*(hy[i][j] - hy[i-1][j] - hx[i][j] + hx[i][j-1]);    
//       }
//     }
  
//     for(int i = 1;i<nx;i++) { 
//       for(int j = 1;j<nx;j++) {
//         ez[i][j] = epsilon_r[i][j] * dz[i][j];    
//       }
//     }
  
//     ez[nx/2][nx/2] = exp(-0.5*(pow((t0 - t)/spread, 2)));
  
//     //hy and Hz
//     for(int i = 0;i<nx-1;i++) {
//       for(int j = 0;j<nx-1;j++) {
//         hx[i][j] += 0.5*(ez[i][j] - ez[i][j+1]);
//       }
//     }

    
//     for(int i = 0;i<nx-1;i++) {
//       for(int j = 0;j<nx-1;j++) {
//         hy[i][j] += 0.5*(ez[i+1][j] - ez[i][j]);
//       }
//     }

//     write_Matrix_To_File(ez, "../data/Ez/Ez_" + std::to_string(t) + ".txt");

//   }


// }



