#include <iostream> 
#include <vector> 
#include <cmath>
#include <string>
#include <fstream>
#include <jsoncpp/json/value.h>
#include <json/json.h>

#include "FDTD_1D.h"

using namespace std;

// Constructor
FDTD_1D::FDTD_1D(string filename) {
  ifstream file(filename);
  Json::Value jroot;
  Json::CharReaderBuilder builder;
  string errs;

  Json::parseFromStream(builder, file, &jroot, &errs);

  frequency = jroot["frequency"].asDouble();
  CFL = jroot["CFL"].asDouble();

  domain.resize(2);
  domain[0] = jroot["domain_x"][0].asDouble();
  domain[1] = jroot["domain_x"][1].asDouble();  
  
  bandwidth_percent = jroot["bandwidth"].asInt();

  source_position = jroot["source_position"][0].asDouble();  

  if (source_position > domain[1] || source_position < domain[0]) {
    cerr << "Source position is out of domain" << endl;
    exit(1);
  }

  nt = jroot["iterations"].asInt();  
  data_capture_interval = jroot["data_capture_interval"].asInt();

  source_type = jroot["source"].asString();

}

// Destructor
FDTD_1D::~FDTD_1D() {}

void FDTD_1D::init() {
  set_spacing(); 
  set_dt();
  set_fields();
  setup_source_function();
  set_PML_layers();
  // set_sigma_graded();
  set_field_coeffs();
  // exit(0);
  print_initial_state();
  
}

void FDTD_1D::print_initial_state() {
  cout<<"CFL = "<<CFL<<endl;
  cout<<"Frequency  = "<<frequency<<endl;
  cout<<"Wavelength = "<<c0/frequency<<endl;
  cout<<"Grid size = "<<delta_z<<endl;
  cout<<"Grid points = "<<nz<<endl;
  cout<<"Number of PML layers = "<<pml_layers<<endl;  
  cout<<"domain = "<<domain[0]<<" to "<<domain[1]<<endl;
  cout<<"source_type = "<<source_type<<endl;
  cout<<"Number of iterations = "<<nt<<endl;
  cout<<"Data capture interval = "<<data_capture_interval<<endl;
  
}

void FDTD_1D::set_spacing() {
  delta_z = c0 / (10 * frequency);  
  double wavelength = c0 / frequency;
  pml_layers =  int(3*wavelength / delta_z);  
  // pml_layers = 25;  
  nz = int((domain[1] - domain[0]) / delta_z) + 2 * pml_layers;
  cout << "nz: " << nz << endl;
}

void FDTD_1D::set_dt() {
  dt = CFL * delta_z / c0;  
  bandwidth = frequency * bandwidth_percent / 100;
  spread = 1 / (bandwidth * 2 * M_PI);
  tau = 3 * spread;
}

void FDTD_1D::set_fields() {
  Ex.resize(nz, 0.0);
  Hy.resize(nz, 0.0); 
  epsilon_r.resize(nz, 1);
  mu_r.resize(nz, 1);  
  sigma_z_e.resize(nz, 0.0);
  kappa_z_e.resize(nz, 1.0);
  pml_a_z_e.resize(nz, 0.0);
  pml_b_z_e.resize(nz, 1.0);
  pml_c_z_e.resize(nz, 0.0);
  sigma_z_h.resize(nz, 0.0);
  kappa_z_h.resize(nz, 1.0);
  pml_a_z_h.resize(nz, 0.0);
  pml_b_z_h.resize(nz, 1.0);
  pml_c_z_h.resize(nz, 0.0);
  psiEx_z.resize(nz, 0);
  psiHy_z.resize(nz, 0);
  e1_coeff.resize(nz, 0);
  e2_coeff.resize(nz, 0);
  h1_coeff.resize(nz, 0);
  h2_coeff.resize(nz, 0);
}


void FDTD_1D::setup_source_function() {
  source_x = pml_layers + int((source_position - domain[0]) / delta_z);

  if (source_type == "gaussian") {
    source_function = [this](int t) {
      return exp(-0.5 * pow((t * dt - tau) / spread, 2));
    };
  } else if (source_type == "sin") {
    source_function = [this](int t) {
      return sin(2 * M_PI * frequency * t * dt);
    };
  } else {
    cerr << "Unknown source type: " << source_type << endl;
    exit(1);
  }
}

void FDTD_1D::set_sigma_graded() {
  
  double grading_order_sigma = 3;
  double Delta = 1e-12;
  
  double sigma_max = -(grading_order_sigma + 1) * log(Delta) * c0 * epsilon_0/(2 * pml_layers * delta_z);
  // Right PML
  for (int i = nz - pml_layers; i < nz; i++) {
    double y = (i - (nz - pml_layers)) / (double)pml_layers;
    sigma_z_e[i] = sigma_max * pow(y, grading_order_sigma);    
  }

  // Left PML
  for (int i = 0; i < pml_layers; i++) {
    double y = (pml_layers - i) / (double)pml_layers;
    sigma_z_e[i] = sigma_max * pow(y, grading_order_sigma);
  }

  for (int i = 0; i < nz - 1; i++) {
    sigma_z_h[i] = (mu_0 / epsilon_0) * 0.5 * (sigma_z_e[i] + sigma_z_e[i + 1]);
  }


  write_Vector_To_File(sigma_z_e, "../data/sigma_z_e_graded.txt");
  write_Vector_To_File(sigma_z_h, "../data/sigma_z_h_graded.txt");

}

void FDTD_1D::set_PML_layers() {
  double grading_order_sigma = 3.0;
  double kappa_max = 5.0, grading_order_kappa = 3.0;
  double a_max = 0.10, grading_order_a = 1.0;

  double Delta = 1e-12;
  double sigma_max = -(grading_order_sigma + 1) * log(Delta) * c0 * epsilon_0 / (2 * pml_layers * delta_z);

  // Left PML (Electric field)
  for (int i = 1; i < pml_layers; i++) {
    double z = (pml_layers - i) / static_cast<double>(pml_layers);
    double z_for_a = 1.0 - (z / static_cast<double>(pml_layers));
    sigma_z_e[i] = sigma_max * pow(z, grading_order_sigma);
    kappa_z_e[i] = 1.0 + (kappa_max - 1.0) * pow(z, grading_order_kappa);
    pml_a_z_e[i] = a_max * pow(z_for_a, grading_order_a);

    double loss = (sigma_z_e[i] / kappa_z_e[i]) + pml_a_z_e[i];
    pml_b_z_e[i] = exp(-loss * dt / epsilon_0);
    double denom = sigma_z_e[i] * kappa_z_e[i] + pow(kappa_z_e[i], 2) * pml_a_z_e[i];
    pml_c_z_e[i] = (denom != 0.0) ? sigma_z_e[i] * (pml_b_z_e[i] - 1.0) / denom : 0.0;
  }

  
  // Right PML (Electric field)
  for (int i = nz - pml_layers + 1; i < nz; i++) {
    double z = (i - (nz - pml_layers)) / static_cast<double>(pml_layers);
    double z_for_a = 1.0 - (z / static_cast<double>(pml_layers));
    sigma_z_e[i] = sigma_max * pow(z, grading_order_sigma);
    kappa_z_e[i] = 1.0 + (kappa_max - 1.0) * pow(z, grading_order_kappa);
    pml_a_z_e[i] = a_max * pow(z_for_a, grading_order_a);

    double loss = (sigma_z_e[i] / kappa_z_e[i]) + pml_a_z_e[i];
    pml_b_z_e[i] = exp(-loss * dt / epsilon_0);
    double denom = sigma_z_e[i] * kappa_z_e[i] + pow(kappa_z_e[i], 2) * pml_a_z_e[i];
    pml_c_z_e[i] = (denom != 0.0) ? sigma_z_e[i] * (pml_b_z_e[i] - 1.0) / denom : 0.0;
  }

  for (int i = 0; i < nz - 1; i++) {
    sigma_z_h[i] = (mu_0 / epsilon_0) * 0.5 * (sigma_z_e[i] + sigma_z_e[i + 1]);
    kappa_z_h[i] = 0.5*(kappa_z_e[i] + kappa_z_e[i+1]);
    pml_a_z_h[i] = 0.5*(pml_a_z_e[i] + pml_a_z_e[i+1]);
    pml_b_z_h[i] = 0.5*(pml_b_z_e[i] + pml_b_z_e[i+1]);
    pml_c_z_h[i] = 0.5*(pml_c_z_e[i] + pml_c_z_e[i+1]);
  }

  // // Left PML (Magnetic field)
  // for (int i = 0; i < pml_layers-1; i++) {
  //   double z = (pml_layers - i - 0.5) / static_cast<double>(pml_layers);
  //   double z_for_a = 1.0 - (z / static_cast<double>(pml_layers));
  //   sigma_z_h[i] = (mu_0 / epsilon_0) * sigma_max * pow(z, grading_order_sigma);
  //   kappa_z_h[i] = 1.0 + (kappa_max - 1.0) * pow(z, grading_order_kappa);
  //   pml_a_z_h[i] = a_max * pow(z_for_a, grading_order_a);

  //   double loss = (sigma_z_h[i] / kappa_z_h[i]) + pml_a_z_h[i];
  //   pml_b_z_h[i] = exp(-loss * dt / epsilon_0);
  //   double denom = sigma_z_h[i] * kappa_z_h[i] + pow(kappa_z_h[i], 2) * pml_a_z_h[i];
  //   pml_c_z_h[i] = (denom != 0.0) ? sigma_z_h[i] * (pml_b_z_h[i] - 1.0) / denom : 0.0;
  // }

  // // Right PML (Magnetic field)
  // for (int i = nz - pml_layers; i < nz-1; i++) {
  //   double z = (i - (nz - pml_layers) - 0.5) / static_cast<double>(pml_layers);
  //   double z_for_a = 1.0 - (z / static_cast<double>(pml_layers));
  //   sigma_z_h[i] = (mu_0 / epsilon_0) * sigma_max * pow(z, grading_order_sigma);
  //   kappa_z_h[i] = 1.0 + (kappa_max - 1.0) * pow(z, grading_order_kappa);
  //   pml_a_z_h[i] = a_max * pow(z_for_a, grading_order_a);

  //   double loss = (sigma_z_h[i] / kappa_z_h[i]) + pml_a_z_h[i];
  //   pml_b_z_h[i] = exp(-loss * dt / epsilon_0);
  //   double denom = sigma_z_h[i] * kappa_z_h[i] + pow(kappa_z_h[i], 2) * pml_a_z_h[i];
  //   pml_c_z_h[i] = (denom != 0.0) ? sigma_z_h[i] * (pml_b_z_h[i] - 1.0) / denom : 0.0;
  // }

  write_Vector_To_File(sigma_z_e, "../data/sigma_z_e_pml.txt");
  // write_Vector_To_File(kappa_z_e, "../data/kappa_z_e.txt");
  // write_Vector_To_File(pml_a_z_e, "../data/pml_a_z_e.txt");
  // write_Vector_To_File(pml_b_z_e, "../data/pml_b_z_e.txt");
  // write_Vector_To_File(pml_c_z_e, "../data/pml_c_z_e.txt");

  write_Vector_To_File(sigma_z_h, "../data/sigma_z_h_pml.txt");
  // write_Vector_To_File(kappa_z_h, "../data/kappa_z_h.txt");
  // write_Vector_To_File(pml_a_z_h, "../data/pml_a_z_h.txt");
  // write_Vector_To_File(pml_b_z_h, "../data/pml_b_z_h.txt");
  // write_Vector_To_File(pml_c_z_h, "../data/pml_c_z_h.txt");
}

void FDTD_1D::set_field_coeffs() {
  
  for(int i = 0;i<nz;i++) {

    double numerator     = epsilon_0 / dt - sigma_z_e[i] / 2.0;
    double denominator   = epsilon_0 / dt + sigma_z_e[i] / 2.0;
    double numerator_h   = mu_0 / dt - sigma_z_h[i] / 2.0;
    double denominator_h = mu_0 / dt + sigma_z_h[i] / 2.0;
    e1_coeff[i] = numerator/denominator;
    e2_coeff[i] = 1/(denominator);
    h1_coeff[i] = numerator_h/denominator_h;
    h2_coeff[i] = 1/(denominator_h);
  }
  
}

void FDTD_1D::writeArrayToFile(const string field, const string data_directory, int t) {
  string filename = data_directory + field + "_" + to_string(t) + ".txt";
  ofstream outFile(filename);
  if (!outFile) {
    cerr << "Error opening file: " << filename << "\n";
    return;
  }

  if (field == "ElectricField") {
    for (double val : Ex) outFile << val << "\n";
  } else if (field == "MagneticField") {
    for (double val : Hy) outFile << val << "\n";
  }

  outFile.close();
}

void FDTD_1D::write_Vector_To_File(vector<double>& field, string filename) {
  ofstream fout(filename);
  if (!fout.is_open()) {
      cerr << "Error: Cannot open file " << filename << endl;
      return;
  }

  for (const auto& val : field) {
      fout << val << "\n";
  }

  fout.close();
}

void FDTD_1D::fdtd_1d_basic() {

  for (int t = 0; t < nt; t++) {    
    for (int i = 1; i < nz; i++) {
      Ex[i] += (dt/(delta_z*epsilon_0)) * (Hy[i - 1] - Hy[i]);
    }

    Ex[source_x] += source_function(t);

    for (int i = 0; i < nz - 1; i++) {
      Hy[i] += (dt/(delta_z*mu_0)) * (Ex[i] - Ex[i + 1]);
    }

    if(t%data_capture_interval == 0){
      writeArrayToFile("ElectricField", "../data/Ex/", t);
      // writeArrayToFile("MagneticField", "../data/Hy/", t);
    }
  }
}

void FDTD_1D::fdtd_1d_graded_pml() {
  
  // set_sigma_graded();
  // set_field_coeffs();

  for (int t = 0; t < nt; t++) {

    for (int i = 1; i < nz - 1; i++) {
      Ex[i] = e1_coeff[i] * Ex[i]  + (e2_coeff[i]/delta_z) * (Hy[i - 1] - Hy[i]);
    }

    Ex[source_x] = source_function(t);

    for (int i = 0; i < nz - 1; i++) {
      Hy[i] = h1_coeff[i] * Hy[i] +  (h2_coeff[i]/delta_z) * (Ex[i] - Ex[i + 1]);
    }

    writeArrayToFile("ElectricField", "../data/Ex/", t);
    // writeArrayToFile("MagneticField", "../data/Hy/", t);
  }
}

void FDTD_1D::fdtd_1d_graded_Cpml() {

  // set_PML_layers();
  // set_field_coeffs();

  // for(int i = 0;i<nz;i++) {
  //   pm
  // }

  

  for (int t = 0; t < nt; t++) {

    // Update Electric Field
    for (int i = 1; i < nz; i++) {
      psiEx_z[i] = pml_b_z_e[i] * psiEx_z[i] + pml_c_z_e[i] * (Hy[i - 1] - Hy[i]) / delta_z;
      // psiEx_z[i] = 0;
      // kappa_z_e[i] = 1;
      Ex[i] = e1_coeff[i] * Ex[i] + e2_coeff[i] * ((Hy[i - 1] - Hy[i]) / (kappa_z_e[i] * delta_z) + psiEx_z[i]);
      // Ex[i] = e1_coeff[i] * Ex[i]  + (e2_coeff[i]/delta_z) * (Hy[i - 1] - Hy[i]);
    }

    // Add source
    Ex[source_x] = source_function(t);

    // Update Magnetic Field
    for (int i = 0; i < nz - 1; i++) {
      psiHy_z[i] = pml_b_z_h[i] * psiHy_z[i] + pml_c_z_h[i] * (Ex[i] - Ex[i + 1]) / delta_z;
      // psiHy_z[i] = 0;
      // kappa_z_h[i] = 1;
      Hy[i] = h1_coeff[i] * Hy[i] + h2_coeff[i] * ((Ex[i] - Ex[i + 1]) / (kappa_z_h[i] * delta_z) + psiHy_z[i]);
      // Hy[i] = h1_coeff[i] * Hy[i] +  (h2_coeff[i]/delta_z) * (Ex[i] - Ex[i + 1]);
    }

    // Save output
    if (t % data_capture_interval == 0) {
      writeArrayToFile("ElectricField", "../data/Ex/", t);
      // writeArrayToFile("MagneticField", "../data/Hy/", t);
    }
  }
}
