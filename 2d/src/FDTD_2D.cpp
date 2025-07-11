#include <iostream> 
#include <vector> 
#include <cmath>
#include <string>
#include <fstream>
#include "FDTD_2D.h"

#include <jsoncpp/json/value.h>
#include <json/json.h>
#include <functional>


using namespace std;

FDTD_2D::FDTD_2D(string filename) {
  ifstream file(filename);
  Json::Value root;
  Json::CharReaderBuilder builder;
  string errs;

  Json::parseFromStream(builder, file, &root, &errs);

  frequency = root["frequency"].asDouble();
  domain[0] = root["domain_x"][0].asDouble();
  domain[1] = root["domain_x"][1].asDouble();
  domain[2] = root["domain_y"][0].asDouble();
  domain[3] = root["domain_y"][1].asDouble();
  source_type = root["source"].asString();
  bandwidth = root["bandwidth"].asInt();
  tau = root["tau"].asInt();
  source_position[0] = root["source_position"][0].asDouble();
  source_position[1] = root["source_position"][1].asDouble();
  globalIteration = 0;
  nt = root["iterations"].asInt();  
  setup_source_function(); 
}

FDTD_2D::~FDTD_2D() = default;

void FDTD_2D::set_time_steps(int nt){
  this->nt = nt;
}

void FDTD_2D::set_spacing(int nx, int ny){
  this->nx = nx;
  this->ny = ny;
  this->nc = nx / 2;

  ex.resize(nx, vector<double>(ny, 0.0));
  ey.resize(nx, vector<double>(ny, 0.0));
  ez.resize(nx, vector<double>(ny, 0.0));
  dx.resize(nx, vector<double>(ny, 0.0));
  dy.resize(nx, vector<double>(ny, 0.0));
  dz.resize(nx, vector<double>(ny, 0.0));
  hx.resize(nx, vector<double>(ny, 0.0));
  hy.resize(nx, vector<double>(ny, 0.0));
  hz.resize(nx, vector<double>(ny, 0.0));
  
  epsilon_r.resize(nx, vector<double>(ny, 1.0));
  mu_r.resize(nx, vector<double>(ny, 1.0));
  sigma.resize(nx, vector<double>(ny, 0.0));
}

void FDTD_2D::set_excitation(double t0, double spread){
  this->t0 = t0;
  this->spread = spread;
}

void FDTD_2D::set_PML_layers(int pml_nx, int pml_ny){
  


  for(int i = nx-pml_nx;i<nx;i++){
    for(int j = ny-pml_ny;j<ny;j++){
      
    }
  }
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

void FDTD_2D::setup_source_function() {
    if (source_type == "gaussian") {
        source_function = [this](int t) {
            return exp(-0.5 * pow((t0 - t) / spread, 2));
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

void FDTD_2D::fdtd_2d_basic_TEz(){
   
  for(int t = 0;t<nt;t++) {
     
    for(int i = 1;i<nx;i++) { 
      for(int j = 1;j<nx;j++) {
        dz[i][j] += 0.5*(hy[i][j] - hy[i-1][j] - hx[i][j] + hx[i][j-1]);    
      }
    }
  
    for(int i = 1;i<nx;i++) { 
      for(int j = 1;j<nx;j++) {
        ez[i][j] = epsilon_r[i][j] * dz[i][j];    
      }
    }
  
    ez[nx/2][nx/2] = source_function(t);
  
    //hy and Hz
    for(int i = 0;i<nx-1;i++) {
      for(int j = 0;j<nx-1;j++) {
        hx[i][j] += 0.5*(ez[i][j] - ez[i][j+1]);
      }
    }

    
    for(int i = 0;i<nx-1;i++) {
      for(int j = 0;j<nx-1;j++) {
        hy[i][j] += 0.5*(ez[i+1][j] - ez[i][j]);
      }
    }

    write_Matrix_To_File(ez, "../../data/Ez/Ez_" + std::to_string(t) + ".txt");

  }

}




void FDTD_2D::fdtd_2d_basic_TMz(){
   
  for(int t = 0;t<nt;t++) {
     
    for(int i = 1;i<nx;i++) { 
      for(int j = 1;j<nx;j++) {
        dz[i][j] += 0.5*(hy[i][j] - hy[i-1][j] - hx[i][j] + hx[i][j-1]);    
      }
    }
  
    for(int i = 1;i<nx;i++) { 
      for(int j = 1;j<nx;j++) {
        ez[i][j] = epsilon_r[i][j] * dz[i][j];    
      }
    }
  
    ez[nx/2][nx/2] = exp(-0.5*(pow((t0 - t)/spread, 2)));
  
    //hy and Hz
    for(int i = 0;i<nx-1;i++) {
      for(int j = 0;j<nx-1;j++) {
        hx[i][j] += 0.5*(ez[i][j] - ez[i][j+1]);
      }
    }

    
    for(int i = 0;i<nx-1;i++) {
      for(int j = 0;j<nx-1;j++) {
        hy[i][j] += 0.5*(ez[i+1][j] - ez[i][j]);
      }
    }

    write_Matrix_To_File(ez, "../data/Ez/Ez_" + std::to_string(t) + ".txt");

  }

}



void FDTD_2D::fdtd_2d_PML(){
  
  for(int t = 0;t<nt;t++) {
     
    for(int i = 1;i<nx;i++) { 
      for(int j = 1;j<nx;j++) {
        dz[i][j] += 0.5*(hy[i][j] - hy[i-1][j] - hx[i][j] + hx[i][j-1]);    
      }
    }
  
    for(int i = 1;i<nx;i++) { 
      for(int j = 1;j<nx;j++) {
        ez[i][j] = epsilon_r[i][j] * dz[i][j];    
      }
    }
  
    ez[nx/2][nx/2] = exp(-0.5*(pow((t0 - t)/spread, 2)));
  
    //hy and Hz
    for(int i = 0;i<nx-1;i++) {
      for(int j = 0;j<nx-1;j++) {
        hx[i][j] += 0.5*(ez[i][j] - ez[i][j+1]);
      }
    }

    
    for(int i = 0;i<nx-1;i++) {
      for(int j = 0;j<nx-1;j++) {
        hy[i][j] += 0.5*(ez[i+1][j] - ez[i][j]);
      }
    }

    write_Matrix_To_File(ez, "../data/Ez/Ez_" + std::to_string(t) + ".txt");

  }


}



