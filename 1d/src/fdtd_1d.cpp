#include "fdtd_1d.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

#include <sstream>
#include <json/json.h>

FDTD_1D::FDTD_1D(const std::string &control_file) {
	load_parameters(control_file);
	setup_grid_and_arrays();
	setup_pml();
	setup_coefficients();
}

void FDTD_1D::load_parameters(const std::string &control_file) {
	std::ifstream ifs(control_file, std::ifstream::binary);
	if (!ifs.is_open()) {
			throw std::runtime_error("Cannot open control file: " + control_file);
	}

	Json::Value root;
	Json::CharReaderBuilder builder;
	std::string errs;
	if (!Json::parseFromStream(builder, ifs, &root, &errs)) {
			throw std::runtime_error("JSON parse error: " + errs);
	}

    
	const auto &s = root["simulation"];
	simulation_size = s["simulation_size"].asDouble();
	frequency = s["frequency"].asDouble();
	cells_per_wavelength = s["cells_per_wavelength"].asInt();
	CFL = s["CFL"].asDouble();
	simulation_time = s["simulation_time"].asDouble();
	data_capture_interval = s["data_capture_interval"].asInt();

	const auto &p = root["pml"];
	pml_size = p["pml_size"].asInt();
	m = p["m"].asDouble();
	a_max = p["a_max"].asDouble();
	kappa_max = p["kappa_max"].asDouble();

	const auto &so = root["source"];
	source_position = so["source_position"].asInt();
	pulse_width = so["pulse_width"].asDouble();
	pulse_delay = 4 * pulse_width;  // derived
	amplitude = so["amplitude"].asDouble();
    
}


void FDTD_1D::setup_grid_and_arrays() {
	// --- Step size (from wavelength & cells per wavelength) ---
	double center_wavelength = c0 / frequency;
	step_size = center_wavelength / cells_per_wavelength;

	// --- Number of spatial cells ---
	N_space_cells = static_cast<int>(simulation_size / step_size);

	// --- Time step & number of steps ---
	dt           = CFL * step_size / c0;
	N_time_steps = static_cast<int>(simulation_time / dt);

	// --- Allocate arrays ---
	Ex.assign(N_space_cells, 0.0);
	Ex_prev.assign(N_space_cells, 0.0);
	Hz.assign(N_space_cells, 0.0);
	Hz_prev.assign(N_space_cells, 0.0);

	P_Ex_y.assign(N_space_cells, 0.0);
	P_Hz_y.assign(N_space_cells, 0.0);
	eps.assign(N_space_cells, 1.0);
	sigma.assign(N_space_cells, 0.0);
	sigma_h.assign(N_space_cells, 0.0);
	refractive_index.assign(N_space_cells, 1.0);
	inv_kappa_dy.assign(N_space_cells, 1.0 / step_size);
	inv_kappa_h_dy.assign(N_space_cells, 1.0 / step_size);

	// --- PML arrays ---
	sigma_pml.assign(pml_size, 0.0);
	sigma_h_pml.assign(pml_size, 0.0);
	kappa_far.assign(pml_size, 0.0);
	kappa_h_far.assign(pml_size, 0.0);
	a_pml.assign(pml_size, 0.0);
	a_h_pml.assign(pml_size, 0.0);
	kappa_left.assign(pml_size, 0.0);
	kappa_h_left.assign(pml_size, 0.0);

	be_y_f.assign(pml_size, 0.0);
	ce_y_f.assign(pml_size, 0.0);
	bh_y_f.assign(pml_size, 0.0);
	ch_y_f.assign(pml_size, 0.0);

	be_y.assign(pml_size, 0.0);
	ce_y.assign(pml_size, 0.0);
	bh_y.assign(pml_size, 0.0);
	ch_y.assign(pml_size, 0.0);

	// --- Field update coeffs ---
	e_coeff1.assign(N_space_cells, 1.0);
	e_coeff2.assign(N_space_cells, 1.0);
	h_coeff1.assign(N_space_cells, 1.0);
	h_coeff2.assign(N_space_cells, 1.0);

	// --- Source setup ---
	omega0 = 2.0 * M_PI * frequency;

	// ---- Debug print ----
	
	std::cout << "simulation_size       = " << simulation_size << "\n";
	std::cout << "frequency             = " << frequency << "\n";
	std::cout << "cells_per_wavelength  = " << cells_per_wavelength << "\n";
	std::cout << "CFL                   = " << CFL << "\n";
	std::cout << "simulation_time       = " << simulation_time << "\n";
	std::cout << "data_capture_interval = " << data_capture_interval << "\n\n";

	std::cout << "pml_size   = " << pml_size << "\n";
	std::cout << "m          = " << m << "\n";
	std::cout << "a_max      = " << a_max << "\n";
	std::cout << "kappa_max  = " << kappa_max << "\n\n";

	std::cout << "source_position = " << source_position << "\n";
	std::cout << "pulse_width     = " << pulse_width << "\n";
	std::cout << "pulse_delay     = " << pulse_delay << "\n";
	std::cout << "amplitude       = " << amplitude << "\n";
	std::cout << "------------------------------------\n";


	std::cout << "[Grid setup]\n";
	std::cout << "N_space_cells = " << N_space_cells << "\n";
	std::cout << "N_time_steps  = " << N_time_steps << "\n";
	std::cout << "step_size     = " << step_size << "\n";
	std::cout << "dt            = " << dt << "\n";
	
}


void FDTD_1D::setup_pml() {
// NOTE: match Python: log(1e-10), not 1e-9
	double n_edge     = refractive_index[N_space_cells - 1 - pml_size];
	double pml_cond_e = -(m + 1.0) * log(1e-10) * c0 * EPSILON_0 / (2.0 * n_edge * pml_size * step_size);

// ---------- Right (+y) PML ----------
	for (int i = 0; i < pml_size; i++) {
		double x = static_cast<double>(i) / pml_size; // 0 → 1
		sigma_pml[i]   = pml_cond_e * pow(x, m);
		kappa_far[i]   = 1.0 + (kappa_max - 1.0) * pow(x, m);
		a_pml[i]       = a_max * ((pml_size - i) / static_cast<double>(pml_size));

		sigma_h_pml[i] = pml_cond_e * pow((i + 0.5) / pml_size, m);
		kappa_h_far[i] = 1.0 + (kappa_max - 1.0) * pow((i + 0.5) / pml_size, m);
		a_h_pml[i]     = a_max * ((pml_size - i + 0.5) / static_cast<double>(pml_size));

		be_y_f[i] = exp(-(sigma_pml[i] / kappa_far[i] + a_pml[i]) * dt / EPSILON_0);
		ce_y_f[i] = sigma_pml[i] * (be_y_f[i] - 1.0) /
								((sigma_pml[i] + kappa_far[i] * a_pml[i]) * kappa_far[i]);

		bh_y_f[i] = exp(-(sigma_h_pml[i] / kappa_h_far[i] + a_h_pml[i]) * dt / EPSILON_0);
		ch_y_f[i] = sigma_h_pml[i] * (bh_y_f[i] - 1.0) /
								((sigma_h_pml[i] + kappa_h_far[i] * a_h_pml[i]) * kappa_h_far[i]);

		// Right edge inverse-kappa
		inv_kappa_dy[N_space_cells - pml_size + i]   = 1.0 / (kappa_far[i] * step_size);
		inv_kappa_h_dy[N_space_cells - pml_size + i] = 1.0 / (kappa_h_far[i] * step_size);
		
	}	

// ---------- Left (−y) PML ----------
	for (int i = 0; i < pml_size; i++) {
		double x = static_cast<double>(pml_size - i) / pml_size; // 1 → 1/pml_size
		sigma_pml[i]   = pml_cond_e * pow(x, m);
		kappa_left[i]  = 1.0 + (kappa_max - 1.0) * pow(x, m);
		a_pml[i]       = a_max * (static_cast<double>(i) / pml_size);

		sigma_h_pml[i] = pml_cond_e * pow((pml_size - i - 0.5) / pml_size, m);
		kappa_h_left[i]= 1.0 + (kappa_max - 1.0) * pow((pml_size - i - 0.5) / pml_size, m);
		a_h_pml[i]     = a_max * ((i + 0.5) / pml_size);

		be_y[i] = exp(-(sigma_pml[i] / kappa_left[i] + a_pml[i]) * dt / EPSILON_0);
		ce_y[i] = sigma_pml[i] * (be_y[i] - 1.0) /
							((sigma_pml[i] + kappa_left[i] * a_pml[i]) * kappa_left[i]);

		bh_y[i] = exp(-(sigma_h_pml[i] / kappa_h_left[i] + a_h_pml[i]) * dt / EPSILON_0);
		ch_y[i] = sigma_h_pml[i] * (bh_y[i] - 1.0) /
							((sigma_h_pml[i] + kappa_h_left[i] * a_h_pml[i]) * kappa_h_left[i]);

		inv_kappa_dy[i]   = 1.0 / (kappa_left[i] * step_size);
		inv_kappa_h_dy[i] = 1.0 / (kappa_h_left[i] * step_size);
	
	}

	
}

void FDTD_1D::setup_coefficients() {
	for (int i = 0; i < N_space_cells; i++) {
		// For E-field update
		double denom_e = EPSILON_0 * eps[i] / dt  + sigma[i]/2.0;
		e_coeff1[i] = (EPSILON_0 * eps[i] / dt - sigma[i]/2.0) / denom_e;
		e_coeff2[i] = (1.0 / denom_e);

		// For H-field update
		double denom_h = MU_0 / dt + sigma_h[i]/2;
		h_coeff1[i] = (MU_0 / dt - sigma_h[i]/2) / denom_h;
		h_coeff2[i] = (1.0 / denom_h);
	}

	std::cout << "[Coefficients setup done]\n";
}

double FDTD_1D::signal(double time) {
	// Gaussian envelope
	double envelope = exp(-pow((time - pulse_delay) / pulse_width, 2));

	// Sinusoidal carrier
	return envelope * amplitude * sin(omega0 * time);
}

void FDTD_1D::output_fields(int timestep) {
// Write Ex to file
	std::ofstream exf("data/Ex/ElectricField_" + std::to_string(timestep) + ".txt");
	for (int i = 0; i < N_space_cells; i++) {
		exf << Ex[i] << '\n';
	}

// Write Hz to file (optional)
	std::ofstream hzf("data/Hz/MagneticField_" + std::to_string(timestep) + ".txt");
	for (int i = 0; i < N_space_cells; i++) {
		hzf << Hz[i] << '\n';
	}	

// Diagnostics
	double maxE = *std::max_element(Ex.begin(), Ex.end());
	double minE = *std::min_element(Ex.begin(), Ex.end());
	std::cout << "Timestep: " << timestep
						<< "  Ex_min = " << minE
						<< "  Ex_max = " << maxE << "\n";
}

void FDTD_1D::run() {
	source_position = int(N_space_cells/2);
	for (int n = 0; n < N_time_steps; n++) {
		// --- Save previous fields ---
		Hz_prev = Hz;
		Ex_prev = Ex;

		// --- Update Magnetic CPML fields ---
		for (int i = N_space_cells - pml_size; i < N_space_cells - 1; i++) {
				int idx = i - (N_space_cells - pml_size);
				P_Hz_y[i] = bh_y_f[idx] * P_Hz_y[i]
									+ ch_y_f[idx] * (Ex[i + 1] - Ex[i]) / step_size;
		}
		for (int i = 0; i < pml_size; i++) {
				P_Hz_y[i] = bh_y[i] * P_Hz_y[i]
									+ ch_y[i] * (Ex[i + 1] - Ex[i]) / step_size;
		}

		// --- Update Magnetic field Hz ---
		for (int j = 0; j < N_space_cells - 1; j++) {
				Hz[j] = h_coeff1[j] * Hz_prev[j]
							+ h_coeff2[j] * (inv_kappa_h_dy[j] * (Ex[j + 1] - Ex[j]) + P_Hz_y[j]);
		}

		// --- Update Electric CPML fields ---
		for (int i = N_space_cells - pml_size; i < N_space_cells; i++) {
				int idx = i - (N_space_cells - pml_size);
				P_Ex_y[i] = be_y_f[idx] * P_Ex_y[i]
									+ ce_y_f[idx] * (Hz[i] - Hz[i - 1]) / step_size;
		}
		for (int i = 1; i < pml_size; i++) {
				P_Ex_y[i] = be_y[i] * P_Ex_y[i]
									+ ce_y[i] * (Hz[i] - Hz[i - 1]) / step_size;
		}

		// --- Update Electric field Ex ---
		for (int j = 1; j < N_space_cells - 1; j++) {
				Ex[j] = e_coeff1[j] * Ex_prev[j]
							+ e_coeff2[j] * (inv_kappa_dy[j] * (Hz[j] - Hz[j - 1]) + P_Ex_y[j]);
		}

		// --- Inject Source ---
		Ex[source_position] += signal((n + 1) * dt);

		// --- Output every 100 steps ---
		if (n % data_capture_interval == 0) {
				output_fields(n);
		}
	}
}
