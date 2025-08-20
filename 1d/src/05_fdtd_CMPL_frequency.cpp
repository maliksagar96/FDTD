#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace std;

double EPSILON_0 = 8.8541878128e-12;
double MU_0      = 1.256637062e-6;
double c0        = 2.99792458e8;
double imp0      = sqrt(MU_0 / EPSILON_0);

double simulation_size = 20e-6;
double frequency = 193414489.032258e6;
double center_wavelength = c0/frequency;
double step_size = center_wavelength/100;

int    N_space_cells   = int(simulation_size / step_size);
double CFL = 1;
double dt              = CFL * step_size / c0;
double simulation_time = 1e-12;
int    N_time_steps    = int(simulation_time / dt);

vector<double> Ex(N_space_cells, 0.0), Ex_prev(N_space_cells, 0.0);
vector<double> Hz(N_space_cells, 0.0), Hz_prev(N_space_cells, 0.0);
vector<double> P_Ex_y(N_space_cells, 0.0);
vector<double> P_Hz_y(N_space_cells, 0.0);
vector<double> eps(N_space_cells, 1.0);
vector<double> sigma(N_space_cells, 0.0);
vector<double> sigma_h(N_space_cells, 0.0);
vector<double> refractive_index(N_space_cells, 1.0);
vector<double> inv_kappa_dy(N_space_cells, 1.0);
vector<double> inv_kappa_h_dy(N_space_cells, 1.0);

int    pml_size  = 10;
double m         = 3.0;
double a_max     = 0.05;
double kappa_max = 5.0;

vector<double> sigma_pml(pml_size, 0.0);
vector<double> sigma_h_pml(pml_size, 0.0);
vector<double> kappa_far(pml_size, 0.0);
vector<double> kappa_h_far(pml_size, 0.0);
vector<double> a_pml(pml_size, 0.0);
vector<double> a_h_pml(pml_size, 0.0);
vector<double> kappa_left(pml_size, 0.0);
vector<double> kappa_h_left(pml_size, 0.0);

vector<double> be_y_f(pml_size, 0.0);
vector<double> ce_y_f(pml_size, 0.0);
vector<double> bh_y_f(pml_size, 0.0);
vector<double> ch_y_f(pml_size, 0.0);

vector<double> be_y(pml_size, 0.0);
vector<double> ce_y(pml_size, 0.0);
vector<double> bh_y(pml_size, 0.0);
vector<double> ch_y(pml_size, 0.0);

double signal(double time, double pulse_width, double pulse_delay, double omega0, double amplitude) {
    double envelope = exp(-pow(((time - pulse_delay) / pulse_width), 2));
    return envelope * amplitude * sin(omega0 * time);
}

int main() {
    
    // Base material & stretched-coordinate arrays
    for (int i = 0; i < N_space_cells; i++) {
        refractive_index[i] = sqrt(eps[i]);
        inv_kappa_dy[i]    /= step_size;     // 1/dy
        inv_kappa_h_dy[i]  /= step_size;     // 1/dy
    }

    // NOTE: match Python: log(1e-10), not 1e-9
    double n_edge     = refractive_index[N_space_cells - 1 - pml_size];
    double pml_cond_e = -(m + 1.0) * log(1e-10) * c0 * EPSILON_0 / (2.0 * n_edge * pml_size * step_size);

    // ---------- Right (+y) PML ----------
    for (int i = 0; i < pml_size; i++) {
        double x = static_cast<double>(i) / pml_size; // 0 -> 1
        sigma_pml[i]   = pml_cond_e * pow(x, m);
        kappa_far[i]   = 1.0 + (kappa_max - 1.0) * pow(x, m);
        a_pml[i]       = a_max * ((pml_size - i) / static_cast<double>(pml_size));

        sigma_h_pml[i] = pml_cond_e * pow((i + 0.5) / pml_size, m);
        kappa_h_far[i] = 1.0 + (kappa_max - 1.0) * pow((i + 0.5) / pml_size, m);
        a_h_pml[i]     = a_max * ((pml_size - i + 0.5) / static_cast<double>(pml_size));

        be_y_f[i] = exp(-(sigma_pml[i] / kappa_far[i] + a_pml[i]) * dt / EPSILON_0);
        ce_y_f[i] = sigma_pml[i] * (be_y_f[i] - 1.0) / ((sigma_pml[i] + kappa_far[i] * a_pml[i]) * kappa_far[i]);

        bh_y_f[i] = exp(-(sigma_h_pml[i] / kappa_h_far[i] + a_h_pml[i]) * dt / EPSILON_0);
        ch_y_f[i] = sigma_h_pml[i] * (bh_y_f[i] - 1.0) / ((sigma_h_pml[i] + kappa_h_far[i] * a_h_pml[i]) * kappa_h_far[i]);

        // Right edge inverse-kappa (fix off-by-one: start at N - pml_size)
        inv_kappa_dy[N_space_cells - pml_size + i]    = 1.0 / (kappa_far[i] * step_size);
        inv_kappa_h_dy[N_space_cells - pml_size + i]  = 1.0 / (kappa_h_far[i] * step_size);
    }

    // ---------- Left (−y) PML ----------
    for (int i = 0; i < pml_size; i++) {
        double x = static_cast<double>(pml_size - i) / pml_size; // 1 -> 1/pml_size
        sigma_pml[i]   = pml_cond_e * pow(x, m);
        kappa_left[i]  = 1.0 + (kappa_max - 1.0) * pow(x, m);
        a_pml[i]       = a_max * (static_cast<double>(i) / pml_size); // start at 0, not (i+1)/pml_size

        sigma_h_pml[i] = pml_cond_e * pow((pml_size - i - 0.5) / pml_size, m);
        kappa_h_left[i]= 1.0 + (kappa_max - 1.0) * pow((pml_size - i - 0.5) / pml_size, m);
        a_h_pml[i]     = a_max * ((i + 0.5) / pml_size);

        be_y[i] = exp(-(sigma_pml[i] / kappa_left[i] + a_pml[i]) * dt / EPSILON_0);
        ce_y[i] = sigma_pml[i] * (be_y[i] - 1.0) / ((sigma_pml[i] + kappa_left[i] * a_pml[i]) * kappa_left[i]);

        bh_y[i] = exp(-(sigma_h_pml[i] / kappa_h_left[i] + a_h_pml[i]) * dt / EPSILON_0);
        ch_y[i] = sigma_h_pml[i] * (bh_y[i] - 1.0) / ((sigma_h_pml[i] + kappa_h_left[i] * a_h_pml[i]) * kappa_h_left[i]);

        inv_kappa_dy[i]   = 1.0 / (kappa_left[i] * step_size);
        inv_kappa_h_dy[i] = 1.0 / (kappa_h_left[i] * step_size);
    }

    // Optional: write sigma_e/h (will be zeros unless you set dispersive media)
    {
        std::ofstream sigma_out("sigma_e.txt");
        std::ofstream sigma_h_out("sigma_h.txt");
        for (int i = 0; i < N_space_cells; i++) {
            sigma_out   << sigma[i]   << '\n';
            sigma_h_out << sigma_h[i] << '\n';
        }
    }

    // Field update coefficients
    vector<double> denominator(N_space_cells, 1.0);
    vector<double> e_coeff1(N_space_cells, 1.0);
    vector<double> e_coeff2(N_space_cells, 1.0);

    vector<double> denominator_h(N_space_cells, 1.0);
    vector<double> h_coeff1(N_space_cells, 1.0);
    vector<double> h_coeff2(N_space_cells, 1.0);

    for (int i = 0; i < N_space_cells; i++) {
        denominator[i]  = EPSILON_0 * eps[i] / dt + sigma[i] / 2.0;
        e_coeff1[i]     = (EPSILON_0 * eps[i] / dt - sigma[i] / 2.0) / denominator[i];
        e_coeff2[i]     = 1.0 / denominator[i];

        denominator_h[i]= MU_0 / dt + sigma_h[i] / 2.0;
        h_coeff1[i]     = (MU_0 / dt - sigma_h[i] / 2.0) / denominator_h[i];
        h_coeff2[i]     = 1.0 / denominator_h[i];
    }

    // Source setup
    
    double omega0            = 2.0 * M_PI * c0 / center_wavelength;
    double pulse_width       = 10e-15;
    double pulse_delay       = 4.0 * pulse_width;

	cout<<"Simulation size = "<<simulation_size<<endl;
    cout<<"Simulation time = "<<simulation_time<<endl;
    cout<<"Step Size = "<<step_size<<endl;
    cout<<"Pulse Width = "<<pulse_width<<endl;
    cout<<"Pulse delay = "<<pulse_delay<<endl;
    cout<<"center_wavelength = "<<center_wavelength<<endl;
    // exit(0);

    int    jsource  = int(N_space_cells / 2);
    double t_offset = refractive_index[jsource] * step_size / (2.0 * c0);
    double Z        = imp0 / refractive_index[jsource];
    double amplitude= 1.0;

    // Time stepping
    for (int n = 0; n < N_time_steps; n++) {
        // Save previous fields
        for (int j = 0; j < N_space_cells; j++) {
            Hz_prev[j] = Hz[j];
            Ex_prev[j] = Ex[j];
        }

        // --- Update Magnetic CPML fields ---
        // Right (+y): indices [N - pml_size, N - 2] like Python [-pml_size:-1]
        for (int i = N_space_cells - pml_size; i < N_space_cells - 1; i++) {
            int idx     = i - (N_space_cells - pml_size);
            P_Hz_y[i]   = bh_y_f[idx] * P_Hz_y[i]
                        + ch_y_f[idx] * (Ex[i + 1] - Ex[i]) / step_size;
        }
        // Left (−y): indices [0, pml_size - 1]
        for (int i = 0; i < pml_size; i++) {
            P_Hz_y[i]   = bh_y[i] * P_Hz_y[i]
                        + ch_y[i] * (Ex[i + 1] - Ex[i]) / step_size;
        }

        // --- Update Magnetic Field Hz ---
        for (int j = 0; j < N_space_cells - 1; j++) {
            Hz[j] = h_coeff1[j] * Hz_prev[j]
                  + h_coeff2[j] * (inv_kappa_h_dy[j] * (Ex[j + 1] - Ex[j]) + P_Hz_y[j]);
        }

        // (Optional magnetic source—commented to match your current choice)
        // Hz[jsource - 1] -= signal((n + 0.5) * dt - t_offset, pulse_width, pulse_delay, omega0, amplitude) / Z;
        // Hz[jsource - 1] += signal((n + 0.5) * dt + t_offset, pulse_width, pulse_delay, omega0, amplitude) / Z;

        // --- Update Electric CPML fields ---
        // Right (+y): indices [N - pml_size, N - 1]
        for (int i = N_space_cells - pml_size; i < N_space_cells; i++) {
            int idx     = i - (N_space_cells - pml_size);
            P_Ex_y[i]   = be_y_f[idx] * P_Ex_y[i]
                        + ce_y_f[idx] * (Hz[i] - Hz[i - 1]) / step_size;
        }
        // Left (−y): indices [1, pml_size - 1] like Python 1:pml_size
        for (int i = 1; i < pml_size; i++) {
            P_Ex_y[i]   = be_y[i] * P_Ex_y[i]
                        + ce_y[i] * (Hz[i] - Hz[i - 1]) / step_size;
        }

        // --- Update Electric Field Ex ---
        for (int j = 1; j < N_space_cells - 1; j++) {
            Ex[j] = e_coeff1[j] * Ex_prev[j]
                  + e_coeff2[j] * (inv_kappa_dy[j] * (Hz[j] - Hz[j - 1]) + P_Ex_y[j]);
        }

        // Soft electric source at jsource
        Ex[jsource] += signal((n + 1) * dt, pulse_width, pulse_delay, omega0, amplitude);

        // --- Diagnostics ---
        if (n % 100 == 0) {
            double max_val = *std::max_element(Ex.begin(), Ex.end());
            double min_val = *std::min_element(Ex.begin(), Ex.end());
            cout << "Time step : " << n << "\n";
            cout << "Emin = " << min_val << ", Emax = " << max_val << "\n";

            std::ofstream exf("data/Ex/ElectricField_" + to_string(n) + ".txt");
            for (int i = 0; i < N_space_cells; i++) exf << Ex[i] << '\n';
        }
    }

    return 0;
}
