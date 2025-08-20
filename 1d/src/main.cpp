#include "fdtd_1d.h"
#include <json/json.h>
#include <fstream>
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " control_file.json" << std::endl;
        return 1;
    }

    std::string control_file = argv[1];

    // Load JSON configuration
    Json::Value root;
    std::ifstream file(control_file, std::ifstream::binary);
    if (!file.is_open()) {
        std::cerr << "Error: could not open control file " << control_file << std::endl;
        return 1;
    }
    file >> root;

    // Initialize and run simulation
    FDTD_1D sim(control_file);   // pass string, not Json::Value

    sim.run();

    std::cout << "Simulation finished successfully!" << std::endl;
    return 0;
}
