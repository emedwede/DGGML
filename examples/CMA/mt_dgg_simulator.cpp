#include <iostream>
#include <chrono>

#include "DggFactory.hpp"
#include "cmaModel.h"
#include "simdjson.h"

int main(int argc, char* argv[])
{

    // Check if filename argument is provided
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }
    std::cout << "Running Microtubule Dynamic Graph Grammar Simulator\n";

    std::string filename = argv[1];
    // The user defines the parser, loads, and parses configuration file.
    simdjson::ondemand::parser parser;
    simdjson::padded_string json = simdjson::padded_string::load(filename);
    simdjson::ondemand::document settings_file = parser.iterate(json);

    auto start = std::chrono::high_resolution_clock::now();
    DGGML::SimulatorInterface<CMA::cmaModel> cma_simulation;
    CMA::cmaModel experiment1;
    experiment1.set_parameters(settings_file);
    //return 0;
    cma_simulation.setModel(experiment1);
    cma_simulation.simulate();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    std::cout << "\n\nSimulation took " << duration.count() / 1000.0 << " seconds\n";

}




