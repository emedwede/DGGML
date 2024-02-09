#include <iostream>
#include <chrono>

#include "DggFactory.hpp"
#include "particleModel.h"
#include "simdjson.h"

int main()
{   

    std::cout << "Running the Particle Dynamic Graph Grammar Simulator\n";

    // The user defines the parser, loads, and parses configuration file.
    //simdjson::ondemand::parser parser;
    //simdjson::padded_string json = simdjson::padded_string::load("settings.json");
    //simdjson::ondemand::document settings_file = parser.iterate(json);


    auto start = std::chrono::high_resolution_clock::now();
    DGGML::SimulatorInterface<MD::particleModel> particle_simulation;
    MD::particleModel experiment1;
    //experiment1.set_parameters(settings_file);
    particle_simulation.setModel(experiment1);
    particle_simulation.simulate();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    std::cout << "\n\nSimulation took " << duration.count() / 1000.0 << " seconds\n";

}




