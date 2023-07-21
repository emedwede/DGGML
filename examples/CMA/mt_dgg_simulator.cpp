#include <iostream>

#include "DggFactory.hpp"
#include "PlantModel.hpp"
#include "cmaModel.h"
#include "simdjson.h"

int main()
{   

    std::cout << "Running Microtubule Dynamic Graph Grammar Simulator\n";

    simdjson::ondemand::parser parser;
    simdjson::padded_string json = simdjson::padded_string::load("settings.json");
    simdjson::ondemand::document settings_file = parser.iterate(json);

    DGGML::DggFactory<DGGML::PlantModel, simdjson::ondemand::document> plant_factory;
    
    auto num_simulations = 1;
   
    //TODO: add some sort of ensemble simulation mode
    for(auto i = 0; i < num_simulations; i++)
        plant_factory.execute(settings_file);

    DGGML::SimulatorInterface<CMA::cmaModel> cma_simulation;

    CMA::cmaModel experiment1;
    experiment1.set_parameters(settings_file);
    experiment1.initialize();
    experiment1.gamma.print();
    experiment1.name = "experiment1";


    cma_simulation.setModel(experiment1);
    cma_simulation.simulate();
}
