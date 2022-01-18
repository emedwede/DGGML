#include <iostream>

#include "DggFactory.hpp"
#include "PlantModel.hpp"

#include "simdjson.h"

struct MyInterface {};

int main()
{   

    std::cout << "Running Microtubule Dynamic Graph Grammar Simulator\n";

    MyInterface interface;

    simdjson::ondemand::parser parser;
    simdjson::padded_string json = simdjson::padded_string::load("settings.json");
    simdjson::ondemand::document settings_file = parser.iterate(json);

    std::cout << uint64_t(settings_file["SETTINGS"]["NUM_STEPS"]) << " number of steps\n";
    Cajete::DggFactory<Cajete::PlantModel, simdjson::ondemand::document> plant_factory;
    
    auto num_simulations = 1;
    
    for(auto i = 0; i < num_simulations; i++)
        plant_factory.execute(settings_file);
}
