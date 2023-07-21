#include <iostream>

#include "catch.hpp"

#include "DggFactory.hpp"

#include "simdjson.h"

TEST_CASE("Running Factory Test", "[factory test]")
{    
    std::cout << "Running factory test (please run with a memory analysis tool)\n";
    
    simdjson::ondemand::parser parser;
    simdjson::padded_string json = simdjson::padded_string::load("settings.json");
    simdjson::ondemand::document settings_file = parser.iterate(json);

    DGGML::DggFactory<DGGML::PlantModel, simdjson::ondemand::document> plant_factory;
    
    auto num_simulations = 3;
    
    for(auto i = 0; i < num_simulations; i++)
        plant_factory.execute(settings_file);
}
