#include <iostream>

#include "catch.hpp"

#include "../src/DggFactory.hpp"
#include "../src/PlantModel.hpp"

struct MyInterface {};

TEST_CASE("Running Factory Test", "[factory test]")
{    
    std::cout << "Running factory test (please run with a memory analysis tool)\n";

    MyInterface interface;

    Cajete::DggFactory<Cajete::PlantModel, MyInterface> plant_factory;
    
    auto num_simulations = 3;
    
    for(auto i = 0; i < num_simulations; i++)
        plant_factory.execute(interface);
}
