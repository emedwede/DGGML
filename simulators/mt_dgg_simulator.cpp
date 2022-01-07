#include <iostream>

#include "DggFactory.hpp"
#include "PlantModel.hpp"

struct MyInterface {};

int main()
{   

    std::cout << "Running Microtubule Dynamic Graph Grammar Simulator\n";

    MyInterface interface;

    Cajete::DggFactory<Cajete::PlantModel, MyInterface> plant_factory;
    
    auto num_simulations = 1;
    
    for(auto i = 0; i < num_simulations; i++)
        plant_factory.execute(interface);
}
