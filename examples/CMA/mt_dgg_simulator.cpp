#include <iostream>

#include "DggFactory.hpp"
#include "PlantModel.hpp"
#include "cmaModel.h"
#include "simdjson.h"

int main()
{   

    std::cout << "Running Microtubule Dynamic Graph Grammar Simulator\n";

    // The user defines the parser, loads, and parses configuration file.
    simdjson::ondemand::parser parser;
    simdjson::padded_string json = simdjson::padded_string::load("settings.json");
    simdjson::ondemand::document settings_file = parser.iterate(json);

//    DGGML::DggFactory<DGGML::PlantModel, simdjson::ondemand::document> plant_factory;
//
//    auto num_simulations = 1;
//
//    //TODO: add some sort of ensemble simulation mode
//    for(auto i = 0; i < num_simulations; i++)
//        plant_factory.execute(settings_file);

    // The user builds their model
    // I think I want something to be like this:
//    using GT = Graph<TypeA, TypeB>;
//    GT g1 = {{Node<TypeA>{0,...params}, Node<TypeB>{1,0,...params}, Node<TypeA>{2,0,...params}}, {Edge(0,1), Edge(1,2), Edge(2,0)}};
//
//    Grammar gamma;
//
//    gamma.lhs(g1).rhs(g2)
//            .with([](GT& lhs) -> double
//                  {
//                      return 2.0;
//                  })
//            .where([](const GT& lhs, GT& rhs) -> void
//                   {
//                       auto lnode1 = std::get<TypeA>(lhs[0]);
//                       auto rnode1 = std::get<TypeA>(rhs[0]);
//                       rnode1.x = lnode1.x;
//                   });
//
//    Simulator experiment1(gamma, cell_complex, initial_state);
//    experiment1.run();

    using GT = DGGML::Plant::graph_type;
    DGGML::Grammar<GT> gamma;
    GT g1, g2;
    DGGML::WithRule<GT> r1;

    r1.name("test rule").lhs(g1).rhs(g2)
    .with([](const GT& l) -> double { return 2.0;})
    .where([](const GT& l, GT& r) -> void {});

    gamma.addRule(r1);

    std::cout << "P: " << gamma.stochastic_rules["test rule"].propensity(g1) << "\n";

    DGGML::SimulatorInterface<CMA::cmaModel> cma_simulation;

    CMA::cmaModel experiment1;
    experiment1.set_parameters(settings_file);
    cma_simulation.setModel(experiment1);
    cma_simulation.simulate();
}
