#include <iostream>

#include "DggFactory.hpp"
#include "PlantModel.hpp"
#include "cmaModel.h"
#include "simdjson.h"

int main()
{   

    std::cout << "Running Microtubule Dynamic Graph Grammar Simulator\n";

    // The user defines the parser, loads, and parses configuration file.
    //simdjson::ondemand::parser parser;
    //simdjson::padded_string json = simdjson::padded_string::load("settings.json");
    //simdjson::ondemand::document settings_file = parser.iterate(json);

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

    Plant::graph_type g3;
    g3.addNode({1, {Plant::Negative{}}});
    g3.addNode({2, {Plant::Intermediate{}}});
    g3.addNode({3, {Plant::Positive{}}});
    g3.addEdge(1, 2);
    g3.addEdge(2, 3);

    for(double& x : g3[1].position) x = 2.2;
    for(double & i : g3[1].getData<Plant::Negative>().velocity) i = 2.5;
    for(double & i : g3[2].getData<Plant::Intermediate>().velocity) i = 1.0;

    Plant::graph_type g4;
    g4.addNode({9, {Plant::Negative{}}});
    g4.addNode({10, {Plant::Intermediate{}}});
    g4.addNode({11, {Plant::Positive{}}});
    g4.addEdge(9, 10);
    g4.addEdge(10, 11);

    for(double& x : g4[9].position) x = 39.6;
    for(double & x : g4[9].getData<Plant::Negative>().velocity) x = 5.5;
    for(double & x : g4[10].getData<Plant::Intermediate>().velocity) x = 6.0;

    auto result = YAGL::subgraph_isomorphism(g3, g4);

    for(auto& item : result)
        for(auto& m : item)
            std::cout << m.first << " -> " << m.second << "\n";

    std::cout << result[0][1] << "\n";
    using GT = Plant::graph_type;
    DGGML::Grammar<GT> gamma;
    using GMT = DGGML::WithRule<GT>::GraphMapType;
    GT g1, g2;
    DGGML::WithRule<GT> r1;

    r1.name("test rule").lhs(g3).rhs(g4)
    .with([](GT& l, GMT& m) -> double {
        return l[m[1]].getData<Plant::Negative>().velocity[0];
    })
    .where([](GT& l, GT& r, GMT& m) -> void { for(double& x : r[m[1]].getData<Plant::Negative>().velocity) x = 1.3; });

    gamma.addRule(r1);

    std::cout << "P: " << gamma.stochastic_rules["test rule"].propensity(g4, result[0]) << "\n";
    std::cout << "Before V: " << g4[9].getData<Plant::Negative>().velocity[0] << "\n";
    gamma.stochastic_rules["test rule"].update(g3, g4, result[0]);
    std::cout << "Updated V: " << g4[9].getData<Plant::Negative>().velocity[0] << "\n";

    DGGML::SimulatorInterface<CMA::cmaModel> cma_simulation;
    CMA::cmaModel experiment1;
    //experiment1.set_parameters(settings_file);
    experiment1.initialize();
    experiment1.gamma.print();
    cma_simulation.setModel(experiment1);
    cma_simulation.simulate();
}




