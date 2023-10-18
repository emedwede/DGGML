#include <random>

#include "catch.hpp"
#include "VtkWriter.hpp"
#include "Grammar.h"
#include "AnalyzedGrammar.hpp"
#include "YAGL_Algorithms.hpp"
#include "ExpandedComplex2D.hpp"
#include "RuleSystem.hpp"
#include "MathUtils.hpp"

#include "../examples/CMA/PlantModel/PlantTypes.hpp"


using GraphType = Plant::graph_type;
using KeyType = GraphType::key_type;

void define_model(DGGML::Grammar<GraphType>& gamma) {

    GraphType g1;
    g1.addNode({1, {Plant::Intermediate{}}});
    g1.addNode({2, {Plant::Positive{}}});
    g1.addEdge(1, 2);

    GraphType g2;
    g2.addNode({1, {Plant::Intermediate{}}});
    g2.addNode({3, {Plant::Intermediate{}}});
    g2.addNode({2, {Plant::Positive{}}});
    g2.addEdge(1, 3);
    g2.addEdge(3, 2);

    DGGML::WithRule<GraphType> r1("with_growth", g1, g2,
                                  [](auto& lhs, auto& m) { return 0.0; },
                                  [](auto& lhs, auto& rhs, auto& m1, auto& m2) {
                                      rhs[m2[3]].position[0] = (lhs[m1[2]].position[0] + lhs[m1[1]].position[0])/2.0;
                                      rhs[m2[3]].position[1] = (lhs[m1[2]].position[1] + lhs[m1[1]].position[1])/2.0;
                                      rhs[m2[3]].position[2] = (lhs[m1[2]].position[2] + lhs[m1[1]].position[2])/2.0;
                                  });
    gamma.addRule(r1);

    GraphType g3;
    g3.addNode({1, {Plant::Negative{}}});
    g3.addNode({2, {Plant::Intermediate{}}});
    g3.addNode({3, {Plant::Positive{}}});
    g3.addEdge(1, 2);
    g3.addEdge(2, 3);

    g3.addNode({4, {Plant::Negative{}}});
    g3.addNode({5, {Plant::Intermediate{}}});
    g3.addNode({6, {Plant::Positive{}}});
    g3.addEdge(4, 5);
    g3.addEdge(5, 6);

    GraphType g4;


    DGGML::WithRule<GraphType> r2("with_interaction", g3, g4,
                                  [](auto& lhs, auto& m) { return 0.0; },
                                  [](auto& lhs, auto& rhs, auto& m1, auto& m2) {});
    gamma.addRule(r2);
}

void initialize_system(GraphType& system_graph)
{
    for(auto i = 0; i < 4; i++) {
        std::size_t k1 = i*3;
        std::size_t k2 = k1+1;
        std::size_t k3 = k1+2;
        //Add two microtubule segments parallel to each other
        system_graph.addNode({k1, {Plant::Negative{0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0},
                                    i+0.1, i+0.1, 0.0}});
        system_graph.addNode({k2, {Plant::Intermediate{0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0},
                                    i+0.5, i+0.1, 0.0}});
        system_graph.addNode({k3, {Plant::Positive{0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0},
                                    i+0.9, i+0.1, 0.0}});
        system_graph.addEdge(k1, k2);
        system_graph.addEdge(k2, k3);
    }
}

TEST_CASE("Phi No Decomposition", "[Phi Test]")
{
    DGGML::Grammar<GraphType> gamma;
    define_model(gamma);
    DGGML::AnalyzedGrammar<GraphType> gamma_analysis(gamma);
    DGGML::ExpandedComplex2D<> grid(1, 1, 8.0, 8.0, false, 1.0);
    DGGML::GridFileWriter grid_writer;
    grid_writer.save({grid.reaction_grid, grid.dim_label}, "test_grid");

    GraphType system_graph;
    initialize_system(system_graph);
    DGGML::VtkFileWriter<GraphType> vtk_writer;
    vtk_writer.save(system_graph, "phi_graph");

    //find matches
    auto c1 = gamma_analysis.unique_components[0];
    auto end_matches = YAGL::subgraph_isomorphism2(c1, system_graph);
    REQUIRE(end_matches.size() == 4);
    auto c2 = gamma_analysis.unique_components[1];
    auto mt_matches = YAGL::subgraph_isomorphism2(c2, system_graph);
    REQUIRE(mt_matches.size() == 4);

    std::map<std::size_t, DGGML::Instance<std::size_t>> component_match_set;
    int k = 0;
    for(auto& item : end_matches)
    {
        std::vector<std::size_t> match;
        for(auto& [key, value] : item)
        {
            match.emplace_back(value);
        }
        DGGML::Instance<KeyType> inst;
        inst.match = match;
        inst.type = 0;
        inst.anchor = match[0];
        component_match_set.insert({k++, inst});
        //rule_system.push_back(inst);
    }
    for(auto& item : mt_matches)
    {
        std::vector<std::size_t> match;
        for(auto& [key, value] : item)
        {
            match.emplace_back(value);
        }
        DGGML::Instance<KeyType> inst;
        inst.match = match;
        inst.type = 1;
        inst.anchor = match[0];
        component_match_set.insert({k++, inst});
        //rule_system.push_back(inst);
    }


    REQUIRE(component_match_set.size() == 8);

    //building the rule instances
    int kk = 0;
    std::map<std::size_t, DGGML::RuleInstType<std::size_t>> rule_instances;
    for(int i = 0; i < component_match_set.size(); i++)
    {
        if(component_match_set[i].type == 0) {
            rule_instances[kk];
            rule_instances[kk].name = "with_growth";
            rule_instances[kk].category = "stochastic";
            rule_instances[kk].components.push_back(i);
            rule_instances[kk].anchor = component_match_set[i].anchor;
            kk++;
        }
        for(int j = 0; j < component_match_set.size(); j++)
        {
            if(j <= i ) continue;
            auto& comp1 = component_match_set[i];
            auto& comp2 = component_match_set[j];
            auto& p1 = system_graph.findNode(comp1.anchor)->second.getData().position;
            auto& p2 = system_graph.findNode(comp2.anchor)->second.getData().position;
            auto d = DGGML::calculate_distance(p1, p2);
            if(comp1.type == 1 && comp2.type == 1 && d <= 2.0)
            {
                rule_instances[kk];
                rule_instances[kk].name = "with_interaction";
                rule_instances[kk].category = "stochastic";
                rule_instances[kk].components.push_back(i);
                rule_instances[kk].components.push_back(j);
                rule_instances[kk].anchor = comp1.anchor;
                kk++;
            }
        }
    }
    //TODO: fix this part
    //REQUIRE(rule_instances.size() == 4);

    for(auto& [key, inst] : rule_instances)
    {
        auto& p = system_graph.findNode(inst.anchor)->second.getData().position;
        int ic, jc;
        grid.reaction_grid.locatePoint(p[0], p[1], ic, jc);
        auto cardinal = grid.reaction_grid.cardinalCellIndex(ic, jc);
        // here would be the anchor list reduction, plus whatever other code needs to be added
        //anchor_list.insert({match.first.anchor, cardinal});
        auto max_cell = grid.cell_label[cardinal];
        //rule_map[max_cell].push_back(key);
        std::cout << "rule is mapped to cell " << max_cell << "\n";
    }


}

TEST_CASE("Phi Decomposition", "[Phi Test]")
{
    //TODO: fix epsilon here, causing errors
    //DGGML::ExpandedComplex2D<> grid(2, 2, 4.0, 4.0, false, 1.0);
    DGGML::GridFileWriter grid_writer;
    //grid_writer.save({grid.reaction_grid, grid.dim_label}, "test_grid");
    //for(auto& item : grid.dim_label) std::cout << item << "\n";
}