#include <random>
#include <list>
#include <algorithm>
#include <numeric>

#include "catch.hpp"
#include "VtkWriter.hpp"
#include "Grammar.h"
#include "AnalyzedGrammar.hpp"
#include "YAGL_Algorithms.hpp"
#include "ExpandedComplex2D.hpp"
#include "ComponentMap.hpp"
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

//creates a random mt centered in a cell with keys k, k+1, k+2
void create_mt(int ic, int jc, std::size_t k, GraphType& system_graph)
{
    double epsilon_min = 0.35;
    double epsilon_max = 0.45;
    std::random_device random_device;
    std::mt19937 random_engine(random_device());
    std::uniform_real_distribution<double> distribution_angle(0.0, 2.0*3.14);
    std::uniform_real_distribution<double> distribution_local(epsilon_min, epsilon_max);

    double x_c, y_c;
    double z_c = 0;
    x_c = double(ic)+0.5;
    y_c = double(jc)+0.5;
    auto theta = distribution_angle(random_engine);
    auto seg_len = distribution_local(random_engine);

    auto x_s = 0.0;
    auto y_s = seg_len;
    auto x_r_t = x_s*cos(theta) + y_s*sin(theta);
    auto y_r_t = -x_s*sin(theta) +y_s*cos(theta);
    auto x_r = x_c + x_r_t;
    auto y_r = y_c + y_r_t;
    auto z_r = 0.0;

    auto x_l = x_c - (x_r - x_c);
    auto y_l = y_c - (y_r - y_c);
    auto z_l = 0.0;

    //compute dist and unit vector
    double p1[3] = {x_r - x_c, y_r - y_c, 0.0};
    double p2[3] = {0.0, 0.0, 0.0};
    double p3[3] = {x_c - x_l, y_c - y_l, 0.0};

    system_graph.addNode({k, {Plant::Negative{0.0, 0.0, 0.0,
                                               0.0, 0.0, 0.0},
                               x_l, y_l, z_l}});
    system_graph.addNode({k+1, {Plant::Intermediate{0.0, 0.0, 0.0,
                                              0.0, 0.0, 0.0},
                              x_c, y_c, z_c}});
    system_graph.addNode({k+2, {Plant::Positive{0.0, 0.0, 0.0,
                                              0.0, 0.0, 0.0},
                              x_r, y_r, z_r}});
    system_graph.addEdge(k, k+1);
    system_graph.addEdge(k+1, k+2);
}

//This initializes a test system hand created, MTs are oriented randomly
//within a reaction cell, but it shouldn't change the results for the reactions
//mapped to a particular geocell
void initialize_system(GraphType& system_graph)
{
    create_mt(6, 7, 1, system_graph);
    create_mt(3, 6, 4, system_graph);
    create_mt(6, 6, 7, system_graph);
    create_mt(0, 4, 10, system_graph);
    create_mt(1, 4, 13, system_graph);
    create_mt(2, 4, 16, system_graph);
    create_mt(3, 4, 19, system_graph);
    create_mt(5, 3, 22, system_graph);
    create_mt(7, 3, 25, system_graph);
    create_mt(1, 2, 28, system_graph);
    //manual placement for this one
    system_graph.addNode({31, {Plant::Negative{0.0, 0.0, 0.0,
                                              0.0, 0.0, 0.0},
                               2.5, 2.5, 0.0}});
    system_graph.addNode({32, {Plant::Intermediate{0.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0},
                                2.5, 2.1, 0.0}});
    system_graph.addNode({33, {Plant::Positive{0.0, 0.0, 0.0,
                                                0.0, 0.0, 0.0},
                                2.5, 1.7, 0.0}});
    system_graph.addEdge(31, 32);
    system_graph.addEdge(32, 33);
    create_mt(7, 2, 34, system_graph);
    create_mt(6, 1, 37, system_graph);
    create_mt(5, 0, 40, system_graph);
}

TEST_CASE("Graph sort", "[Graph sort]")
{
    std::random_device random_device;
    std::mt19937 random_engine(random_device());
    std::uniform_real_distribution<double> distribution_real(0.0, 10.0);
    std::uniform_real_distribution<double> distribution_angle(0.0, 2.0*3.14);
    std::uniform_real_distribution<double> distribution_local(0.5, 0.7);

    YAGL::Graph<std::size_t, SpatialNode3D<int>> g;
    for(auto i = 0; i < 20; i+=2) {
        double x1 = distribution_real(random_engine);
        double x2 = distribution_real(random_engine);
        double theta = distribution_angle(random_engine);
        double length = distribution_local(random_engine);
        auto x3 = x1+(0.0*cos(theta) + length*sin(theta));
        auto x4 = x2+(-0.0*sin(theta) +length*cos(theta));
        g.addNode({i, {{0}, {x1, x2, 0.0}}});
        g.addNode({i+1, {{0}, {x3, x4, 0.0}}});
        g.addEdge(i, i+1);
        auto d = sqrt((x3-x1)*(x3-x1)+(x4-x2)*(x4-x2));
        std::cout << d << "\n";
    }

    for(auto& [k, v] : g.getNodeSetRef())
    {
        std::cout << k << ": { " << v.getData().position[0] << ", " << v.getData().position[1] << " }\n";
    }
    std::cout << YAGL::connected_components(g) << "\n";
}
TEST_CASE("Phi No Decomposition", "[Phi Test]")
{
    DGGML::Grammar<GraphType> gamma;
    define_model(gamma);
    DGGML::AnalyzedGrammar<GraphType> gamma_analysis(gamma);
    DGGML::ExpandedComplex2D<> grid(2, 2, 12.0, 12.0, false, 1.0);
    DGGML::GridFileWriter grid_writer;
    grid_writer.save({grid.reaction_grid, grid.dim_label}, "test_grid");

    GraphType system_graph;
    initialize_system(system_graph);
    REQUIRE(system_graph.numNodes() == 42);
    REQUIRE(YAGL::connected_components(system_graph) == 14);
    DGGML::VtkFileWriter<GraphType> vtk_writer;
    vtk_writer.save(system_graph, "phi_graph");

    //find matches
    auto c1 = gamma_analysis.unique_components[0];
    auto end_matches = YAGL::subgraph_isomorphism2(c1, system_graph);
    REQUIRE(end_matches.size() == 14);
    auto c2 = gamma_analysis.unique_components[1];
    auto mt_matches = YAGL::subgraph_isomorphism2(c2, system_graph);
    REQUIRE(mt_matches.size() == 14);

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

    int kk = 0;
    std::map<std::size_t, DGGML::RuleInstType<std::size_t>> rule_instances;
    for(int i = 0; i < component_match_set.size(); i++) {
        if (component_match_set[i].type == 0) {
            rule_instances[kk];
            rule_instances[kk].name = "with_growth";
            rule_instances[kk].category = "stochastic";
            rule_instances[kk].components.push_back(i);
            rule_instances[kk].anchor = component_match_set[i].anchor;
            kk++;
        }
    }

    //should find first 14 growing ends
    REQUIRE(rule_instances.size() == 14);

    //TODO: this should really be a cell list class, to reduce boilerplate code
    std::vector<std::list<std::size_t>> cell_list;
    cell_list.resize(grid.reaction_grid.totalNumCells());
    for(auto& [key, inst] : component_match_set)
    {
        auto& p = system_graph.findNode(inst.anchor)->second.getData().position;
        int ic, jc;
        grid.reaction_grid.locatePoint(p[0], p[1], ic, jc);
        auto c = grid.reaction_grid.cardinalCellIndex(ic, jc);
        cell_list[c].push_back(key);
    }

    std::size_t sum = 0.0; for(auto& cell : cell_list) sum += cell.size();
    std::cout << sum << "\n";

    for(auto i = 0; i < grid.reaction_grid._nx; i++)
    {
        for(auto j = 0; j < grid.reaction_grid._ny; j++) {
            auto c = grid.reaction_grid.cardinalCellIndex(i, j);
            for (auto &k1: cell_list[c]) {
                //TODO: there may be issues because this formulation does not ignore
                //ghosted reaction cells
                auto cell_range = 1;
                auto imin = (i - cell_range > 0) ? i - cell_range : 0;
                auto imax = (i + cell_range + 1) < grid.reaction_grid._nx ? i + cell_range + 1 : grid.reaction_grid._nx;

                auto jmin = (j - cell_range > 0) ? j - cell_range : 0;
                auto jmax = (j + cell_range + 1) < grid.reaction_grid._ny ? j + cell_range + 1 : grid.reaction_grid._ny;
                for (auto ii = imin; ii < imax; ii++) {
                    for (auto jj = jmin; jj < jmax; jj++) {
                        auto n = grid.reaction_grid.cardinalCellIndex(ii, jj);
                        if (n != c) {
                            for(auto& k2 : cell_list[n])
                            {
                                if(k1 == k2) continue;
                                auto& comp1 = component_match_set[k1];
                                auto& comp2 = component_match_set[k2];
                                //auto& p1 = system_graph.findNode(comp1.anchor)->second.getData().position;
                                //auto& p2 = system_graph.findNode(comp2.anchor)->second.getData().position;
                                //auto d = DGGML::calculate_distance(p1, p2);
                                if(comp1.type == 1 && comp2.type == 1) {//&& d <= 1.5) {
                                    rule_instances[kk];
                                    rule_instances[kk].name = "with_interaction";
                                    rule_instances[kk].category = "stochastic";
                                    rule_instances[kk].components.push_back(k1);
                                    rule_instances[kk].components.push_back(k2);
                                    rule_instances[kk].anchor = comp1.anchor;
                                    kk++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //should find 16 more pair wise interaction rule instances
    //the number would be different if we imposed a reaction radius
    //here we just search neighboring cells
    REQUIRE(rule_instances.size() == 30);

    using cplex_key_t = typename DGGML::CartesianComplex2D<>::graph_type::key_type;
    std::map<cplex_key_t, std::vector<std::size_t>> rule_map;

    //building rule map sets
    for(auto& [key, value] : grid.graph.getNodeSetRef()) {
        if(value.getData().interior)
            rule_map.insert({key, {}});
    }

    for(auto& [key, inst] : rule_instances)
    {
        auto& p = system_graph.findNode(inst.anchor)->second.getData().position;
        int ic, jc;
        grid.reaction_grid.locatePoint(p[0], p[1], ic, jc);
        auto cardinal = grid.reaction_grid.cardinalCellIndex(ic, jc);
        // here would be the anchor list reduction, plus whatever other code needs to be added
        //anchor_list.insert({match.first.anchor, cardinal});
        auto max_cell = grid.cell_label[cardinal];
        rule_map[max_cell].push_back(key);
        //std::cout << "rule is mapped to cell " << max_cell << "\n";
    }

    for(auto& [k, v] : rule_map)
    {
        auto d = grid.graph.findNode(k)->second.getData().type;
        std::cout << d << "D cell " << k << ": " << v.size() << "\n";
    }

}

//TEST_CASE("Phi Decomposition", "[Phi Test]")
//{
//    //TODO: fix epsilon here, causing errors
//    DGGML::Grammar<GraphType> gamma;
//    define_model(gamma);
//    DGGML::AnalyzedGrammar<GraphType> gamma_analysis(gamma);
//    DGGML::ExpandedComplex2D<> grid(2, 2, 4.0, 4.0, false, 1.0);
//    DGGML::VtkFileWriter<typename DGGML::ExpandedComplex2D<>::types::graph_type> writer;
//    writer.save(grid.getGraph(), "test_cell_complex");
//    DGGML::GridFileWriter grid_writer;
//    grid_writer.save({grid.reaction_grid, grid.dim_label}, "test_grid");
//
//    GraphType system_graph;
//    initialize_system(system_graph);
//    REQUIRE(system_graph.numNodes() == 42);
//    REQUIRE(YAGL::connected_components(system_graph) == 14);
//    DGGML::VtkFileWriter<GraphType> vtk_writer;
//    vtk_writer.save(system_graph, "phi_graph");
//
//}