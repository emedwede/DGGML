#include <iostream>

#include "catch.hpp"

#include "YAGL_Graph.hpp"
#include "YAGL_Algorithms.hpp"

#include "RuleSystem.hpp"
#include "Grammar.h"
#include "AnalyzedGrammar.hpp"
#include "ExpandedComplex2D.hpp"
#include "VtkWriter.hpp"
#include "MathUtils.hpp"

#include "../examples/CMA/PlantModel/PlantTypes.hpp"


using GraphType = Plant::graph_type;
using KeyType = typename GraphType::key_type;
using NodeType = typename GraphType::node_type;

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
                                  [](auto& lhs, auto& rhs, auto& m) {});
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
                                  [](auto& lhs, auto& rhs, auto& m) {});
    gamma.addRule(r2);
}

//initializes a simple MT system for testing
void initialize_system(GraphType& system_graph)
{
    //Add two microtubule segments parallel to each other
    system_graph.addNode({100,{Plant::Negative{0.0, 0.0, 0.0,
                                               0.0, 0.0, 0.0},
                               0.1, 0.1, 0.0}});
    system_graph.addNode({101,{Plant::Intermediate{0.0, 0.0, 0.0,
                                                   0.0, 0.0, 0.0},
                               0.5, 0.1, 0.0}});
    system_graph.addNode({102,{Plant::Positive{0.0, 0.0, 0.0,
                                               0.0, 0.0, 0.0},
                               0.9, 0.1, 0.0}});
    system_graph.addEdge(100, 101);
    system_graph.addEdge(101, 102);

    system_graph.addNode({103,{Plant::Negative{0.0, 0.0, 0.0,
                                               0.0, 0.0, 0.0},
                               0.1, 1.1, 0.0}});
    system_graph.addNode({104,{Plant::Intermediate{0.0, 0.0, 0.0,
                                                   0.0, 0.0, 0.0},
                               0.5, 1.1, 0.0}});
    system_graph.addNode({105,{Plant::Positive{0.0, 0.0, 0.0,
                                               0.0, 0.0, 0.0},
                               0.9, 1.1, 0.0}});
    system_graph.addEdge(103, 104);
    system_graph.addEdge(104, 105);


    system_graph.addNode({106,{Plant::Negative{0.0, 0.0, 0.0,
                                               0.0, 0.0, 0.0},
                               2.1, 2.1, 0.0}});
    system_graph.addNode({107,{Plant::Intermediate{0.0, 0.0, 0.0,
                                                   0.0, 0.0, 0.0},
                               2.5, 2.5, 0.0}});
    system_graph.addNode({108,{Plant::Positive{0.0, 0.0, 0.0,
                                               0.0, 0.0, 0.0},
                               2.9, 2.9, 0.0}});
    system_graph.addEdge(106, 107);
    system_graph.addEdge(107, 108);
}

void save_state(GraphType& system_graph, DGGML::ExpandedComplex2D<>& grid)
{
    DGGML::VtkFileWriter<typename DGGML::ExpandedComplex2D<>::types::graph_type> writer;
    writer.save(grid.getGraph(), "incremental_cell_complex");
    DGGML::GridFileWriter grid_writer;
    grid_writer.save({grid.reaction_grid, grid.dim_label}, "incremental_grid");
    DGGML::VtkFileWriter<GraphType> vtk_writer;
    vtk_writer.save(system_graph, "incremental_graph");
}

//TODO: start with simple single component rules and then do multi-component tests later
/* TODO: ideas for unifying the indexed rule collection or whatever it ends up being called.
     *  Minimum requirements: rules matches -> connected component matches -> ? maybe nodes.
     *  Bonus minimum: geocell -> rules matches
     *  Consideration: what is the complexity of just doing a forward search for rules containing a node i.e.
     *  is node in any of the connected components of any of the rule matches assign to a geocells subcell?
     */
TEST_CASE("Incremental Update Test", "[incremental-update-test]")
{
    std::cout << "\n" << "Running the incremental update test" << "\n";
    DGGML::Grammar<GraphType> gamma;
    define_model(gamma);
    DGGML::AnalyzedGrammar<GraphType> gamma_analysis(gamma);

    GraphType system_graph;
    initialize_system(system_graph);

    DGGML::ExpandedComplex2D<> grid(1, 1, 3.0, 3.0, false, 1.0);
    save_state(system_graph, grid);

    //find matches
    auto c1 = gamma_analysis.unique_components[0];
    std::cout << c1.numNodes() << "\n";
    auto end_matches = YAGL::subgraph_isomorphism2(c1, system_graph);
    REQUIRE(end_matches.size() == 3);
    auto c2 = gamma_analysis.unique_components[1];
    std::cout << c2.numNodes() << "\n";
    auto mt_matches = YAGL::subgraph_isomorphism2(c2, system_graph);
    REQUIRE(mt_matches.size() == 3);

    //create the ordering
    std::vector<std::size_t> end_ordering;
    for(auto& item : end_matches[0])
        end_ordering.push_back(item.first);
    std::vector<std::size_t> mt_ordering;
    for(auto& item : mt_matches[0])
        mt_ordering.push_back(item.first);

    //build the set of unique component matches
    // (TODO: should change name of rule_system to be something like components matches)
    //DGGML::RuleSystem<KeyType> rule_system;
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

    for(auto& [key, inst] : component_match_set)
    {
        if(inst.type == 0)
        {
            std::cout << "End Match, key: " << key << "\n";
            for(auto i = 0; i < end_ordering.size(); i++)
            {
                std::cout << "{ " << end_ordering[i] << " -> " << inst.match[i] << " }\n";
            }
            std::cout << "\n";
        }
        if(inst.type == 1)
        {
            std::cout << "MT Match, key: " << key << "\n";
            for(auto i = 0; i < mt_ordering.size(); i++)
            {
                std::cout << "{ " << mt_ordering[i] << " -> " << inst.match[i] << " }\n";
            }
            std::cout << "\n";
        }
    }

    REQUIRE(component_match_set.size() == 6);

    //building the rule instances
    int kk = 0;
    std::map<std::size_t, DGGML::RuleInstType<KeyType>> rule_instances;
    for(int i = 0; i < component_match_set.size(); i++)
    {
        if(component_match_set[i].type == 0) {
            rule_instances[kk];
            rule_instances[kk].name = "with_rule";
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
            if(comp1.type == 1 && comp2.type == 1 && d <= 1.0)
            {
                rule_instances[kk];
                rule_instances[kk].name = "incremental_rule";
                rule_instances[kk].category = "stochastic";
                rule_instances[kk].components.push_back(i);
                rule_instances[kk].components.push_back(j);
                rule_instances[kk].anchor = comp1.anchor;
                kk++;
            }
        }
    }

    REQUIRE(rule_instances.size() == 4);

    for(auto& [key, rinst] : rule_instances)
    {
        std::cout << "Component keys for rule " << key << ": ";
        for(auto& comp : rinst.components)
        {
            std::cout << "{ ";
            for(auto& key2 : component_match_set[comp].match)
                std::cout << key2 << " ";
            std::cout << "} ";
        }
        std::cout << "\n";
    }
    // hierarchy: (geocells) -> rule instances -> components instances -> nodes
    //Test: remove matches with node 101
    //slow, but easy case to code - assume we know nothing and must search everywhere
    std::set<std::size_t> node_rule_invalidations;
    std::set<std::size_t> node_component_invalidations;
    for(auto& [k1, rinst] : rule_instances)
    {
        bool found = false;
        for(auto& k2: rinst.components)
        {
            for (auto& k3 : component_match_set[k2].match)
            {
                if(k3 == 101) // found so must be marked for invalidation
                {
                    node_rule_invalidations.insert(k1);
                    node_component_invalidations.insert(k2);
                    found = true;
                    break; //here we can break after first of, since rule components do not overlap
                }
            }
            if(found) break;
        }
    }
    std::cout << "Rules to invalidate: ";
    for(auto id : node_rule_invalidations) std::cout << id << " ";
    std::cout << "\n";
    std::cout << "Components to invalidate: ";
    for(auto id : node_component_invalidations) std::cout << id << " ";
    std::cout << "\n\n";

    //can definitely combine these sets with the node invalidation ones
    std::set<std::size_t> edge_rule_invalidations;
    std::set<std::size_t> edge_component_invalidations;
    //Test: search for edge (101, 102) to remove
    for(auto& [k1, rinst] : rule_instances)
    {
        bool found = false;
        for(auto& k2 : rinst.components)
        {
            for(auto i = 0; i < component_match_set[k2].match.size(); i++)
            {
                auto k3 = component_match_set[k2].match[i];
                for(auto j = 0; j < component_match_set[k2].match.size(); j++)
                {
                    auto k4 = component_match_set[k2].match[j];
                    auto& c = gamma_analysis.unique_components[component_match_set[k2].type];
                    //match contains the nodes
                    std::vector<std::size_t> ordering;
                    if(component_match_set[k2].type == 0) ordering = end_ordering;
                    else ordering = mt_ordering;
                    if(k3 == 101 && k4 == 102 && c.adjacent(ordering[i], ordering[j]))
                    {
                        std::cout << "{ " << end_ordering[i] << " -> " << k3 << " }\n";
                        std::cout << "{ " << end_ordering[j] << " -> " << k4 << " }\n";
                        std::cout << "rule " << k1 << ", component " << k2 << " contains the edge\n";
                        edge_rule_invalidations.insert(k1);
                        edge_component_invalidations.insert(k2);
                        found = true;
                        break;
                    }
                }
                if(found) break;
            }
            if(found) break;
        }
    }
    std::cout << "Rules to invalidate: ";
    for(auto id : edge_rule_invalidations) std::cout << id << " ";
    std::cout << "\n";
    std::cout << "Components to invalidate: ";
    for(auto id : edge_component_invalidations) std::cout << id << " ";
    std::cout << "\n\n";


    //Say we selected a rule to fire and we had to modify the system_graph by
    //adding node 200
    //removing edge 104, 105
    //adding edge 104, 106 and edge 105, 106
    //How does that change the rule_instance set?
    system_graph.addNode({200,{Plant::Intermediate{0.0, 0.0, 0.0,
                                                   0.0, 0.0, 0.0},
                               0.5, 0.0, 0.0}});
    system_graph.removeEdge(104, 105);
    system_graph.addEdge(104, 200);
    system_graph.addEdge(105, 200);
    std::vector<std::size_t> added_nodes = {200};
    std::vector<std::pair<std::size_t, std::size_t>> added_edges = {{104, 200}, {105, 200}};
    //REQUIRE(YAGL::connected_components(system_graph) == 2);

    std::cout << system_graph << "\n";

    //goal, take the candidate nodes, do a search the depth of the height of the tallest rooted spanning tree
    //for each candidate and create a set of nodes used to induce a graph that we will search for new components
    std::set<std::size_t> candidate_nodes = {104, 105, 200};
    std::set<std::size_t> inducers;
    for(auto n : candidate_nodes)
    {
        auto res = YAGL::recursive_dfs(system_graph, n, 2);
        for(auto& item : res)
            inducers.insert(item);
    }

    auto candidate_graph = YAGL::induced_subgraph(system_graph, inducers);
    std::cout << candidate_graph << "\n";

    //search the candidate graph for components
    auto candidate_matches = YAGL::subgraph_isomorphism2(c1, candidate_graph);
    std::cout << candidate_matches.size() << "\n";

    //now we must refine the candidate matches and reject the ones that do no contain a newly added node or edge
    decltype(candidate_matches) accepted_matches;
    for(auto& m : candidate_matches)
    {
        bool accepted = false;
        auto& c = gamma_analysis.unique_components[0];
        for(auto& [v1, v2] : m)
        {
            for(auto& n : added_nodes) {
                if (v2 == n) {
                    accepted = true;
                    accepted_matches.push_back(m);
                    break;
                }
            }
            if(accepted) break;
        }
        if(accepted) continue;
        for(auto& [v1, v2] : m)
        {
            for(auto& [u1, u2] : m)
            {
                for(auto& [n1, n2] : added_edges) {
                    if (v2 == n1 && u2 == n2 && c.adjacent(v1, u1)) {
                        accepted = true;
                        accepted_matches.push_back(m);
                        break;
                    }
                    if(accepted) break;
                }
            }
            if(accepted) break;
        }
        if(accepted) continue;
    }

    std::cout << "Accepted: \n\n";
    for(auto& m : accepted_matches)
    {
        std::cout << "Match: \n";
        for(auto& [v1, v2] : m)
        {
            std::cout << "{ " << v1 << " -> " << v2 << " }\n";
        }
        std::cout << "\n";
    }
    //given a series of new nodes/edges being added, how do we know where to search for new components?
}

//void print_matches(std::vector<std::vector<key_type>>& matches)
//{
//    std::cout << "Found " << matches.size() << " Matches: \n";
//    for(auto& m : matches)
//    {
//        std::cout << "\t{ ";
//        for(auto& v : m)
//        {
//            std::cout << v << " ";
//        } std::cout << "}\n";
//    }
//}
//
//void print_matches(std::unordered_map<key_type, std::vector<key_type>>& matches)
//{
//    std::cout << "Found " << matches.size() << " Matches: \n";
//    for(auto& [key, m] : matches)
//    {
//        std::cout << "\t{ ";
//        for(auto& v : m)
//        {
//            std::cout << v << " ";
//        } std::cout << "}\n";
//    }
//}
//
//template <typename GraphType>
//void microtubule_scatter(GraphType& graph, std::size_t num_mt)
//{
//    double epsilon_min = 0.5;
//    double epsilon_max = 1.0;
//    std::random_device random_device;
//    std::mt19937 random_engine(random_device());
//
//    std::uniform_real_distribution<double> distribution_global_x(2, 10);
//
//    std::uniform_real_distribution<double> distribution_global_y(2, 10);
//
//    std::uniform_real_distribution<double> distribution_local(epsilon_min, epsilon_max);
//
//    std::uniform_real_distribution<double> distribution_angle(0.0, 2.0*3.14);
//    using node_type = typename GraphType::node_type;
//
//    std::size_t segments = 3;
//    for(auto i = 0; i < num_mt ; i++)
//    {
//        double x_c, y_c;
//        double z_c = 0;
//        x_c = distribution_global_x(random_engine);
//        y_c = distribution_global_y(random_engine);
//        auto theta = distribution_angle(random_engine);
//        auto seg_len = distribution_local(random_engine);
//
//        auto x_s = 0.0;
//        auto y_s = seg_len;
//        auto x_r_t = x_s*cos(theta) + y_s*sin(theta);
//        auto y_r_t = -x_s*sin(theta) +y_s*cos(theta);
//        auto x_r = x_c + x_r_t;
//        auto y_r = y_c + y_r_t;
//        auto z_r = 0.0;
//
//        auto x_l = x_c - (x_r - x_c);
//        auto y_l = y_c - (y_r - y_c);
//        auto z_l = 0.0;
//
//        //compute dist and unit vector
//        double p1[3] = {x_r - x_c, y_r - y_c, 0.0};
//        double p2[3] = {0.0, 0.0, 0.0};
//        double p3[3] = {x_c - x_l, y_c - y_l, 0.0};
//        double u1[3], u2[3];
//
//        DGGML::set_unit_vector(p1, p2, u1);
//        DGGML::set_unit_vector(p3, p2, u2);
//
//        node_type node_l(i*segments,
//                {{x_l, y_l, z_l},
//                 {0.0, 0.0, 0.0},
//                 DGGML::Plant::mt_type::negative,
//                 {-1, -1, -1},
//                 {u2[0], u2[1], u2[2]}});
//
//        node_type node_c(i*segments+1,
//                {{x_c, y_c, z_c},
//                 {0.0, 0.0, 0.0},
//                 DGGML::Plant::mt_type::intermediate,
//                 {-1, -1, -1},
//                 {u1[0], u1[1], u1[2]}});
//
//        node_type node_r(i*segments+2,
//                {{x_r, y_r, z_r},
//                 {0.0, 0.0, 0.0},
//                 DGGML::Plant::mt_type::positive,
//                 {-1, -1, -1},
//                 {u1[0], u1[1], u1[2]}});
//
//        graph.addNode(node_l);
//        graph.addNode(node_c);
//        graph.addNode(node_r);
//
//        graph.addEdge(node_l, node_c);
//        graph.addEdge(node_r, node_c);
//    }
//}
//
//
////Simple first attempt a polymerizing
//template <typename GraphType>
//std::pair<std::set<key_type>, std::set<key_type>> test_rewrite(GraphType& graph, std::vector<key_type>& match)
//{
//    //if(match.size() != 2) return;
//    auto i = match[0]; auto j = match[1];
//
//    //TODO: need a unique key generator
//    typename GraphType::key_type key = graph.numNodes()+1;
//
//    while(graph.findNode(key) != graph.node_list_end()) key++; //TODO: fix, very greedy
//    double x3[3];
//
//    auto& x1 = graph.findNode(i)->second.getData().position;
//    auto& x2 = graph.findNode(j)->second.getData().position;
//    auto& u1 = graph.findNode(i)->second.getData().unit_vec;
//    auto gamma = 0.75;
//    for(auto iter = 0; iter < 3; iter++)
//    {
//        x3[iter] = x2[iter] - ((x2[iter]-x1[iter]) * gamma);
//    }
//    std::random_device random_device; std::mt19937 random_engine(random_device());
//    std::uniform_real_distribution<double> distribution_angle(-3.14/8.0, 3.14/8.0);
//    double theta = distribution_angle(random_engine);
//    auto u10_rot = u1[0]*cos(theta) + u1[1]*sin(theta);
//    auto u11_rot = -u1[0]*sin(theta) +u1[1]*cos(theta);
//    u1[0] = u10_rot; u1[1] = u11_rot;
//
//    graph.addNode({key,
//            {{x3[0], x3[1], x3[2]},
//             {0, 0, 0},
//             DGGML::Plant::mt_type::intermediate,
//             {0, 0, 0},
//             {u1[0], u1[1], u1[2]}}});
//
//    graph.removeEdge(i, j);
//    graph.addEdge(i, key);
//    graph.addEdge(j, key);
//
//    return std::pair<std::set<key_type>, std::set<key_type>>{{i, j}, {i, j, key}}; //matches to invalidate and induce
//}
//
//auto random_rewrite(std::size_t n)
//{
//    //randomly select a rewrite to occur
//    std::random_device random_device;
//    std::mt19937 random_engine(random_device());
//
//    std::uniform_int_distribution<std::size_t> rand_rewrite(0, n-1);
//
//    auto selected = rand_rewrite(random_engine);
//    std::cout << "Rule selected to rewrite: " << selected << "\n";
//    return selected;
//}
//
//TEST_CASE("Incremental Update Test", "[update-test]")
//{
//    graph_type graph;
//
//    std::size_t n = 10; // num mt
//    microtubule_scatter(graph, n);
//
//    std::vector<DGGML::Plant::mt_key_type> bucket;
//
//    for(const auto& [key, value] : graph.getNodeSetRef())
//        bucket.push_back(key);
//
//    // compute all the matches
//    auto matches = DGGML::Plant::microtubule_growing_end_matcher(graph, bucket);
//    //print_matches(matches);
//    REQUIRE(matches.size() == n);
//
//    //randomly select a rewrite to occur
//    auto selected = random_rewrite(n);
//
//    std::cout << "Rule selected to rewrite: " << selected << "\n";
//
//    std::unordered_map<std::size_t, std::vector<key_type>> gi_map;
//    for(std::size_t i = 0; i < matches.size(); i++)
//    {
//        gi_map.insert({i, matches[i]});
//    }
//    print_matches(gi_map);
//    //each node needs to know what match it belongs to
//    std::unordered_map<key_type, std::vector<std::size_t>> gi_inverse;
//    for(const auto& key : bucket)
//       gi_inverse.insert({key, {}});
//
//    for(const auto& [key, match] : gi_map)
//    {
//        for(const auto& item : match)
//        {
//            auto search = gi_inverse.find(item);
//            if(search != gi_inverse.end())
//            {
//                //shoud probably insert key into set not vector
//                search->second.push_back(key);
//            }
//        }
//    }
//
//    REQUIRE(gi_map.size() == n);
//
//    auto lambda_print = [](std::unordered_map<key_type, std::vector<std::size_t>>& gi_inverse)
//    {
//        for(const auto& [key, value] : gi_inverse)
//        {
//            std::cout << "Node " << key << ": {";
//            for(const auto& node : value)
//                std::cout << " " << node << " ";
//            std::cout << "}\n";
//        }
//    };
//    lambda_print(gi_inverse);
//    // rewrite the graph and gather the nodes the rule invalidates
//    auto [invalidations, inducers] = test_rewrite(graph, matches[selected]);
//
//    //first test the graph was rewritten
//    REQUIRE(YAGL::connected_components(graph) == n);
//    // only a single new unique node is added
//    REQUIRE(graph.numNodes() == n*3+1);
//
//    REQUIRE(invalidations.size() == 2);
//
//    //invalidate old matches
//    for(const auto& i : invalidations)
//    {
//        const auto& rules = gi_inverse.find(i);
//        if(rules != gi_inverse.end())
//        {
//            //TODO: need to remove matches from inverse map in a less destructive way
//            for(const auto& j : rules->second)
//            {
//                gi_map.erase(j);
//            }
//            gi_inverse.erase(i);
//        }
//    }
//
//    REQUIRE(gi_map.size() == n - 1);
//    print_matches(gi_map);
//
//    auto temp = inducers;
//    // get all nodes one hop away and add those to inducers
//    for(const auto& i : temp)
//    {
//        auto res = YAGL::recursive_dfs(graph, i, 2);
//        for(auto& j : res)
//            inducers.insert(j);
//    }
//    REQUIRE(inducers.size() == 4);
//
//    //induce a subgraph to search for new matches
//    auto subgraph = YAGL::induced_subgraph(graph, inducers);
//
//    REQUIRE(subgraph.numNodes() == 4);
//
//    std::vector<DGGML::Plant::mt_key_type> sub_bucket;
//
//    for(const auto& [key, value] : subgraph.getNodeSetRef())
//        sub_bucket.push_back(key);
//
//    // compute all the matches
//    auto incremental_matches = DGGML::Plant::microtubule_growing_end_matcher(graph, sub_bucket);
//    //print_matches(matches);
//    REQUIRE(incremental_matches.size() == 1);
//    for(auto& key : sub_bucket)
//    {
//        if(gi_inverse.find(key) == gi_inverse.end())
//            gi_inverse.insert({key, {}});
//    }
//    int iii = 1;
//    for(const auto& match : incremental_matches)
//    {
//        auto key = n + iii;
//        gi_map.insert({key, match}); iii++;
//        for(const auto& item : match)
//        {
//            auto search = gi_inverse.find(item);
//            if(search != gi_inverse.end())
//            {
//                search->second.push_back(key);
//            }
//        }
//    }
//    REQUIRE(gi_map.size() == n);
//    print_matches(gi_map);
//    lambda_print(gi_inverse);
//}
//
//TEST_CASE("RuleSystem Test", "[rule-system-test]")
//{
//    graph_type graph;
//
//    std::size_t n = 10; // num mt
//    microtubule_scatter(graph, n);
//
//    std::vector<DGGML::Plant::mt_key_type> bucket;
//
//    for(const auto& [key, value] : graph.getNodeSetRef())
//        bucket.push_back(key);
//
//    // compute all the matches
//    auto matches = DGGML::Plant::microtubule_growing_end_matcher(graph, bucket);
//    //print_matches(matches);
//    REQUIRE(matches.size() == n);
//
//    DGGML::RuleSystem<std::size_t> rule_system;
//    REQUIRE(rule_system.size() == 0);
//
//    for(auto& item : matches)
//        rule_system.push_back({std::move(item), DGGML::Rule::G});
//    REQUIRE(rule_system.size() == matches.size());
//
//    matches = DGGML::Plant::microtubule_retraction_end_matcher(graph, bucket);
//    REQUIRE(matches.size() == n);
//
//    for(auto& item : matches)
//        rule_system.push_back({std::move(item), DGGML::Rule::R});
//
//    REQUIRE(rule_system.size() == 2*n);
//
//    rule_system.print_index();
//
//    auto grow = std::count_if(rule_system.begin(), rule_system.end(),
//            [](const auto& item) {
//        if(item.first.type == DGGML::Rule::R)
//            return true;
//        else return false;
//    });
//
//    REQUIRE(grow == n);
//
//    REQUIRE(rule_system.inverse_index.size() == 30);
//
//    auto selected = random_rewrite(n);
//
//    auto lambda_print = [](typename DGGML::RuleSystem<std::size_t>::inverse_type& gi_inverse)
//    {
//        for(const auto& [key, value] : gi_inverse)
//        {
//            std::cout << "Node " << key << ": {";
//            for(const auto& node : value)
//                std::cout << " " << node << " ";
//            std::cout << "}\n";
//        }
//    };
//    lambda_print(rule_system.inverse_index);
//
//    REQUIRE(rule_system[0].type == DGGML::Rule::G);
//
//    // rewrite the graph and gather the nodes the rule invalidates
//    auto [invalidations, inducers] = test_rewrite(graph, rule_system[selected].match);
//
//    //first test the graph was rewritten
//    REQUIRE(YAGL::connected_components(graph) == n);
//    // only a single new unique node is added
//    REQUIRE(graph.numNodes() == n*3+1);
//
//    REQUIRE(invalidations.size() == 2);
//
//    rule_system.print_index();
//    //invalidate old matches
//    for(const auto& i : invalidations)
//        rule_system.invalidate(i);
//
//    REQUIRE(rule_system.size() == 2*n - 2);
//
//    rule_system.print_index();
//
//
//    auto temp = inducers;// get all nodes one hop away and add those to inducers
//    for(const auto& i : temp)
//    {
//        auto res = YAGL::recursive_dfs(graph, i, 2);
//        for(auto& j : res)
//            inducers.insert(j);
//    }
//    REQUIRE(inducers.size() == 4);
//
//    //induce a subgraph to search for new matches
//    auto subgraph = YAGL::induced_subgraph(graph, inducers);
//
//    REQUIRE(subgraph.numNodes() == 4);
//
//    std::vector<DGGML::Plant::mt_key_type> sub_bucket;
//
//    for(const auto& [key, value] : subgraph.getNodeSetRef())
//        sub_bucket.push_back(key);
//
//    // compute all the matches
//    auto incremental_matches = DGGML::Plant::microtubule_growing_end_matcher(graph, sub_bucket);
//    REQUIRE(incremental_matches.size() == 1);
//
//    for(auto& item : incremental_matches)
//        rule_system.push_back({std::move(item), DGGML::Rule::G});
//
//    incremental_matches = DGGML::Plant::microtubule_retraction_end_matcher(graph, sub_bucket);
//    REQUIRE(incremental_matches.size() == 1);
//
//    for(auto& item : incremental_matches)
//        rule_system.push_back({std::move(item), DGGML::Rule::R});
//
//    REQUIRE(rule_system.size() == 2*n);
//
//    lambda_print(rule_system.inverse_index);
//
//    rule_system.print_index();
//
//    auto num_grow = rule_system.count(DGGML::Rule::G);
//    REQUIRE(num_grow == n);
//
//    auto num_retract = rule_system.count(DGGML::Rule::R);
//    REQUIRE(num_retract == n);
//}
//
//#include "CartesianComplex2D.hpp"
//
//TEST_CASE("CellList Test", "[cell-list-test]")
//{
//    graph_type graph;
//
//    std::size_t n = 10; // num mt
//    microtubule_scatter(graph, n);
//
//    DGGML::CartesianComplex2D cplex(2, 2, 6.0, 6.0);
//    using cplex_key_t = typename DGGML::CartesianComplex2D<>::graph_type::key_type;
//    using plant_key_t = DGGML::Plant::mt_key_type;
//
//    std::vector<cplex_key_t> bucket2d;
//    for(auto& item : cplex.graph.getNodeSetRef())
//    {
//        if(item.second.getData().type == 0)
//            bucket2d.push_back(item.first);
//    }
//    for(auto& item : bucket2d) std::cout << item << " ";
//    std::cout << "\n";
//    REQUIRE(bucket2d.size() == cplex.get2dTotalCellCount());
//
//    std::unordered_map<plant_key_t, cplex_key_t> cell_list;
//
//    for(auto& [key, value] : graph.getNodeSetRef())
//    {
//        auto& node_data = value.getData();
//        double xp = node_data.position[0];
//        double yp = node_data.position[1];
//        int ic, jc;
//
//        cplex.coarse_grid.locatePoint(xp, yp, ic, jc);
//        cplex.coarse_cell_to_fine_lattice(ic, jc);
//        auto cardinal = cplex.fine_grid.cardinalLatticeIndex(ic, jc);
//        std::cout << cardinal << " ";
//        cell_list.insert({key, cardinal});
//    }
//    std::cout << "\n";
//    REQUIRE(cell_list.size() == graph.numNodes());
//
//    for(auto& c : bucket2d)
//    {
//        auto total = std::count_if(cell_list.begin(), cell_list.end(), [&c](auto& item) {
//            if(item.second == c)
//                return true;
//            else return false;
//        });
//        std::cout << c << ": " << total << "\n";
//    }
//
//    std::vector<DGGML::Plant::mt_key_type> bucket;
//
//    for(const auto& [key, value] : graph.getNodeSetRef())
//        bucket.push_back(key);
//
//    auto matches = DGGML::Plant::microtubule_growing_end_matcher(graph, bucket);
//
//    DGGML::RuleSystem<std::size_t> rule_system;
//
//    using rule_key_t = std::size_t;
//    std::map<cplex_key_t, std::vector<rule_key_t>> rule_map;
//    for(const auto& item : bucket2d)
//        rule_map.insert({item, {}});
//
//    for(auto& item : matches)
//        rule_system.push_back({std::move(item), DGGML::Rule::G});
//
//    std::cout << "Num matches: " << rule_system.size() << "\n";
//    //map the rules to particular cell
//    for(const auto& match : rule_system)
//    {
//        auto& instance = match.first.match;
//        for(auto& l : instance)
//            std::cout << cell_list.find(l)->second << " ";
//        std::cout << "\n";
//        auto it =
//        std::adjacent_find(instance.begin(), instance.end(), [&cell_list](auto& lhs, auto& rhs){
//            if(cell_list.find(rhs)->second != cell_list.find(lhs)->second)
//                return true;
//            else
//                return false;
//        });
//
//        if(it == instance.end())
//        {
//            std::cout << "maps to highest dim cell\n";
//            auto cell_id = cell_list.find(instance.front())->second;
//            rule_map[cell_id].push_back(match.second);
//        }
//        else
//            std::cout << "maps to lower dim cell\n";
//    }
//    for(const auto& c : rule_map)
//    {
//        std::cout << "{ Cell ID: " << c.first << ", Rule IDs: { ";
//        for(auto& i : c.second)
//            std::cout << i << " ";
//        std::cout << "} }\n";
//    }
//}
