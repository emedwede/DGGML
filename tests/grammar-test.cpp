#include <iostream>
#include <memory>
#include <algorithm>

#include "catch.hpp"

#include "YAGL_Graph.hpp"
#include "YAGL_Algorithms.hpp"
//#include "Utlities/MemoryManager.hpp"
#include "../examples/CMA/PlantModel/PlantTypes.hpp"
#include "Grammar.h"
#include "AnalyzedGrammar.hpp"
#include "ComponentMatchMap.hpp"
#include "RuleMatchMap.hpp"
#include "IncrementalUpdate.hpp"
#include "MathUtils.hpp"

#include <unordered_map> 

using GraphType = Plant::graph_type;

struct KeyGenerator
{
    std::size_t current_key;

    KeyGenerator(std::size_t current_key) : current_key(current_key) {}

    std::size_t getKey()
    {
        return current_key++;
    }
};

//this will need extensive error checking since it's user defined?
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
                                  [](auto& lhs, auto& m) { return 2.0; },
                                  [](auto& lhs, auto& rhs, auto& m1, auto& m2) {
        auto d = DGGML::calculate_distance(lhs[m1[1]].position, lhs[m1[2]].position);
        //std::cout << d << "\n";
        //auto& data1 = std::get<Plant::Intermediate>(lhs[m[1]].data);
        std::cout << "doing some update calculations\n";
        rhs[m2[3]].position[0] = (lhs[m1[2]].position[0] - lhs[m1[1]].position[0])/2;
        rhs[m2[3]].position[1] = (lhs[m1[2]].position[1] - lhs[m1[1]].position[1])/2;
        rhs[m2[3]].position[2] = (lhs[m1[2]].position[2] - lhs[m1[1]].position[2])/2;
    });
    gamma.addRule(r1);

    GraphType g3;
    g3.addNode({3, {Plant::Intermediate{}}});
    g3.addNode({4, {Plant::Positive{}}});
    g3.addEdge(3, 4);

    DGGML::SolvingRule<GraphType> r2("solving_growth", g3, g3,[](auto& lhs) {return 2.0;});
    gamma.addRule(r2);

    GraphType g4;
    g4.addNode({1, {Plant::Negative{}}});
    g4.addNode({2, {Plant::Intermediate{}}});
    g4.addNode({3, {Plant::Positive{}}});
    g4.addEdge(1, 2);
    g4.addEdge(2, 3);

    g4.addNode({4, {Plant::Negative{}}});
    g4.addNode({5, {Plant::Intermediate{}}});
    g4.addNode({6, {Plant::Positive{}}});
    g4.addEdge(4, 5);
    g4.addEdge(5, 6);

    GraphType g5;

    DGGML::WithRule<GraphType> r3("with_interaction", g4, g5,
                                  [](auto& lhs, auto& m) { return 0.0; },
                                  [](auto& lhs, auto& rhs, auto& m1, auto& m2) {});
    gamma.addRule(r3);

}

//initializes a simple MT system for testing
void initialize_system(GraphType& system_graph)
{
    //Add two microtubule segments parallel to each other
    system_graph.addNode({100,{Plant::Negative{0.0, 0.0, 0.0,
                                           0.0, 0.0, 0.0},
                           0.0, 0.0, 0.0}});
    system_graph.addNode({101,{Plant::Intermediate{0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0},
                             0.5, 0.0, 0.0}});
    system_graph.addNode({102,{Plant::Positive{0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0},
                             1.0, 0.0, 0.0}});
    system_graph.addEdge(100, 101);
    system_graph.addEdge(101, 102);

    system_graph.addNode({103,{Plant::Negative{0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0},
                             0.0, 2.0, 0.0}});
    system_graph.addNode({104,{Plant::Intermediate{0.0, 0.0, 0.0,
                                                 0.0, 0.0, 0.0},
                             0.5, 2.0, 0.0}});
    system_graph.addNode({105,{Plant::Positive{0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0},
                             1.0, 2.0, 0.0}});
    system_graph.addEdge(103, 104);
    system_graph.addEdge(104, 105);
}

TEST_CASE("Basic Grammar Test", "[basic-grammar-test]")
{
    std::cout << "\n" << "Running basic grammar test" << "\n";
    DGGML::Grammar<GraphType> gamma;
    define_model(gamma);

    REQUIRE(gamma.stochastic_rules.at("with_growth").lhs_graph.numNodes() == 2);
    REQUIRE(gamma.stochastic_rules.at("with_growth").lhs_graph.numEdges() == 1);

    REQUIRE(gamma.stochastic_rules.at("with_growth").rhs_graph.numNodes() == 3);
    REQUIRE(gamma.stochastic_rules.at("with_growth").rhs_graph.numEdges() == 2);

    gamma.print();

}


TEST_CASE("Basic Grammar Analysis Test", "[basic-grammar-analysis-test]")
{
    std::cout << "\n" << "Running basic grammar analysis test" << "\n";
    DGGML::Grammar<GraphType> gamma;
    define_model(gamma);

    DGGML::AnalyzedGrammar<GraphType> gamma_analysis(gamma);

    REQUIRE(gamma_analysis.with_rules.at("with_growth").lhs_graph.numNodes() == 2);
    REQUIRE(gamma_analysis.with_rules.at("with_growth").lhs_graph.numEdges() == 1);

    REQUIRE(gamma_analysis.with_rules.at("with_growth").rhs_graph.numNodes() == 3);
    REQUIRE(gamma_analysis.with_rules.at("with_growth").rhs_graph.numEdges() == 2);

    REQUIRE(gamma_analysis.unique_components.size() == 2);
    REQUIRE(gamma_analysis.with_rewrites.at("with_growth").node_set_create.size() == 1);
    REQUIRE(gamma_analysis.with_rewrites.at("with_growth").node_set_destroy.empty());
    REQUIRE(gamma_analysis.with_rewrites.at("with_growth").edge_set_create.size() == 2);
    REQUIRE(gamma_analysis.with_rewrites.at("with_growth").edge_set_destroy.size() == 1);
    gamma.print();

}

//testable version
void perform_rewrite(DGGML::RuleMatch<std::size_t>& inst,
                      std::map<std::size_t, DGGML::ComponentMatch<std::size_t>>& component_match_set,
                      KeyGenerator& gen,
                     DGGML::AnalyzedGrammar<GraphType>& gamma_analysis, GraphType& system_graph)
{
    std::cout << "\n" << "Beginning rewrite code ..." << "\n";

    std::string rname = inst.name;

    //construct a vertex map for the lhs to the rule instance
    std::map<std::size_t, std::size_t> lhs_vertex_map;
    std::vector<std::size_t> left, mid, right;
    for(auto& m : gamma_analysis.ccuv_mappings[rname]) {
        for (auto &[k, v]: m) {
            left.push_back(k);
            mid.push_back(v);
            std::cout << k << " " << v << "\n";
        }
    }

    for(auto& c : inst.components)
    {
        for(auto& k : component_match_set[c].match)
        {
            right.push_back(k); std::cout << k << "\n";
        }
    }
    REQUIRE((right.size() == mid.size() && mid.size() == left.size()));

    //print out the mapping info, so we know how a lhs numbering maps to an instance numbering
    std::cout << "\nMappings: { LHS Key -> Minimal Component Key -> Rule ComponentMatch Key }\n";
    for(auto i = 0; i < left.size(); i++)
    {
        lhs_vertex_map[left[i]] = right[i];
        std::cout << "{ " << left[i] << " -> " << mid[i] << " -> " << right[i] << " }\n";
    }

    auto& rewrite = gamma_analysis.with_rewrites.at(rname);
    rewrite.print_node_sets(rname);
    std::cout << "\n";
    rewrite.print_edge_sets(rname);

    auto lhs_match = YAGL::induced_subgraph(system_graph, right);
    auto lhs_match_copy = lhs_match;
    auto rhs_rule_copy = gamma_analysis.with_rules.at(rname).rhs_graph;

    //we can build rhs of the vertex map by making a copy of the lhs and deleting
    auto rhs_vertex_map = lhs_vertex_map;

    for(auto& k : rewrite.node_set_destroy)
    {
        rhs_vertex_map.erase(k);
        lhs_match_copy.removeNode(lhs_match_copy.findNode(lhs_vertex_map[k])->second);
    }

    std::cout << lhs_match_copy << "\n";
    REQUIRE(YAGL::connected_components(lhs_match_copy) == 1); //nothing has changed
    REQUIRE(rhs_vertex_map.size() == 2);

    for(auto& k : rewrite.node_set_create)
    {
        //doing this way helps cheese my way into creating the correct types
        auto n = rhs_rule_copy.findNode(k)->second;
        decltype(n) node(gen.getKey(), n.getData());
        lhs_match_copy.addNode(node);
        rhs_vertex_map[k] = node.getKey();
    }
    std::cout << lhs_match_copy << "\n";
    REQUIRE(YAGL::connected_components(lhs_match_copy) == 2); //new node added
    REQUIRE(rhs_vertex_map.size() == 3);

    for(auto& [k1, k2] : rhs_vertex_map)
        std::cout << k1 << " " << k2 << "\n";

    for(auto& [u, v] : rewrite.edge_set_destroy)
    {
        //need to use the removal list edges and map those to the correct edges for the match
        lhs_match_copy.removeEdge(lhs_vertex_map[u], lhs_vertex_map[v]);
    }
    std::cout << lhs_match_copy << "\n";
    REQUIRE(YAGL::connected_components(lhs_match_copy) == 3); //edge removed => 3 disconnected nodes

    for(auto& [u, v] : rewrite.edge_set_create)
    {
        // I think I need to have a vertex mapping for the rhs, which is incomplete until
        // all new nodes are created since their keys are uniquely generated
        lhs_match_copy.addEdge(rhs_vertex_map[u], rhs_vertex_map[v]);
    }
    for(auto& [k1, k2] : rhs_vertex_map)
        std::cout << k1 << " " << k2 << "\n";
    std::cout << lhs_match_copy << "\n";
    REQUIRE(YAGL::connected_components(lhs_match_copy) == 1); //two new edges added to connect up the nodes

    //TODO: I also need to do some work to make sure any nodes that change only type
    // are changed. This is because the node rewrite set currently finds set difference
    // on keys not types, below is a first attempt, but there may be a more efficient way
    for(auto& [k, v] : lhs_match_copy.getNodeSetRef())
        std::cout << v.getData().type << "\n";
    for(auto& [k, _] : rhs_rule_copy.getNodeSetRef())
    {
        auto& v1 = lhs_match_copy[rhs_vertex_map[k]];

        auto& v2 = rhs_rule_copy[k];

        if(v1.type != v2.type)
        {
            v1.setData(v2.data);
            //copy assignment doesn't work since it won't invoke updating type
            //v1.data = v2.data;
        }
    }
    std::cout << "Here\n";
    for(auto& [k, v] : lhs_match_copy.getNodeSetRef())
        std::cout << v.getData().type << "\n";

    //I think we actually need a map for the lhs, and the rhs
    gamma_analysis.with_rules.at(rname).update(lhs_match, lhs_match_copy, lhs_vertex_map, rhs_vertex_map); //h

    std::cout << system_graph << "\n";
    //update the system graph
    for(auto& k : rewrite.node_set_destroy)
        system_graph.removeNode(system_graph.findNode(lhs_vertex_map[k])->second);
    std::cout << system_graph << "\n";
    for(auto& [u, v] : rewrite.edge_set_destroy)
        system_graph.removeEdge(lhs_vertex_map[u], lhs_vertex_map[v]);
    std::cout << system_graph << "\n";

    for(auto& [k, v] : lhs_match_copy.getNodeSetRef())
    {
        auto res = system_graph.findNode(k);
        if(res != system_graph.node_list_end())
        {
            std::cout << "Node " << k << " was found\n";
            //just update the data
            res->second = v;
        }
        else
        {
            std::cout << "Node " << k << " was not found\n";
            system_graph.addNode(v);
        }
    }
    std::cout << system_graph << "\n";
    for(auto& [u, v] : rewrite.edge_set_create)
        system_graph.addEdge(rhs_vertex_map[u], rhs_vertex_map[v]);
    std::cout << system_graph << "\n";
    REQUIRE(YAGL::connected_components(system_graph) == 2);
}

TEST_CASE("Basic Rewrite Test", "[basic-rewrite-test]")
{
    std::cout << "\n" << "Running basic rewrite test" << "\n";
    //build grammar and analyze it
    DGGML::Grammar<GraphType> gamma;
    define_model(gamma);
    DGGML::AnalyzedGrammar<GraphType> gamma_analysis(gamma);
    gamma.print();

    //initialize the system graph
    GraphType system_graph;
    initialize_system(system_graph);
    REQUIRE(system_graph.numNodes() == 6);
    REQUIRE(system_graph.numEdges() == 4);
    REQUIRE(YAGL::connected_components(system_graph) == 2);

    //only one unique component, so we can do this
    auto c = gamma_analysis.unique_components.begin()->second;
    auto results = YAGL::subgraph_isomorphism2(c, system_graph);
    REQUIRE(results.size() == 2);

    //verify the matches by output
    for(auto i = 0; i < results.size(); i++)
    {
        std::cout << "\nMatch " << i << ":\n";
        for(auto& [k1, k2] : results[i])
            std::cout << "{ " << k1 << " <- " << k2 << " }\n";
    }

    std::map<std::size_t, DGGML::ComponentMatch<std::size_t>> component_match_set;
    int k = 0;
    for(auto& item : results)
    {
        std::vector<std::size_t> match;
        for(auto& [key, value] : item)
        {
            match.emplace_back(value);
        }
        DGGML::ComponentMatch<std::size_t> inst;
        inst.match = match;
        inst.type = 0;
        inst.anchor = match[0];
        component_match_set.insert({k++, inst});
        //rule_system.push_back(inst);
    }

    std::map<std::size_t, DGGML::RuleMatch<std::size_t>> rule_matches;
    for(int i = 0; i < component_match_set.size(); i++) {
        rule_matches[i];
        rule_matches[i].name = "with_growth";
        rule_matches[i].category = "stochastic";
        rule_matches[i].components.push_back(i);
        rule_matches[i].anchor = component_match_set[i].anchor;
    }

    //we'll select the first match to use for our rewrite
    auto m1 = results[0];

    //simple key generator for newly created nodes, in the main code, something like it should be used to ensure that
    //keys are always generated uniquely in all phases of the code
    KeyGenerator gen(200);

    perform_rewrite(rule_matches[0], component_match_set, gen, gamma_analysis, system_graph);

//    auto ccuvm = gamma_analysis.ccuv_mappings["with_growth"][0];
//
//    //print out the mapping info, so we know how a lhs numbering maps to an instance numbering
//    std::cout << "\nMappings: { LHS Key -> Minimal Component Key -> Rule ComponentMatch Key }\n";
//    for(auto& [k, v] : gamma_analysis.lhs_connected_components.at("with_growth"))
//    {
//        for(auto& [id, n] : v.getNodeSetRef())
//        {
//            std::cout << "{ " << id << " -> " << ccuvm[id] << " -> " << m1[ccuvm[id]] << " }\n";
//        }
//    }
//
//    std::cout << "\n" << "Beginning rewrite code ..." << "\n";
//
//    auto& rewrite = gamma_analysis.with_rewrites.at("with_growth");
//    rewrite.print_node_sets("with_growth");
//    std::cout << "\n";
//    rewrite.print_edge_sets("with_growth");
//
//    std::vector<std::size_t> inst_keys; for(auto& [k, v] : m1) inst_keys.push_back(v);
//    auto lhs_match = YAGL::induced_subgraph(system_graph, inst_keys); //picks up extra edges, but should be fine
//    auto lhs_match_copy = lhs_match;
//    auto rhs_rule_copy = gamma_analysis.with_rules.at("with_growth").rhs_graph;
//
//    //simple key generator for newly created nodes, in the main code, something like it should be used to ensure that
//    //keys are always generated uniquely in all phases of the code
//    KeyGenerator gen(200);
//
//    for(auto& k : rewrite.node_set_destroy)
//        lhs_match_copy.removeNode(lhs_match_copy.findNode(m1[ccuvm[k]])->second);
//
//    std::cout << lhs_match_copy << "\n";
//    REQUIRE(YAGL::connected_components(lhs_match_copy) == 1); //nothing has changed
//
//    decltype(ccuvm) rhs_vertex_map;
//
//    //hacky way to build, start with what's left after removal, and add in what's next
//    for(auto& [k, v] : lhs_match_copy.getNodeSetRef()) {
//        for(auto& item : m1) {
//            if (item.second == k)
//            {
//                for(auto& jtem: ccuvm)
//                {
//                    if(jtem.second == item.first)
//                    {
//                        rhs_vertex_map[jtem.first] = k;
//                        break;
//                    }
//                }
//                break;
//            }
//        }
//    }
//    for(auto& [k1, k2] : rhs_vertex_map)
//        std::cout << k1 << " " << k2 << "\n";
//
//    for(auto& k : rewrite.node_set_create)
//    {
//        //doing this way helps cheese my way into creating the correct types
//        auto n = rhs_rule_copy.findNode(k)->second;
//        decltype(n) node(gen.getKey(), n.getData());
//        lhs_match_copy.addNode(node);
//        rhs_vertex_map[k] = node.getKey();
//    }
//    std::cout << lhs_match_copy << "\n";
//    REQUIRE(YAGL::connected_components(lhs_match_copy) == 2); //new node added
//
//    for(auto& [k1, k2] : rhs_vertex_map)
//        std::cout << k1 << " " << k2 << "\n";
//
//    for(auto& [u, v] : rewrite.edge_set_destroy)
//    {
//        //need to use the removal list edges and map those to the correct edges for the match
//        lhs_match_copy.removeEdge(m1[ccuvm[u]], m1[ccuvm[v]]);
//    }
//    std::cout << lhs_match_copy << "\n";
//    REQUIRE(YAGL::connected_components(lhs_match_copy) == 3); //edge removed => 3 disconnected nodes
//
//    for(auto& [u, v] : rewrite.edge_set_create)
//    {
//        // I think I need to have a vertex mapping for the rhs, which is incomplete until
//        // all new nodes are created since their keys are uniquely generated
//        lhs_match_copy.addEdge(rhs_vertex_map[u], rhs_vertex_map[v]);
//    }
//    for(auto& [k1, k2] : rhs_vertex_map)
//        std::cout << k1 << " " << k2 << "\n";
//    std::cout << lhs_match_copy << "\n";
//    REQUIRE(YAGL::connected_components(lhs_match_copy) == 1); //two new edges added to connect up the nodes
//
//    //TODO: I also need to do some work to make sure any nodes that change only type
//    // are changed. This is because the node rewrite set currently finds set difference
//    // on keys not types, below is a first attempt, but there may be a more efficient way
//    for(auto& [k, v] : lhs_match_copy.getNodeSetRef())
//        std::cout << v.getData().type << "\n";
//    for(auto& [k, _] : rhs_rule_copy.getNodeSetRef())
//    {
//        auto& v1 = lhs_match_copy[rhs_vertex_map[k]];
//
//        auto& v2 = rhs_rule_copy[k];
//
//        if(v1.type != v2.type)
//        {
//            v1.setData(v2.data);
//            //copy assignment doesn't work since it won't invoke updateing type
//            //v1.data = v2.data;
//        }
//    }
//    for(auto& [k, v] : lhs_match_copy.getNodeSetRef())
//        std::cout << v.getData().type << "\n";
//
//    //compose maps
//    decltype(ccuvm) h;
//    for(auto& [k1, v1] : ccuvm)
//    {
//        h[k1] = m1[v1];
//    }
//    for(auto& [k1, v1] : h)
//        std::cout << k1 << " --> " << v1 << "\n";
//    gamma_analysis.with_rules.at("with_growth").update(lhs_match, lhs_match_copy, h);//rhs_vertex_map);
//
//    std::cout << system_graph << "\n";
//    //update the system graph
//    for(auto& k : rewrite.node_set_destroy)
//        system_graph.removeNode(system_graph.findNode(m1[ccuvm[k]])->second);
//    std::cout << system_graph << "\n";
//    for(auto& [u, v] : rewrite.edge_set_destroy)
//        system_graph.removeEdge(m1[ccuvm[u]], m1[ccuvm[v]]);
//    std::cout << system_graph << "\n";
//
//    for(auto& [k, v] : lhs_match_copy.getNodeSetRef())
//    {
//        auto res = system_graph.findNode(k);
//        if(res != system_graph.node_list_end())
//        {
//            std::cout << "Node " << k << " was found\n";
//            //just update the data
//            res->second = v;
//        }
//        else
//        {
//            std::cout << "Node " << k << " was not found\n";
//            system_graph.addNode(v);
//        }
//    }
//    std::cout << system_graph << "\n";
//    for(auto& [u, v] : rewrite.edge_set_create)
//        system_graph.addEdge(rhs_vertex_map[u], rhs_vertex_map[v]);
//    std::cout << system_graph << "\n";
//    REQUIRE(YAGL::connected_components(system_graph) == 2);

}

TEST_CASE("Interaction Rewrite Test", "[basic-rewrite-test]")
{
    std::cout << "\n" << "Running interaction rewrite test" << "\n";
    //build grammar and analyze it
    DGGML::Grammar<GraphType> gamma;
    define_model(gamma);
    DGGML::AnalyzedGrammar<GraphType> gamma_analysis(gamma);
    gamma.print();

    //initialize the system graph
    GraphType system_graph;
    initialize_system(system_graph);
    REQUIRE(system_graph.numNodes() == 6);
    REQUIRE(system_graph.numEdges() == 4);
    REQUIRE(YAGL::connected_components(system_graph) == 2);

    //find matches
    auto c1 = gamma_analysis.unique_components[0];
    std::cout << c1.numNodes() << "\n";
    auto end_matches = YAGL::subgraph_isomorphism2(c1, system_graph);
    REQUIRE(end_matches.size() == 2);
    auto c2 = gamma_analysis.unique_components[1];
    std::cout << c2.numNodes() << "\n";
    auto mt_matches = YAGL::subgraph_isomorphism2(c2, system_graph);
    REQUIRE(mt_matches.size() == 2);

    std::map<std::size_t, DGGML::ComponentMatch<std::size_t>> component_match_set;
    int k = 0;
    for(auto& item : end_matches)
    {
        std::vector<std::size_t> match;
        for(auto& [key, value] : item)
        {
            match.emplace_back(value);
        }
        DGGML::ComponentMatch<std::size_t> inst;
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
        DGGML::ComponentMatch<std::size_t> inst;
        inst.match = match;
        inst.type = 1;
        inst.anchor = match[0];
        component_match_set.insert({k++, inst});
        //rule_system.push_back(inst);
    }

    //building the rule instances
    int kk = 0;
    std::map<std::size_t, DGGML::RuleMatch<std::size_t>> rule_matches;
    for(int i = 0; i < component_match_set.size(); i++)
    {
        if(component_match_set[i].type == 0) {
            rule_matches[kk];
            rule_matches[kk].name = "with_growth";
            rule_matches[kk].category = "stochastic";
            rule_matches[kk].components.push_back(i);
            rule_matches[kk].anchor = component_match_set[i].anchor;
            kk++;
        }
        for(int j = 0; j < component_match_set.size(); j++)
        {
            if(j <= i ) continue;
            auto& comp1 = component_match_set[i];
            auto& comp2 = component_match_set[j];
            auto& p1 = system_graph.findNode(comp1.anchor)->second.getData().position;
            auto& p2 = system_graph.findNode(comp2.anchor)->second.getData().position;
            //auto d = DGGML::calculate_distance(p1, p2);
            if(comp1.type == 1 && comp2.type == 1)// && d <= 1.0)
            {
                rule_matches[kk];
                rule_matches[kk].name = "with_interaction";
                rule_matches[kk].category = "stochastic";
                rule_matches[kk].components.push_back(i);
                rule_matches[kk].components.push_back(j);
                rule_matches[kk].anchor = comp1.anchor;
                kk++;
            }
        }
    }

    //simple key generator for newly created nodes, in the main code, something like it should be used to ensure that
    //keys are always generated uniquely in all phases of the code
    DGGML::KeyGenerator<std::size_t> gen(200);
    std::cout << "Here\n";
    for(auto& r : rule_matches) {
        std::cout << r.second.name << ": { ";
        for (auto &item: r.second.components) std::cout << item << " ";
        std::cout << "}\n";
    }
    DGGML::perform_invalidations_and_rewrite(rule_matches[1], component_match_set, gen, gamma_analysis, system_graph);
    std::cout << system_graph << "\n";
    std::cout << YAGL::connected_components(system_graph) << "\n";
    //REQUIRE(system_graph.numNodes() == 0);
}
