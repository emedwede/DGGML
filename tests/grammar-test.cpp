#include <iostream>

#include "catch.hpp"

#include "YAGL_Graph.hpp"
#include "YAGL_Algorithms.hpp"
#include "MemoryManager.hpp"
#include "PlantTypes.hpp"
#include "PlantGrammar.hpp"

#include <unordered_map> 

struct MyInterface {};

using key_type = DGGML::Plant::mt_key_type;
using data_type = DGGML::Plant::MT_NodeData;

using graph_type = YAGL::Graph<key_type, data_type>;
using node_type = YAGL::Node<key_type, data_type>;

struct LHS
{
    graph_type graph;
    std::vector<graph_type> components;
    
    LHS() {}

    LHS(graph_type& graph) : graph(graph) 
    {
        //build the list of components 
        std::unordered_set<key_type> visited;
        std::size_t count = 0; //no connected components found to start
        for(auto i = graph.node_list_begin(); i != graph.node_list_end(); i++)
        {
            auto v = i->first;
            //node hasn't been visited so it must be the start of a new connected component 
            if(visited.find(v) == visited.end())
            {
                std::vector<key_type> path;
                //we could use whatever search method we feel like 
                YAGL::impl_iterative_bfs(graph, v, visited, path);
                components.push_back(YAGL::induced_subgraph(graph, path));
                count++;
            }
        }
    }
};

//should LHS and RHS be the same object?
struct RHS
{
    graph_type graph;
    std::vector<graph_type> components;
    
    RHS() {}

    RHS(graph_type& graph) : graph(graph) {}
};

struct Rule
{
    LHS lhs;
    RHS rhs;
    std::string name;
    
    Rule() {}

    Rule(std::string name, graph_type& lhs_graph, graph_type& rhs_graph) : name(name), lhs(lhs_graph), rhs(rhs_graph) {}
};

struct Grammar 
{
    std::map<std::string, Rule> rule_set;
    std::map<std::string, std::vector<std::size_t>> rule_component;
    std::map<std::size_t, std::vector<std::string>> component_rule;
    std::vector<graph_type> minimal_set;

    void addRule(Rule& r)
    {
        if(rule_set.find(r.name) == rule_set.end())
        {
            rule_set[r.name] = r;
            incremental_build_min_set(r);
        }
        else
            std::cout << "name already exists in rule_set\n";
    }

    void incremental_build_min_set(Rule& r)
    {
        
        for(auto& c : r.lhs.components)
        {
            auto idx = 0;
            bool found = false;
            for(auto i = 0; i < minimal_set.size(); i++)
            {
                auto matches = YAGL::subgraph_isomorphism2(c, minimal_set[i]);
                if(!matches.empty())
                {
                    idx = i;
                    found = true;
                    break;
                }
            }
            if(!found)
            {
                //add forward and backward references
                rule_component[r.name].push_back(minimal_set.size());
                component_rule[minimal_set.size()].push_back(r.name);
                minimal_set.push_back(c);
            }
            else //found
            {
                rule_component[r.name].push_back(idx);
                component_rule[idx].push_back(r.name);
            }
        }
    }

    void print() 
    {
        std::cout << "Rules: ";
        for(auto& rule : rule_set)
            std::cout << "\t" << rule.first << "\n";
    }
};

void print_mapping(Grammar& gamma)
{
    std::cout << "\nRule Component Mapppings:\n";
    for(auto& [key, value] : gamma.rule_set)
    {
        std::cout << "\t" << key << ": ";
        for(auto& item : gamma.rule_component[key])
        {
            std::cout << item << " ";
        } std::cout << "\n";
    }

    std::cout << "\nComponent Rule Mapppings:\n";
    for(auto key = 0; key < gamma.minimal_set.size(); key++)
    {
        std::cout << "\t" << key << ": ";
        for(auto& item : gamma.component_rule[key])
        {
            std::cout << item << " ";
        } std::cout << "\n";
    }

}

//this will need extensive error checking since it's user defined?
void define_model(Grammar& gamma) {
    //graph for a growing MT LHS
    graph_type g1;
    g1.addNode({0, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g1.addNode({1, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::positive}});
    g1.addEdge(0, 1);
    
    //graph for a growing MT RHS
    graph_type g2;
    g2.addNode({0, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g2.addNode({1, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g2.addNode({2, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::positive}});
    g2.addEdge(0, 1);
    g2.addEdge(1, 2);
    
    //create the grammar rule 
    Rule r1("growing", g1, g2);
    
    //graph for a catastrophe LHS
    graph_type g3;
    g3.addNode({0, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g3.addNode({1, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g3.addNode({2, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g3.addNode({3, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g3.addNode({4, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::positive}});
    g3.addEdge(0, 1);
    g3.addEdge(1, 2);
    g3.addEdge(3, 4);

    //graph for catastrophe RHS
    graph_type g4;
    g4.addNode({0, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g4.addNode({1, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g4.addNode({2, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g4.addNode({3, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g4.addNode({4, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::negative}});
    g4.addEdge(0, 1);
    g4.addEdge(1, 2);
    g4.addEdge(3, 4);
    
    //create the grammar rule 
    Rule r2("catastrophe", g3, g4);
    
    //graph for a retraction MT LHS
    graph_type g5;
    g5.addNode({0, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g5.addNode({1, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::negative}});
    g5.addEdge(0, 1);
    
    //graph for a retraction MT RHS
    graph_type g6;
    g6.addNode({0, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g6.addNode({1, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g6.addNode({2, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::negative}});
    g6.addEdge(0, 1);
    g6.addEdge(1, 2);
    
    //create the grammar rule 
    Rule r3("retraction", g5, g6);
    
    //graph for a zipper LHS
    graph_type g7;
    g7.addNode({0, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g7.addNode({1, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g7.addNode({2, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g7.addNode({3, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g7.addNode({4, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::positive}});
    g7.addEdge(0, 1);
    g7.addEdge(1, 2);
    g7.addEdge(3, 4);

    //graph for a zipper RHS
    graph_type g8;
    g8.addNode({0, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g8.addNode({1, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::zipper}});
    g8.addNode({2, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g8.addNode({3, {{0, 0, 0}, {0, 0, 0}, DGGML::Plant::intermediate}});
    g8.addEdge(0, 1);
    g8.addEdge(1, 2);
    g8.addEdge(3, 1);
    
    Rule r4("zipper", g7, g8);

    //build the Grammar 
    gamma.addRule(r1);
    gamma.addRule(r2);
    gamma.addRule(r3);
    gamma.addRule(r4);

    //if runtime just run 
    //else we need a two phase compile or some compile time way to generate optimized code?
}

using grammar_type = std::map<std::string, LHS>;

TEST_CASE("Grammar Test", "[grammar-test]")
{

    Grammar gamma;
    define_model(gamma);
    gamma.print();
    
    REQUIRE(gamma.rule_set["growing"].lhs.components.size() == 1);
    REQUIRE(gamma.rule_set["catastrophe"].lhs.components.size() == 2);
    REQUIRE(gamma.rule_set["retraction"].lhs.components.size() == 1); 
    REQUIRE(gamma.rule_set["zipper"].lhs.components.size() == 2); 

    REQUIRE(gamma.minimal_set.size() == 3);
    
    print_mapping(gamma);
}

