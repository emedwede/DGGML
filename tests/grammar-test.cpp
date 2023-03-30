#include <iostream>

#include "catch.hpp"

#include "YAGL_Graph.hpp"
#include "YAGL_Algorithms.hpp"
#include "MemoryManager.hpp"
#include "PlantTypes.hpp"
#include "PlantGrammar.hpp"

#include <unordered_map> 

struct MyInterface {};

using key_type = Cajete::Plant::mt_key_type; 
using data_type = Cajete::Plant::MT_NodeData;

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

using grammar_type = std::map<std::string, LHS>;

TEST_CASE("Grammar Test", "[grammar-test]")
{
    //map for the grammar
    grammar_type gamma; 
   
    //graph for a growing MT
    graph_type g1;
    g1.addNode({0, {{0, 0, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    g1.addNode({1, {{0, 0, 0}, {0, 0, 0}, Cajete::Plant::positive}});
    g1.addEdge(0, 1);

    //graph for a retraction MT
    graph_type g2;
    g2.addNode({0, {{0, 0, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    g2.addNode({1, {{0, 0, 0}, {0, 0, 0}, Cajete::Plant::negative}});
    g2.addEdge(0, 1);

    //graph for a collision
    graph_type g3;
    g3.addNode({0, {{0, 0, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    g3.addNode({1, {{0, 0, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    g3.addNode({2, {{0, 0, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    g3.addNode({3, {{0, 0, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    g3.addNode({4, {{0, 0, 0}, {0, 0, 0}, Cajete::Plant::negative}});
    g3.addEdge(0, 1);
    g3.addEdge(1, 2);
    g3.addEdge(3, 4);

    REQUIRE(g1.numNodes() == 2);
    REQUIRE(g1.numEdges() == 1);
    REQUIRE(YAGL::connected_components(g1) == 1);

    REQUIRE(g2.numNodes() == 2);
    REQUIRE(g2.numEdges() == 1);
    REQUIRE(YAGL::connected_components(g2) == 1);

    REQUIRE(g3.numNodes() == 5);
    REQUIRE(g3.numEdges() == 3);
    REQUIRE(YAGL::connected_components(g3) == 2);
    
    auto r1 = LHS(g1);
    auto r2 = LHS(g2);
    auto r3 = LHS(g3);

    gamma.insert({"growing", r1});
    gamma.insert({"retraction", r2});
    gamma.insert({"collision", r3});

    //check the number of rules 
    REQUIRE(gamma.size() == 3); 

    REQUIRE(gamma["growing"].components.size() == 1);
    REQUIRE(gamma["retraction"].components.size() == 1);
    REQUIRE(gamma["collision"].components.size() == 2);
    
    std::vector<graph_type> minimal_set;
    
    //add a component to start comparing 
    minimal_set.push_back(gamma["growing"].components[0]);

    //find the minimal_set of components
    for(auto& [key, value] : gamma)
    {
       for(auto& c : value.components)
       {
            bool found = false;
            for(auto& s : minimal_set)
            {
                auto matches = YAGL::subgraph_isomorphism2(c, s);
                //not there so we should add it
                if(!matches.empty()) 
                {
                    found = true;
                    break;
                }
            }
            if(!found)
            {
                minimal_set.push_back(c);
            }
        }
    }

    REQUIRE(minimal_set.size() == 3);
}

