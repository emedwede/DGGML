#include <iostream>

#include "catch.hpp"

#include "YAGL_Graph.hpp"
#include "YAGL_Algorithms.hpp"
#include "VtkWriter.hpp"
#include "MemoryManager.hpp"
#include "PlantTypes.hpp"
#include "PlantGrammar.hpp"
#include "ExpandedComplex2D.hpp"

#include <vector> 
#include <set>
#include <unordered_map> 
#include <string>

struct MyInterface {};

using key_type = Cajete::Plant::mt_key_type; 
using data_type = Cajete::Plant::MT_NodeData;

using graph_type = YAGL::Graph<key_type, data_type>;
using node_type = YAGL::Node<key_type, data_type>;

using writer_type = Cajete::VtkFileWriter<graph_type>;

void initialize_plant_test_graph(graph_type& graph) 
{
    graph.addNode({1, {{8, 6, 0}, {0, 0, 0}, Cajete::Plant::positive}});
    graph.addNode({2, {{8, 5, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({3, {{5, 5, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({4, {{4, 4, 0}, {0, 0, 0}, Cajete::Plant::junction}});
    graph.addNode({5, {{3, 5, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({6, {{3, 6, 0}, {0, 0, 0}, Cajete::Plant::negative}});
    graph.addNode({7, {{3, 4, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({8, {{2, 3, 0}, {0, 0, 0}, Cajete::Plant::negative}});
    graph.addNode({9, {{1, 5, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({10, {{0, 4, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({11, {{0, 3, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({12, {{0, 1, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({13, {{1, 0, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({14, {{3, 1, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({15, {{4, 2, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({16, {{6, 3, 0}, {0, 0, 0}, Cajete::Plant::junction}});
    graph.addNode({17, {{5, 4, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({18, {{7, 4, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({19, {{9, 4, 0}, {0, 0, 0}, Cajete::Plant::positive}});
    graph.addNode({20, {{7, 2, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({21, {{9, 2, 0}, {0, 0, 0}, Cajete::Plant::negative}});
    
    graph.addEdge(1,2);
    graph.addEdge(2,3);
    graph.addEdge(3,4);
    graph.addEdge(4,5);
    graph.addEdge(5,6);
    graph.addEdge(4,7);
    graph.addEdge(7,8);
    graph.addEdge(7,9);
    graph.addEdge(9,10);
    graph.addEdge(10,11);
    graph.addEdge(11,12);
    graph.addEdge(12,13);
    graph.addEdge(13,14);
    graph.addEdge(14,15);
    graph.addEdge(15,16);
    graph.addEdge(16,17);
    graph.addEdge(17,4);
    graph.addEdge(16,18);
    graph.addEdge(18,19);
    graph.addEdge(16,20);
    graph.addEdge(20,21);
}

void print_matches(std::vector<std::vector<key_type>>& matches)
{
    std::cout << "Found " << matches.size() << " Matches: \n";
    for(auto& m : matches)
    {
        std::cout << "\t{ ";
        for(auto& v : m)
        {
            std::cout << v << " ";
        } std::cout << "}\n";
    }
}

void create_junction_graph(graph_type& graph)
{
    //graph.addNode()
    graph.addNode({1, {{0, 0, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({2, {{0, 0, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({3, {{0, 0, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({4, {{0, 0, 0}, {0, 0, 0}, Cajete::Plant::intermediate}});
    graph.addNode({5, {{0, 0, 0}, {0, 0, 0}, Cajete::Plant::junction}});
    graph.addNode({6, {{0, 0, 0}, {0, 0, 0}, Cajete::Plant::junction}});

    graph.addEdge(5, 1);
    graph.addEdge(5, 2);
    graph.addEdge(5, 3);
    graph.addEdge(5, 4);
    graph.addEdge(1, 6);

}
auto heuristic_search_junctions_extended(graph_type& graph)
{
    std::vector<std::vector<key_type>> matches;

    for(auto iter = graph.node_list_begin(); iter != graph.node_list_end(); iter++)
    {
        auto i = iter->first;
        auto itype = iter->second.getData().type;
        if(itype == Cajete::Plant::junction)
        {
            for(auto jter = graph.out_neighbors_begin(i); jter != graph.out_neighbors_end(i); jter++)
            {
                auto j = *jter;
                auto jtype = graph.findNode(j)->second.getData().type; 
                if(jtype == Cajete::Plant::intermediate)
                {
                    for(auto kter = graph.out_neighbors_begin(j);
                        kter != graph.out_neighbors_end(j); kter++)
                    {
                        auto k = *kter;
                        auto ktype = graph.findNode(k)->second.getData().type; 
                        if(k != i)
                        {
                            if(ktype == Cajete::Plant::junction)
                            {
                                for(auto lter = graph.out_neighbors_begin(i);
                                         lter != graph.out_neighbors_end(i); lter++)
                                {
                                    auto l = *lter;
                                    auto ltype = graph.findNode(l)->second.getData().type; 
                                    if(l != j)
                                    {
                                        if(ltype == Cajete::Plant::intermediate)
                                        {
                                            for(auto mter = graph.out_neighbors_begin(i);
                                                    mter != graph.out_neighbors_end(i); mter++)
                                            {
                                                auto m = *mter;
                                                auto mtype = graph.findNode(m)->second.getData().type;
                                                if(m != j && m != l)
                                                {
                                                    if(mtype == Cajete::Plant::intermediate)
                                                    {
                                                        for(auto nter =graph.out_neighbors_begin(i);
                                                                nter != graph.out_neighbors_end(i);
                                                                nter++)
                                                        {
                                                            auto n = *nter;
                                                            auto ntype =
                                                                graph.findNode(n)->second.getData().type;
                                                            if(n != j && n != l && n != m)
                                                            {
                                                                if(ntype == Cajete::Plant::intermediate)
                                                                {
                                                                    std::vector<key_type> found;
                                                                    found.push_back(i);
                                                                    found.push_back(j);
                                                                    found.push_back(k);
                                                                    found.push_back(l);
                                                                    found.push_back(m);
                                                                    found.push_back(n);
                                                                    matches.push_back(found);
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return matches;
}

template<typename GraphType, typename BucketType>
void push_to_bucket(GraphType& graph, BucketType& bucket)
{
    for(auto iter = graph.node_list_begin(); iter != graph.node_list_end(); iter++)
    {
        bucket.push_back(iter->first);
    }
}

TEST_CASE("Heuristic Growing End Matcher Test", "[heuristic-matcher-test]")
{
    graph_type graph;
    
    initialize_plant_test_graph(graph);
    
    std::vector<Cajete::Plant::mt_key_type> bucket;

    push_to_bucket(graph, bucket);

    auto matches = Cajete::Plant::microtubule_growing_end_matcher(graph, bucket);
    print_matches(matches);
    REQUIRE(matches.size() == 2);
    
    auto old_node_size = graph.numNodes();
    auto old_edge_size = graph.numEdges();

    for(auto match : matches)
    {
        Cajete::Plant::microtubule_growing_end_polymerize_rewrite(graph, match);
    }

    REQUIRE(graph.numNodes() == old_node_size + 2); //since we add a split node
    REQUIRE(graph.numEdges() == old_edge_size + 2); //since we split an edge

    writer_type vtk_writer;

    vtk_writer.save(graph, "polymerize_rewrite");

}

TEST_CASE("Heuristic Retraction End Matcher Test", "[heuristic-matcher-test]")
{
    graph_type graph;
    initialize_plant_test_graph(graph);
    
    std::vector<Cajete::Plant::mt_key_type> bucket;

    push_to_bucket(graph, bucket);

    auto matches = Cajete::Plant::microtubule_retraction_end_matcher(graph, bucket);
    print_matches(matches);
    REQUIRE(matches.size() == 3);

    
}

TEST_CASE("Heuristic Wildcard Matcher Test", "[heuristic-matcher-test]")
{
    graph_type graph;
    initialize_plant_test_graph(graph);

    std::vector<Cajete::Plant::mt_key_type> bucket;

    push_to_bucket(graph, bucket);

    auto matches = Cajete::Plant::wildcard_intermediate_wildcard_matcher(graph, bucket);
    print_matches(matches);

    REQUIRE(matches.size() == 32);
}

TEST_CASE("Heuristic Extended Junction Matcher Test", "[heuristic-matcher-test]")
{  
   
    writer_type* vtk_writer = 
        Cajete::MemoryManager::allocate_std<writer_type>(1);
    
    graph_type graph;
    
    initialize_plant_test_graph(graph);
    
    std::string filename = "plant_matcher_graph_0";    
    
    vtk_writer->save(graph, filename);
 
    auto matches = heuristic_search_junctions_extended(graph);
    print_matches(matches);

    REQUIRE(matches.size() == 12);
    
    //For fun check the number of connected components
    REQUIRE(YAGL::connected_components(graph) == 1); 

    //remove a node to ensure we match properly
    auto node = graph.findNode(7)->second;
    graph.removeNode(node);
    node = graph.findNode(8)->second;
    graph.removeNode(node);
    filename = "plant_matcher_graph_1";     
    vtk_writer->save(graph, filename);

    matches = heuristic_search_junctions_extended(graph);
    print_matches(matches);
    REQUIRE(matches.size() == 6);
    
    //For fun check the number of connected components
    REQUIRE(YAGL::connected_components(graph) == 1); 

    //finally remove node 17 to ensure no matches
    node = graph.findNode(17)->second;
    graph.removeNode(node);
    filename = "plant_matcher_graph_2";     
    vtk_writer->save(graph, filename);

    matches = heuristic_search_junctions_extended(graph);
    print_matches(matches);
    REQUIRE(matches.size() == 0);
    
    // removing node 17 should disconnect the graph
    REQUIRE(YAGL::connected_components(graph) == 2); 
    
    Cajete::MemoryManager::deallocate_std(vtk_writer);

}

TEST_CASE("Labeled Subgraph Isomorphism Matcher Test", "[subgraph-iso-matcher-test]")
{
    graph_type g1, g2;
    
    create_junction_graph(g1);
    initialize_plant_test_graph(g2);
   
    //should find 12 permutations, so the 6 automorphisms for each match
    auto matches = YAGL::subgraph_isomorphism(g1, g2);
    
    REQUIRE(matches.size() == 12); 
    
    //sloppy way of ordering the permutations and hashing them to 
    //find the number of unique matches
    std::set<std::string> match_count;
    for(auto& item : matches)
    {
        std::set<key_type> ordering;
        std::cout << "Permutation : { ";
        for(auto& [key, value] : item)
        {
            ordering.insert(value);
            std::cout << "[ "<< key << " -> " << value << "] ";
        } std::cout << " } \n";
        std::string str{"+"};
        for(auto& s : ordering)
        {
            str += std::to_string(s) + "+";
        }
        match_count.insert(str);
    }
    
    // there should be two unique matches
    REQUIRE(match_count.size() == 2);

    std::cout << "Unique matches: ";
    for(auto& s : match_count)
        std::cout << s << " ";
    std::cout << "\n";
}


