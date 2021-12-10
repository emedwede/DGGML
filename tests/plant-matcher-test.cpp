#include <iostream>

#include "catch.hpp"

#include "YAGL_Graph.hpp"

#include "VtkWriter.hpp"
#include "MemoryManager.hpp"
#include "PlantTypes.hpp"

#include <vector> 

#include <string>

struct MyInterface {};

using key_type = Cajete::Plant::mt_key_type; 
using data_type = Cajete::Plant::MT_NodeData;

using graph_type = YAGL::Graph<key_type, data_type>;
using node_type = YAGL::Node<key_type, data_type>;

using writer_type = Cajete::VtkFileWriter<graph_type>;

enum colors 
{
    red, 
    blue,
    green,
    pink
};

void initialize_plant_test_graph(graph_type& graph) 
{
    graph.addNode({1, {{8, 6, 0}, {0, 0, 0}, pink}});
    graph.addNode({2, {{8, 5, 0}, {0, 0, 0}, blue}});
    graph.addNode({3, {{5, 5, 0}, {0, 0, 0}, blue}});
    graph.addNode({4, {{4, 4, 0}, {0, 0, 0}, red}});
    graph.addNode({5, {{3, 5, 0}, {0, 0, 0}, blue}});
    graph.addNode({6, {{3, 6, 0}, {0, 0, 0}, green}});
    graph.addNode({7, {{3, 4, 0}, {0, 0, 0}, blue}});
    graph.addNode({8, {{2, 3, 0}, {0, 0, 0}, green}});
    graph.addNode({9, {{1, 5, 0}, {0, 0, 0}, blue}});
    graph.addNode({10, {{0, 4, 0}, {0, 0, 0}, blue}});
    graph.addNode({11, {{0, 3, 0}, {0, 0, 0}, blue}});
    graph.addNode({12, {{0, 1, 0}, {0, 0, 0}, blue}});
    graph.addNode({13, {{1, 0, 0}, {0, 0, 0}, blue}});
    graph.addNode({14, {{3, 1, 0}, {0, 0, 0}, blue}});
    graph.addNode({15, {{4, 2, 0}, {0, 0, 0}, blue}});
    graph.addNode({16, {{6, 3, 0}, {0, 0, 0}, red}});
    graph.addNode({17, {{5, 4, 0}, {0, 0, 0}, blue}});
    graph.addNode({18, {{7, 4, 0}, {0, 0, 0}, blue}});
    graph.addNode({19, {{9, 4, 0}, {0, 0, 0}, pink}});
    graph.addNode({20, {{7, 2, 0}, {0, 0, 0}, blue}});
    graph.addNode({21, {{9, 2, 0}, {0, 0, 0}, green}});
    
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
    for(auto m : matches)
    {
        std::cout << "\t{ ";
        for(auto v : m)
        {
            std::cout << v << " ";
        } std::cout << "}\n";
    }
}

auto heuristic_search_junctions_extended(graph_type& graph)
{
    std::vector<std::vector<key_type>> matches;

    for(auto iter = graph.node_list_begin(); iter != graph.node_list_end(); iter++)
    {
        auto i = iter->first;
        auto itype = iter->second.getData().type;
        if(itype == red)
        {
            for(auto jter = graph.out_neighbors_begin(i); jter != graph.out_neighbors_end(i); jter++)
            {
                auto j = *jter;
                auto jtype = graph.findNode(j)->second.getData().type; 
                if(jtype == blue)
                {
                    for(auto kter = graph.out_neighbors_begin(j);
                        kter != graph.out_neighbors_end(j); kter++)
                    {
                        auto k = *kter;
                        auto ktype = graph.findNode(k)->second.getData().type; 
                        if(k != i)
                        {
                            if(ktype == red)
                            {
                                for(auto lter = graph.out_neighbors_begin(i);
                                         lter != graph.out_neighbors_end(i); lter++)
                                {
                                    auto l = *lter;
                                    auto ltype = graph.findNode(l)->second.getData().type; 
                                    if(l != j)
                                    {
                                        if(ltype == blue)
                                        {
                                            for(auto mter = graph.out_neighbors_begin(i);
                                                    mter != graph.out_neighbors_end(i); mter++)
                                            {
                                                auto m = *mter;
                                                auto mtype = graph.findNode(m)->second.getData().type;
                                                if(m != j && m != l)
                                                {
                                                    if(mtype == blue)
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
                                                                if(ntype == blue)
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

TEST_CASE("Heuristic Extended Junction Matcher Test", "[heuristic-matcher-test]")
{  
   
    writer_type* vtk_writer = 
        Cajete::MemoryManager::allocate_std<writer_type>(1);
    
    graph_type graph;
    
    initialize_plant_test_graph(graph);
    
    std::string filename = "plant_matcher_graph";    
    
    vtk_writer->save(graph, filename);

    Cajete::MemoryManager::deallocate_std(vtk_writer);

    auto matches = heuristic_search_junctions_extended(graph);
    print_matches(matches);

    REQUIRE(matches.size() == 12);

}
