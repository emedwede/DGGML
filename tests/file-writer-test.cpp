#include <iostream>

#include "catch.hpp"

#include "VtkWriter.hpp"
#include "MemoryManager.hpp"

#include "YAGL_Graph.hpp"

#include <string>

struct MyInterface {};

TEST_CASE("Running Vtk FileWriter Test", "[vtk test]")
{  
    using key_type = std::size_t; using data_type = double;
    using graph_type = YAGL::Graph<key_type, data_type>;
    using node_type = YAGL::Node<key_type, data_type>;
    using writer_type = Cajete::VtkFileWriter<graph_type>;

    writer_type* vtk_writer = 
        Cajete::MemoryManager::allocate_std<writer_type>(1);
    
    graph_type graph;
    
    // make a triangle graph
    node_type node_a(303, 0.0);
    node_type node_b(101, 1.0);
    node_type node_c(609, 2.0);
    
    graph.addNode(node_a);
    graph.addNode(node_b);
    graph.addNode(node_c);
    graph.addEdge(node_a, node_b);
    graph.addEdge(node_b, node_c);
    graph.addEdge(node_c, node_a);

    std::string filename = "test_step_0";    
    
    vtk_writer->save(graph, filename);

    Cajete::MemoryManager::deallocate_std(vtk_writer);
}
