#include <iostream>

#include "catch.hpp"

#include "Utlities/VtkWriter.hpp"
#include "Utlities/MemoryManager.hpp"

#include "YAGL_Graph.hpp"

#include <string>

struct MyInterface {};

struct PhysicsData
{
    double position[3];
    int type;
};

TEST_CASE("Running Vtk FileWriter Test", "[vtk test]")
{  
    using key_type = std::size_t; using data_type = PhysicsData;
    using graph_type = YAGL::Graph<key_type, data_type>;
    using node_type = YAGL::Node<key_type, data_type>;
    using writer_type = DGGML::VtkFileWriter<graph_type>;

    writer_type* vtk_writer = 
        DGGML::MemoryManager::allocate_std<writer_type>(1);
    
    graph_type graph;
    
    // make a triangle graph
    node_type node_a(303, {{0.0, 0.0, 0.0}, 2});
    node_type node_b(101, {{0.0, 1.0, 0.0}, 2});
    node_type node_c(609, {{1.0, 1.0, 0.0}, 2});
    node_type node_d(506, {{3.0, 3.0, 3.0}, 2});    

    graph.addNode(node_a);
    graph.addNode(node_b);
    graph.addNode(node_c);
    graph.addNode(node_d);

    graph.addEdge(node_a, node_b);
    graph.addEdge(node_b, node_c);
    graph.addEdge(node_c, node_a);
    graph.addEdge(node_d, node_a);
    graph.addEdge(node_d, node_b);
    graph.addEdge(node_d, node_c);

    std::string filename = "test_step_0";    
    
    vtk_writer->save(graph, filename);

    DGGML::MemoryManager::deallocate_std(vtk_writer);
}
