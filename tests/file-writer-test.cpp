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
    using writer_type = Cajete::VtkFileWriter<graph_type>;

    writer_type* vtk_writer = 
        Cajete::MemoryManager::allocate_std<writer_type>(1);
    
    graph_type graph;

    std::string filename = "test_step_0";    
    
    vtk_writer->save(graph, filename);

    Cajete::MemoryManager::deallocate_std(vtk_writer);
}
