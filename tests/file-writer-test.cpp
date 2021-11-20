#include <iostream>

#include "catch.hpp"

#include "../src/VtkWriter.hpp"
#include "../src/MemoryManager.hpp"
#include <string>

struct MyInterface {};

TEST_CASE("Running Vtk FileWriter Test", "[vtk test]")
{   
    using data_type = float;
    using writer_type = Cajete::VtkFileWriter<data_type>;

    writer_type* vtk_writer = 
        Cajete::MemoryManager::allocate_std<writer_type>(1);
    
    data_type data = 2.6;

    std::string filename = "test_step_0";    
    
    vtk_writer->save(data, filename);

    Cajete::MemoryManager::deallocate_std(vtk_writer);
}
