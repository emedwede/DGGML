#include <iostream>

#include "catch.hpp"

#include "../src/VtkWriter.hpp"
#include "../src/MemoryManager.hpp"

struct MyInterface {};

TEST_CASE("Running Vtk FileWriter Test", "[vtk test]")
{   
    using writer_type = Cajete::VtkFileWriter;

    writer_type* vtk_writer = 
        Cajete::MemoryManager::allocate_std<writer_type>(1);

    vtk_writer->save();

    Cajete::MemoryManager::deallocate_std(vtk_writer);
}
