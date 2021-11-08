#include <iostream>

#include "catch.hpp"
#include "../src/MemoryManager.hpp"

struct DataType {
    bool filled;
    int type;
    double energy;
    double position[3];
};

TEST_CASE(" Standard MemoryManager Host Test", "[standard memory host test]")
{    
    std::cout << "Running host allocation test (please run with a memory analysis tool)\n";
    
    std::size_t size;
    
    size = 1;
    auto ptr_a = Cajete::MemoryManager::allocate_std<int>(size);

    Cajete::MemoryManager::deallocate_std(ptr_a);

    REQUIRE(ptr_a == nullptr);

    size = 10;
    
    // Creates an array of structs (AoS)
    auto ptr_b = Cajete::MemoryManager::allocate_std<DataType>(size);

    Cajete::MemoryManager::deallocate_std(ptr_b);

    REQUIRE(ptr_b == nullptr);

}

TEST_CASE(" Smart MemoryManager Host Test", "[smart memory host test]")
{    
    std::cout << "Running host allocation test (please run with a memory analysis tool)\n";
    
    std::size_t size;
    
    size = 1;
    auto ptr_a = Cajete::MemoryManager::allocate_smart<int>(size);

    REQUIRE(ptr_a.use_count() == 1);

    size = 10;
    
    // Creates an array of structs (AoS)
    auto ptr_b = Cajete::MemoryManager::allocate_smart<DataType>(size);

    REQUIRE(ptr_b.use_count() == 1);

}
