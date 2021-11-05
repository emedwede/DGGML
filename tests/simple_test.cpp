#include <iostream>

#include "catch.hpp"

TEST_CASE("Simple Test", "[simple test]")
{    
    std::cout << "Running simple test\n";
    REQUIRE(1 == 1);
}
