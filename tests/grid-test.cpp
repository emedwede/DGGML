#include <iostream>

#include "catch.hpp"

#include "BrickGrid2D.hpp"
#include "CartesianGrid2D.hpp"

#include <iostream>

TEST_CASE("Running Cartesian 2D Grid Test", "[grid test]")
{    
    std::cout << "Running the Cartesian 2D Grid Test\n";

    Cajete::CartesianGrid2D grid;
    
    double min_x = 0.0, min_y = 0.0, max_x = 4.0, max_y = 2.0,
           delta_x = 2.0, delta_y = 2.0;
    grid.init(min_x, min_y, max_x, max_y, delta_x, delta_y);

    REQUIRE( grid.totalNumCells() == 2 );
}

TEST_CASE("Running the Brick Grid 2D Test", "[grid test]")
{
    std::cout << "Running the 2D Brick Grid Test\n";

    double min_x = 0.0, min_y = 0.0, max_x = 4.0, max_y = 2.0,
           delta_x = 2.0, delta_y = 2.0;
    
    Cajete::BrickGrid2D grid(min_x, min_y, max_x, max_y, delta_y, delta_y);

    std::cout << "Total number of staggered grid cells: " << grid.totalNumCells() << "\n";
}
