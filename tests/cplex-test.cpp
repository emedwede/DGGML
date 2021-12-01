#include <iostream>

#include "catch.hpp"

#include "CellComplex.hpp"

#include <iostream>

TEST_CASE("Cell Complex from Cartesian 2D Grid Test", "[cplex test]")
{    
    std::cout << "Running the Cell Complex from Cartesian 2D Grid Test\n";

    double min_x = 0.0, min_y = 0.0, max_x = 3.0, max_y = 6.0,
           delta_x = 1.5, delta_y = 2.0;
    
    Cajete::BrickGrid2D grid(min_x, min_y, max_x, max_y, delta_x, delta_y);
    
    std::cout << grid;
    
    double fat_rad = 0.1;

    Cajete::CellComplex2D cplex2D(grid, fat_rad);
}
