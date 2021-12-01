#include <iostream>

#include "catch.hpp"

#include "BrickGrid2D.hpp"
#include "CartesianGrid2D.hpp"

#include <iostream>

TEST_CASE("Cell Complex from Cartesian 2D Grid Test", "[cplex test]")
{    
    std::cout << "Running the Cell Complex from Cartesian 2D Grid Test\n";

    Cajete::CartesianGrid2D grid;
    
    double min_x = 0.0, min_y = 0.0, max_x = 4.0, max_y = 2.0,
           delta_x = 2.0, delta_y = 2.0;
    grid.init(min_x, min_y, max_x, max_y, delta_x, delta_y);

}
