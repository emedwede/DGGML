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

    double min_x = 0.0, min_y = 0.0, max_x = 3.0, max_y = 6.0,
           delta_x = 1.5, delta_y = 2.0;
    
    Cajete::BrickGrid2D grid(min_x, min_y, max_x, max_y, delta_x, delta_y);
    
    std::cout << grid;
    
    REQUIRE( grid.compute_num_0D_zones() == 6 );
    REQUIRE( grid.compute_num_1D_zones() == 12 );
    REQUIRE( grid.totalNumCells() == 7 );

    int ic, jc;
    double px = 2.3, py = 2.3;
    
    grid.locatePoint(px, py, ic, jc);
    int cell_id = grid.cardinalCellIndex(ic, jc);

    REQUIRE( cell_id == 4 );

    px = 1.0; py = 1.0;
    grid.locatePoint(px, py, ic, jc);
    cell_id = grid.cardinalCellIndex(ic, jc);

    REQUIRE (cell_id == 0);
}
