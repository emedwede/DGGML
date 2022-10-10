#include <iostream>

#include "catch.hpp"

#include "BrickGrid2D.hpp"
#include "CartesianGrid2D.hpp"
#include "VtkWriter.hpp"

#include <iostream>

TEST_CASE("Running Cartesian 2D Grid Test", "[grid test]")
{    
    std::cout << "Running the Cartesian 2D Grid Test\n";

    Cajete::CartesianGrid2D grid;
    
    double min_x = 0.0, min_y = 0.0, max_x = 4.0, max_y = 2.0;
    std::size_t nx = 2, ny = 2;

    grid.init(min_x, min_y, max_x, max_y, nx, ny);
    
    std::cout << grid << std::endl;

    REQUIRE( grid.totalNumCells() ==  nx*ny );
    REQUIRE (grid.totalNumPoints() == (nx+1)*(ny+1));

    int ic, jc, cell_id; double px, py;

    cell_id = 2;

    grid.ijCellIndex(cell_id, ic, jc);

    REQUIRE( ic == 0 );
    REQUIRE( jc == 1 );
    
    REQUIRE( grid.cardinalCellIndex(ic, jc) == cell_id );
    
    grid.cardinalCellToPoint(px, py, 2);

    REQUIRE( px == 1.0 );
    REQUIRE( py == 1.5 );

    px = 1.4; py = 1.2;

    grid.locatePoint(px, py, ic, jc);
    
    REQUIRE( grid.cardinalCellIndex(ic, jc) == cell_id );
    
    int lattice_id = 4;
    int ilp, jlp; //i and j lattice points indices 
    grid.ijLatticeIndex(lattice_id, ilp, jlp);
    
    REQUIRE( ilp == 1 );
    REQUIRE( jlp == 1 );
    
    REQUIRE( grid.cardinalLatticeIndex(ilp, jlp) == lattice_id);

    grid.cardinalLatticeToPoint(px, py, 4);

    REQUIRE( px == 2.0 );
    REQUIRE( py == 1.0 );
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

TEST_CASE("CartesianGrid2D Writer Test", "[grid writer test]")
{    
    std::cout << "Running the Cartesian 2D Grid Writer Test\n";

    Cajete::CartesianGrid2D grid;
    
    double min_x = 0.0, min_y = 0.0, max_x = 16.0, max_y = 9.0;
    std::size_t nx = 100, ny = 100;

    grid.init(min_x, min_y, max_x, max_y, nx, ny);
    std::cout << grid; 
    Cajete::GridFileWriter writer;
    writer.save(grid, "grid_viz");
}
 
