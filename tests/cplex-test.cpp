#include <iostream>

#include "catch.hpp"

#include "CartesianComplex2D.hpp"

#include "VtkWriter.hpp"

#include <iostream>

TEST_CASE("Cell Complex from Cartesian 2D Grid Test", "[cplex test]")
{    
    std::cout << "Running the Cell Complex from Cartesian 2D Grid Test\n";
    
    std::size_t n = 4, m = 3; //number of cells in the x and y direction
    double dx = 2.0, dy = 2.0; //size of each grid cell 

    Cajete::CartesianComplex2D cplex2D(n, m, dx, dy);

    std::cout << cplex2D;
    
    REQUIRE( cplex2D.get1dInteriorCellCount() == 17 );
    REQUIRE( cplex2D.get1dExteriorCellCount() == 14 );
    //Save the cell complex graph
    Cajete::VtkFileWriter<typename Cajete::CartesianComplex2D<>::graph_type> writer;
    
    writer.save(cplex2D.getGraph(), "cplex_graph");
}
