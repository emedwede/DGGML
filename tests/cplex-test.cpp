#include <iostream>

#include "catch.hpp"

#include "CartesianComplex2D.hpp"

#include "VtkWriter.hpp"

#include <iostream>

TEST_CASE("Cell Complex from Cartesian 2D Grid Test", "[cplex-test]")
{    
    std::cout << "Running the Cell Complex from Cartesian 2D Grid Test\n";
    
    std::size_t n = 4, m = 3; //number of cells in the x and y direction
    double dx = 2.0, dy = 2.0; //size of each grid cell 

    Cajete::CartesianComplex2D cplex2D(n, m, dx, dy);

    std::cout << cplex2D << std::endl;
    
    REQUIRE(cplex2D.get0dInteriorCellCount() == 6);
    REQUIRE(cplex2D.get0dExteriorCellCount() == 14);
    REQUIRE(cplex2D.get0dTotalCellCount() == 20);


    REQUIRE(cplex2D.get1dInteriorCellCount() == 17);
    REQUIRE(cplex2D.get1dExteriorCellCount() == 14);
    REQUIRE(cplex2D.get1dTotalCellCount() == 31);

    REQUIRE(cplex2D.get2dExteriorCellCount() == 0);
    REQUIRE(cplex2D.get2dInteriorCellCount() == 12);
    REQUIRE(cplex2D.get2dTotalCellCount() == 12);

    REQUIRE(cplex2D.getTotalCellCount() == 63);

    //Save the cell complex graph
    //Cajete::VtkFileWriter<typename Cajete::CartesianComplex2D<>::graph_type> writer;
    
    //writer.save(cplex2D.getGraph(), "cplex_graph");
}

TEST_CASE("Cell Complex can just be a single domain", "[cplex-test]")
{
    Cajete::CartesianComplex2D cplex2D(1, 1, 6.0, 4.0);
    
    std::cout << cplex2D << std::endl;
    
    REQUIRE(cplex2D.get0dInteriorCellCount() == 0);
    REQUIRE(cplex2D.get0dExteriorCellCount() == 4);
    REQUIRE(cplex2D.get0dTotalCellCount() == 4);


    //It's a rectangle and we have four edges on the boundary
    REQUIRE(cplex2D.get1dExteriorCellCount() == 4);
    REQUIRE(cplex2D.get1dInteriorCellCount() == 0);
    
    //Single rectangle means one 2D interior
    REQUIRE(cplex2D.get2dExteriorCellCount() == 0);
    REQUIRE(cplex2D.get2dInteriorCellCount() == 1);
    
    REQUIRE(cplex2D.getTotalCellCount() == 9);
}

TEST_CASE("Cell Complex can just be a single domain and have ghost cells", "[cplex-test]")
{
    Cajete::CartesianComplex2D cplex2D(1, 1, 6.0, 4.0, true);
    
    std::cout << cplex2D << std::endl;
   
    //Single rectangle means one 2D interior and 9 exterior boundary cells
    REQUIRE(cplex2D.get2dExteriorCellCount() == 8);
    REQUIRE(cplex2D.get2dInteriorCellCount() == 1);
    REQUIRE(cplex2D.get2dTotalCellCount() == 9); 

    //It's a rectangle and we have four edges on the boundary
    REQUIRE(cplex2D.get1dExteriorCellCount() == 24);
    REQUIRE(cplex2D.get1dInteriorCellCount() == 0);
    REQUIRE(cplex2D.get1dTotalCellCount() == 24);

    REQUIRE(cplex2D.get0dExteriorCellCount() == 16);
    REQUIRE(cplex2D.get0dInteriorCellCount() == 0);
    REQUIRE(cplex2D.get0dTotalCellCount() == 16);

    REQUIRE(cplex2D.getTotalCellCount() == 49);

} 

TEST_CASE("Any Cell Complex can have ghost cells", "[cplex-test]")
{
    Cajete::CartesianComplex2D cplex2D(4, 3, 6.0, 4.0, true);
    
    std::cout << cplex2D << std::endl;
   
    //Single rectangle means one 2D interior and 9 exterior boundary cells
    REQUIRE(cplex2D.get2dExteriorCellCount() == 18);
    REQUIRE(cplex2D.get2dInteriorCellCount() == 12);
    REQUIRE(cplex2D.get2dTotalCellCount() == 30); 

    //It's a rectangle and we have four edges on the boundary
    REQUIRE(cplex2D.get1dExteriorCellCount() == 54);
    REQUIRE(cplex2D.get1dInteriorCellCount() == 17);
    REQUIRE(cplex2D.get1dTotalCellCount() == 71);

    REQUIRE(cplex2D.get0dExteriorCellCount() == 36);
    REQUIRE(cplex2D.get0dInteriorCellCount() == 6);
    REQUIRE(cplex2D.get0dTotalCellCount() == 42);

    REQUIRE(cplex2D.getTotalCellCount() == 143);

} 
