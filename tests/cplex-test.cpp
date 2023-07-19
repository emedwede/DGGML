#include <iostream>

#include "catch.hpp"

#include "CartesianComplex2D.hpp"
#include "ExpandedComplex2D.hpp"

#include "VtkWriter.hpp"

#include <iostream>

template <typename CplexType>
std::size_t interior_integrity_check(CplexType& cplex2D)
{
   auto g = cplex2D.getGraph();

   std::size_t interior_nodes = 0;
   for(auto iter = g.node_list_begin(); iter != g.node_list_end(); iter++)
   {
       auto node_data = iter->second.getData();
       
       if(node_data.interior == true) interior_nodes++;
   }
   return interior_nodes;
}

template <typename CplexType>
std::size_t exterior_integrity_check(CplexType& cplex2D)
{
   auto g = cplex2D.getGraph();

   std::size_t exterior_nodes = 0;
   for(auto iter = g.node_list_begin(); iter != g.node_list_end(); iter++)
   {
       auto node_data = iter->second.getData();

       if(node_data.interior == false) exterior_nodes++;
   }
   return exterior_nodes;
}


template<typename CplexType>
std::size_t check_corner_connectivity(CplexType& cplex2D)
{
    //Test the corners 
    auto g = cplex2D.getGraph();

    int completed_corners = 0;
    for(auto iter = g.node_list_begin(); iter != g.node_list_end(); iter++) 
    {
        auto id = iter->first;
        auto node = iter->second;
        auto corners = node.getData().corners;
        int found_corners = 0;
        if(node.getData().type == 0) //check the 0D nbrs
        {
            for(auto jter = g.out_neighbors_begin(id); jter != g.out_neighbors_end(id); jter++)
            {
                auto jnode = g.findNode(*jter);
                if(jnode->second.getData().type == 1)
                {
                    auto jd = jnode->first;
                    for(auto kter = g.out_neighbors_begin(jd); kter != g.out_neighbors_end(jd); kter++)
                    {
                        auto knode = g.findNode(*kter);
                        if(knode->second.getData().type == 2)
                        {
                            for(auto i = 0; i < 4; i++)
                            {
                                if(knode->second.getData().position[0] == corners[i][0] && knode->second.getData().position[1] == corners[i][1])
                                {
                                    found_corners++;
                                }
                            }
                        }
                    }
                }
            }
        }
        //should be eight since each 0D vertex gets visited twice!
        if(found_corners == 8) completed_corners++;
    }

    return completed_corners;
}
TEST_CASE("Cell Complex from Cartesian 2D Grid Test", "[cplex-test]")
{    
    std::cout << "Running the Cell Complex from Cartesian 2D Grid Test\n";
    
    std::size_t n = 4, m = 3; //number of cells in the x and y direction
    double dx = 2.0, dy = 2.0; //size of each grid cell 

    DGGML::CartesianComplex2D cplex2D(n, m, dx, dy);

    std::cout << cplex2D << std::endl;
    
    YAGL::Node<int, double> my_node({2, 2.6});
    
    std::cout << "Data: " << my_node.getData() << "\n";
    auto my_data = my_node.getData();
    my_data = 3.6;
    std::cout << "Data: " << my_node.getData() << "\n";

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
    
    
    REQUIRE(check_corner_connectivity(cplex2D) == cplex2D.get2dTotalCellCount());
    
    REQUIRE(interior_integrity_check(cplex2D) == cplex2D.getTotalInteriorCellCount());
    REQUIRE(exterior_integrity_check(cplex2D) == cplex2D.getTotalExteriorCellCount());

    //Save the cell complex graph
    //DGGML::VtkFileWriter<typename DGGML::CartesianComplex2D<>::graph_type> writer;
    
    //writer.save(cplex2D.getGraph(), "cplex_graph");
}

TEST_CASE("Cell Complex can just be a single domain", "[cplex-test]")
{
    DGGML::CartesianComplex2D cplex2D(1, 1, 6.0, 4.0);
    
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
    REQUIRE(check_corner_connectivity(cplex2D) == cplex2D.get2dTotalCellCount());
}

TEST_CASE("Cell Complex can just be a single domain and have ghost cells", "[cplex-test]")
{
    DGGML::CartesianComplex2D cplex2D(1, 1, 6.0, 4.0, true);
    
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
    REQUIRE(check_corner_connectivity(cplex2D) == cplex2D.get2dTotalCellCount());
} 

TEST_CASE("Any Cell Complex can have ghost cells", "[cplex-test]")
{
    DGGML::CartesianComplex2D cplex2D(4, 3, 6.0, 4.0, true);
    
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
    REQUIRE(check_corner_connectivity(cplex2D) == cplex2D.get2dTotalCellCount());
}

TEST_CASE("Testing the expanded cell complex", "[cplex-test]") 
{
    DGGML::ExpandedComplex2D<> cplex2D(2, 2, 3.0, 2.0, true);
    
    std::cout << cplex2D << std::endl;
    
    std::cout << "Total number of cells: " << cplex2D.getTotalCellCount() << "\n";
    //REQUIRE(cplex2D.getTotalCellCount() == 9);
    
    auto& g = cplex2D.getGraph();
    for(auto iter = g.node_list_begin(); iter != g.node_list_end(); iter++)
    {
        auto data = iter->second.getData();

        std::cout << "\n{Type: " << data.type << ", Ghosted: " << data.ghosted << ", Interior: " << data.interior << "} \t";
        std::cout << "Corners: { ";
        
        for(auto i = 0; i < 4; i++) 
        {
                std::cout << "{ " << data.corners[i][0]
                << ", " << iter->second.getData().corners[i][1] << "} ";
        } std::cout << "}\n";
    }
    std::cout << "Reaction Grid: \n" << cplex2D.reaction_grid;
    DGGML::GridFileWriter writer;
    writer.save({cplex2D.reaction_grid, cplex2D.dim_label}, "labeled_grid_viz");
    DGGML::VtkFileWriter<typename DGGML::ExpandedComplex2D<>::types::graph_type> complex_writer;
    complex_writer.save(cplex2D.getGraph(), "cplex_viz");
}
