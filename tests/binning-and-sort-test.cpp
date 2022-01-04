#include <iostream>

#include "catch.hpp"

#include "ExpandedComplex2D.hpp"
#include "MemoryManager.hpp"
#include "Binner.hpp"
#include "PlantTypes.hpp"
#include "PlantUtils.hpp"
#include "YAGL_Graph.hpp"

template <typename T>
void print_corners(T& corners)
{
    std::cout << "{ " << corners[Cajete::lower_left][0] << ", " << corners[Cajete::lower_left][1] << "} ";
    std::cout << "{ " << corners[Cajete::lower_right][0] << ", " << corners[Cajete::lower_right][1] << "} ";
    std::cout << "{ " << corners[Cajete::upper_left][0] << ", " << corners[Cajete::upper_left][1] << "} ";
    std::cout << "{ " << corners[Cajete::upper_right][0] << ", " << corners[Cajete::upper_right][1] << "} ";
    std::cout << "}\n";
}

template <typename T>
void print_keys(T* keys, std::size_t n, std::size_t dim)
{
    std::cout << "Keys for dim " << dim << ": \n";
    for(auto i = 0; i < n; i++)
    {
        std::cout << keys[i] << " ";
    } std::cout << "\n\n";
}

TEST_CASE("Binner Init Test", "[binning-test]")
{    
    Cajete::Binner<int> bins;
    std::cout << bins << std::endl;

    REQUIRE( bins.numBin() == 0 );

    Cajete::ExpandedComplex2D<> geoplex2D;
    
    std::cout << "Generating the expanded cell complex\n";
    geoplex2D.init(1, 1, 15.0, 15.0, false); //ghosted
    std::cout << geoplex2D;
    
    bins.reset_and_realloc(geoplex2D.getTotalCellCount());

    REQUIRE(bins.numBin() == geoplex2D.getTotalCellCount());
    
    std::size_t bin_sum = 0;
    for(auto i = 0; i < bins.numBin(); i++)
    {
        bin_sum += bins.binSize(i);
    }

    REQUIRE(bin_sum == 0);

    YAGL::Graph<Cajete::Plant::mt_key_type, Cajete::Plant::MT_NodeData> system_graph; 
    std::cout << "Initializing the system graph\n";
    Cajete::Plant::microtubule_unit_scatter(system_graph, geoplex2D, 5); 
        
    auto& fine_grid = geoplex2D.getFineGrid();
    auto& coarse_grid = geoplex2D.getCoarseGrid();
    auto& geoplex_graph = geoplex2D.getGraph();
    
    std::size_t counts2D[geoplex2D.get2dTotalCellCount()];
    std::size_t counts1D[geoplex2D.get1dTotalCellCount()];
    std::size_t counts0D[geoplex2D.get0dTotalCellCount()];

    using key_type = typename Cajete::ExpandedComplex2D<>::graph_type::key_type;
    
    key_type keys2D[geoplex2D.get2dTotalCellCount()];
    key_type keys1D[geoplex2D.get1dTotalCellCount()];
    key_type keys0D[geoplex2D.get0dTotalCellCount()];
    
    std::size_t c2, c1, c0; c2 = c1 = c0 = 0;
    //creates a list of 2D, 1D, 0D cells
    for(auto iter = geoplex_graph.node_list_begin(); iter != geoplex_graph.node_list_end(); iter++)
    {
        auto id = iter->first;
        auto itype = iter->second.getData().type;
        if(itype == 0)
        {
            keys2D[c2] = id;
            c2++;
        }
        if(itype == 1)
        {
            keys1D[c1] = id;
            c1++;
        }
        if(itype == 2)
        {
            keys0D[c0] = id;
            c0++;
        }
    }

    REQUIRE( c2 == geoplex2D.get2dTotalCellCount() );
    REQUIRE( c1 == geoplex2D.get1dTotalCellCount() );
    REQUIRE( c0 == geoplex2D.get0dTotalCellCount() );
    
    print_keys(keys2D, c2, 2);
    print_keys(keys1D, c1, 1);
    print_keys(keys0D, c0, 0);

    //Test code for sorting, TODO: move into classes
    for(auto iter = system_graph.node_list_begin(); iter != system_graph.node_list_end(); iter++)
    {
        auto id = iter->first;
        auto inode_data = iter->second.getData();
        
        double xp = inode_data.position[0];
        double yp = inode_data.position[1];
        int ic, jc;

        coarse_grid.locatePoint(xp, yp, ic, jc);
        std::cout << "\nCoarse grid location: " << ic << " " << jc << "\n";
        geoplex2D.coarse_cell_to_fine_lattice(ic, jc);
        auto cardinal = fine_grid.cardinalLatticeIndex(ic, jc);
        bins.add_count(cardinal);
        
        for(auto jter = geoplex_graph.out_neighbors_begin(cardinal); jter != geoplex_graph.out_neighbors_end(cardinal); jter++)
        {
            auto jd = *jter;
            auto jnode = geoplex_graph.findNode(jd)->second;
            auto jnode_data = jnode.getData();
            if(jnode_data.type == 1 && jnode_data.interior) //check to see which 1D zone, if any the point belongs
            {
                std::cout << jd << " is an edge of 2D cell " << cardinal << "and it's interior flag " << jnode_data.interior << "\n";
                auto& corners = jnode_data.corners;
                print_corners(corners); std::cout << std::endl;
                auto xmin = corners[0][0]; auto xmax = corners[0][0]; auto ymin = corners[0][1]; auto ymax = corners[0][1];
                for(auto i = 1; i < 4; i++)
                {
                    if(corners[i][0] < xmin) xmin = corners[i][0];
                    if(corners[i][0] > xmax) xmax = corners[i][0];
                    if(corners[i][1] < ymin) ymin = corners[i][1];
                    if(corners[i][1] > ymax) ymax = corners[i][1];
                }
                if(xmin <= xp && xp <= xmax) 
                {
                    if(ymin <= yp && yp <= ymax)
                    {
                        std::cout << "\n------------------------------------------------------\n";
                        std::cout << "A point has been found in expanded lattice point: " << jd;
                        std::cout << "\n------------------------------------------------------\n";
                    }
                }
            }
        }
        std::cout << "Node " << id << " with " << "Points: { " << xp << ", " << yp << " } ";
        std::cout << "belong to 2D expanded lattice point " << cardinal;
        std::cout << " with corners: \n";
        auto corners = geoplex_graph.findNode(cardinal)->second.getData().corners;
        print_corners(corners);
    }

    bin_sum = 0;
    for(auto i = 0; i < bins.numBin(); i++)
    {
        bin_sum += bins.binSize(i);
    }

    REQUIRE(bin_sum == 15);
}
