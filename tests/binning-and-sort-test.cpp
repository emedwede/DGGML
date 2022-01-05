#include <iostream>

#include "catch.hpp"

#include "ExpandedComplex2D.hpp"
#include "MemoryManager.hpp"
#include "Histobucket.hpp"
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

/*TEST_CASE("Histobucket Init Test", "[binning-test]")
{    
    Cajete::Histobucket<int> geosysmap2D;
    std::cout << geosysmap2D << std::endl;

    REQUIRE( geosysmap2D.numBin() == 0 );

    Cajete::ExpandedComplex2D<> geoplex2D;
    
    std::cout << "Generating the expanded cell complex\n";
    geoplex2D.init(1, 1, 15.0, 15.0, false); //ghosted
    std::cout << geoplex2D;
    
    geosysmap2D.reset_and_realloc(geoplex2D.get2dTotalCellCount());

    REQUIRE(geosysmap2D.numBin() == geoplex2D.get2dTotalCellCount());
    
    std::size_t bin_sum = 0;
    for(auto i = 0; i < geosysmap2D.numBin(); i++)
    {
        bin_sum += geosysmap2D.binSize(i);
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
    
    std::size_t offsets2D[geoplex2D.get2dTotalCellCount()];
    std::size_t offsets1D[geoplex2D.get1dTotalCellCount()];
    std::size_t offsets0D[geoplex2D.get0dTotalCellCount()];

    std::size_t capacities2D[geoplex2D.get2dTotalCellCount()];
    std::size_t capacities1D[geoplex2D.get1dTotalCellCount()];
    std::size_t capacities0D[geoplex2D.get0dTotalCellCount()];

    using key_type = typename Cajete::ExpandedComplex2D<>::graph_type::key_type;
     
    key_type keys2D[geoplex2D.get2dTotalCellCount()];
    key_type keys1D[geoplex2D.get1dTotalCellCount()];
    key_type keys0D[geoplex2D.get0dTotalCellCount()];
    
    key_type sort_map[geoplex2D.getTotalCellCount()]; 
    
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
    
    for(auto i = 0; i < geoplex2D.get2dTotalCellCount(); i++)
    {
        sort_map[keys2D[i]] = i;
        counts2D[i] = 0;
    }
    
    for(auto i = 0; i < geoplex2D.get1dTotalCellCount(); i++)
    {
        sort_map[keys1D[i]] = i;
        counts1D[i] = 0;
    }
    
    for(auto i = 0; i < geoplex2D.get0dTotalCellCount(); i++)
    {
        sort_map[keys0D[i]] = i;
        counts0D[i] = 0;
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
        geosysmap2D.add_count(sort_map[cardinal]);
        counts2D[sort_map[cardinal]]++;
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
    for(auto i = 0; i < geosysmap2D.numBin(); i++)
    {
        bin_sum += geosysmap2D.binSize(i);
    }
    std::size_t count_sum = 0;
    for(auto i = 0; i < geoplex2D.get2dTotalCellCount(); i++) 
    {
        count_sum += counts2D[i];
    }
    REQUIRE(count_sum == 15);
    REQUIRE(bin_sum == 15);
}*/

template <typename T>
std::size_t bin_summation(const T& hist)
{
    std::size_t bin_sum = 0;
    for(auto i = 0; i < hist.numBin(); i++)
    {
        bin_sum += hist.binSize(i);
    }
    return bin_sum;
}

TEST_CASE("Histobucket 2D Test", "[binning-test]")
{    
    using key_type = typename Cajete::ExpandedComplex2D<>::graph_type::key_type;
     
    Cajete::Histobucket<key_type> geosysmap2D;
    std::cout << geosysmap2D << std::endl;

    REQUIRE( geosysmap2D.numBin() == 0 );

    Cajete::ExpandedComplex2D<> geoplex2D;
    geoplex2D.init(2, 2, 15.0, 15.0, true); //ghosted
    
    geosysmap2D.reset_and_realloc(geoplex2D.get2dTotalCellCount());
    REQUIRE(geosysmap2D.numBin() == geoplex2D.get2dTotalCellCount()); 
    REQUIRE(bin_summation(geosysmap2D) == 0);

    YAGL::Graph<Cajete::Plant::mt_key_type, Cajete::Plant::MT_NodeData> system_graph; 
    Cajete::Plant::microtubule_unit_scatter(system_graph, geoplex2D, 30); 
        
    auto& fine_grid = geoplex2D.getFineGrid();
    auto& coarse_grid = geoplex2D.getCoarseGrid();
    auto& geoplex_graph = geoplex2D.getGraph();
    
    std::size_t counts2D[geoplex2D.get2dTotalCellCount()];
    std::size_t offsets2D[geoplex2D.get2dTotalCellCount()];
    std::size_t capacities2D[geoplex2D.get2dTotalCellCount()];

    key_type keys2D[geoplex2D.get2dTotalCellCount()];
    key_type sort_map[geoplex2D.getTotalCellCount()]; 
    
    std::size_t c2 = 0;
    for(auto iter = geoplex_graph.node_list_begin(); iter != geoplex_graph.node_list_end(); iter++)
    {
        auto id = iter->first;
        auto itype = iter->second.getData().type;
        if(itype == 0) { keys2D[c2] = id; c2++; }
    }
    
    for(auto i = 0; i < geoplex2D.get2dTotalCellCount(); i++)
    {
        sort_map[keys2D[i]] = i; counts2D[i] = 0;
    }
    
    REQUIRE( c2 == geoplex2D.get2dTotalCellCount() );
    
    print_keys(keys2D, c2, 2);

    for(auto iter = system_graph.node_list_begin(); iter != system_graph.node_list_end(); iter++)
    {
        auto id = iter->first;
        auto inode_data = iter->second.getData();
        
        double xp = inode_data.position[0];
        double yp = inode_data.position[1];
        int ic, jc;

        coarse_grid.locatePoint(xp, yp, ic, jc);
        
        geoplex2D.coarse_cell_to_fine_lattice(ic, jc);
        
        auto cardinal = fine_grid.cardinalLatticeIndex(ic, jc);
        
        geosysmap2D.add_count(sort_map[cardinal]);
        
        counts2D[sort_map[cardinal]]++;
    }
    
    //find the max size
    std::size_t max_size = geosysmap2D.maxSize();

    //set the 2D capcities and binned offsets 
    std::size_t reserve_factor = 3;
    std::size_t max_cap = max_size*reserve_factor;
    for(auto i = 0; i < geosysmap2D.numBin(); i++)
    {
        geosysmap2D.setBinCapacity(max_cap, i);
        geosysmap2D.setBinOffset(i*max_cap, i);
        capacities2D[i] = max_cap;
        offsets2D[i] = i*max_cap;
    }
    key_type bucket_data[geosysmap2D.numBin()*max_cap];
   
    //reset counts2D 
    for(auto i = 0; i < geosysmap2D.numBin(); i++) {counts2D[i] = 0;}

    for(auto iter = system_graph.node_list_begin(); iter != system_graph.node_list_end(); iter++)
    {
        auto id = iter->first;
        auto inode_data = iter->second.getData();
        
        double xp = inode_data.position[0];
        double yp = inode_data.position[1];
        int ic, jc;

        coarse_grid.locatePoint(xp, yp, ic, jc);
        
        geoplex2D.coarse_cell_to_fine_lattice(ic, jc);
        
        auto cardinal = fine_grid.cardinalLatticeIndex(ic, jc);
        
        auto pid = sort_map[cardinal];
        auto loc = offsets2D[pid]+counts2D[pid];
        counts2D[pid]++;
        bucket_data[loc] = id;
    }
    
    for(auto i = 0; i < geosysmap2D.numBin(); i++) {
        std::cout << "Nodes in bin " << keys2D[i] << ": { ";
        if(counts2D[i] > 0)
        {
            auto start = offsets2D[i];
            auto end = offsets2D[i] + counts2D[i];
            for(auto j = start; j < end; j++) 
            {
                std::cout << bucket_data[j] << " ";
            } 
        } std::cout << "}\n";
    }
    
    for(auto i = 0; i < c2; i++)
    {
        auto key = keys2D[i];
        auto ghosted = geoplex2D.getGraph().findNode(key)->second.getData().ghosted;
        if(ghosted) 
        {
            std::cout << "We do not simulate cell " << key << ", since it is ghosted\n";
        }
        else 
        {
            std::cout << "Cell " << key << " is not ghosted, so we simulated it!\n";

            auto pid = sort_map[key];
            auto start = offsets2D[pid];
            auto end = offsets2D[pid] + counts2D[pid];
            
            for(auto j = start; j < end; j++) {
                std::cout << "Performing a pattern search on object " << bucket_data[j] << "\n";
            }

        }
    }
    std::size_t count_sum = 0;
    for(auto i = 0; i < geoplex2D.get2dTotalCellCount(); i++) 
    {
        count_sum += counts2D[i];
    }
    REQUIRE(count_sum == 15);
    REQUIRE(bin_summation(geosysmap2D) == 15);
}
