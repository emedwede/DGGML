#include <iostream>

#include "catch.hpp"

#include "ExpandedComplex2D.hpp"
#include "CartesianHashFunctions.hpp"
#include "MemoryManager.hpp"
#include "Histobucket.hpp"
#include "PlantTypes.hpp"
#include "PlantUtils.hpp"
#include "YAGL_Graph.hpp"
#include <map>
#include <unordered_map>
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
/*
TEST_CASE("Histobucket 2D Test", "[binning-test]")
{    
    using key_type = typename Cajete::ExpandedComplex2D<>::graph_type::key_type;
    using node_type = typename YAGL::Graph<Cajete::Plant::mt_key_type, Cajete::Plant::MT_NodeData>::node_type;
    using cplex_type = Cajete::ExpandedComplex2D<>;
    
    Cajete::ExpandedComplex2D<> geoplex2D;
    geoplex2D.init(2, 2, 15.0, 15.0, false); //ghosted
    
    Cajete::Histobucket<key_type> geosysmap2D;
    std::cout << geosysmap2D << std::endl;

    REQUIRE( geosysmap2D.numBin() == 0 );

    geosysmap2D.reset_and_realloc(geoplex2D.get2dTotalCellCount());
    REQUIRE(geosysmap2D.numBin() == geoplex2D.get2dTotalCellCount()); 
    REQUIRE(geosysmap2D.totalSize() == 0);
    
    std::size_t ppmt = 3; std::size_t num_mt = 50;
    YAGL::Graph<Cajete::Plant::mt_key_type, Cajete::Plant::MT_NodeData> system_graph; 
    Cajete::Plant::microtubule_unit_scatter(system_graph, geoplex2D, num_mt); 
        
    auto& geoplex_graph = geoplex2D.getGraph();
    
    std::size_t counts2D[geoplex2D.get2dTotalCellCount()];
    std::size_t offsets2D[geoplex2D.get2dTotalCellCount()];
    std::size_t capacities2D[geoplex2D.get2dTotalCellCount()];
    
    std::map<key_type, std::size_t> sort_map2D;
    std::map<key_type, std::vector<key_type>> buckets;

    key_type keys2D[geoplex2D.get2dTotalCellCount()];
    key_type sort_map[geoplex2D.getTotalCellCount()]; 
    
    std::size_t c2 = 0; std::size_t c3 = 0;
    for(auto iter = geoplex_graph.node_list_begin(); iter != geoplex_graph.node_list_end(); iter++)
    {
        auto id = iter->first;
        auto itype = iter->second.getData().type;
        if(itype == 0) { 
            keys2D[c2] = id; c2++;
            sort_map2D.insert({id, c3++});
            buckets.insert({id, {}});
        }
    }
    
    auto search = sort_map2D.find(16);
    if(search != sort_map2D.end())
    {
        std::cout << "Key found, has value: " << search->second << "\n";
    } else { std::cout << "Key not found\n"; }
    for(auto item : sort_map2D)
    {
        std::cout << "{ " << item.first << ", " << item.second << " } ";
    } std::cout << "\n";

    for(auto i = 0; i < geoplex2D.get2dTotalCellCount(); i++)
    {
        sort_map[keys2D[i]] = i; counts2D[i] = 0;
    }
    
    REQUIRE( c2 == geoplex2D.get2dTotalCellCount() );
    
    print_keys(keys2D, c2, 2);
    
    for(auto iter = system_graph.node_list_begin(); iter != system_graph.node_list_end(); iter++)
    {
        auto id = iter->first;
        auto node = iter->second;
        
        auto cardinal = Cajete::cartesian_complex_expanded_hash2D(node, geoplex2D);
        auto search = sort_map2D.find(cardinal);
        if(search != sort_map2D.end())
            geosysmap2D.incrementBin(search->second);
        if(buckets.find(cardinal) != buckets.end())
            buckets.find(cardinal)->second.push_back(id);
        counts2D[sort_map[cardinal]]++;
    }
    
    auto bucket_tots = 0;
    for(auto item : buckets)
    {
        bucket_tots += item.second.size();
    }
    std::cout << "Bucket tots: " << bucket_tots << "\n";
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
    
    //geosysmap2D.build_buckets();


    //reset counts2D 
    for(auto i = 0; i < geosysmap2D.numBin(); i++) {counts2D[i] = 0;}
    
    std::size_t count0D = 0;
    std::size_t alt2D = 0;
    for(auto iter = system_graph.node_list_begin(); iter != system_graph.node_list_end(); iter++)
    {
        auto id = iter->first;
        auto node = iter->second; 

        auto cardinal = Cajete::cartesian_complex_expanded_hash2D(node, geoplex2D);
        Cajete::cartesian_complex_expanded_hash1D(node, geoplex2D);
        auto id0D = Cajete::cartesian_complex_expanded_hash0D(node, geoplex2D);
        if(geoplex_graph.findNode(id0D)->second.getData().type == 2)
        {
            count0D++;
        } else { alt2D++; }
        auto pid = sort_map[cardinal];
        
        auto loc = offsets2D[pid]+counts2D[pid];
        counts2D[pid]++;
        bucket_data[loc] = id;
    }
    REQUIRE(count0D+alt2D == num_mt*ppmt);
    std::cout << "OD objects: " << count0D << "\n";
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
    
    std::cout << "Counts: ";
    std::size_t count_sum = 0;
    for(auto i = 0; i < geoplex2D.get2dTotalCellCount(); i++) 
    {
        count_sum += counts2D[i];
        std::cout << counts2D[i] << " ";
    } std::cout << std::endl;
    REQUIRE(count_sum == num_mt*ppmt);
    REQUIRE(geosysmap2D.totalSize() == num_mt*ppmt);
}*/

TEST_CASE("VectorMap ND Test", "[binning-test]")
{    
    using key_type = typename Cajete::ExpandedComplex2D<>::graph_type::key_type;
    using cplex_type = Cajete::ExpandedComplex2D<>;
    
    cplex_type geoplex2D;
    geoplex2D.init(2, 2, 15.0, 15.0, true); //ghosted
    auto& geoplex_graph = geoplex2D.getGraph();
    
    std::size_t ppmt = 3; std::size_t num_mt = 5'000; std::size_t total_objects = ppmt*num_mt;
    YAGL::Graph<Cajete::Plant::mt_key_type, Cajete::Plant::MT_NodeData> system_graph; 
    std::cout << "Initializing system graph\n";
    Cajete::Plant::microtubule_unit_scatter(system_graph, geoplex2D, num_mt); 
    
    std::cout << "Sorting the graph\n";
    std::map<key_type, std::vector<key_type>> bucketsND[3];
    std::size_t complementND[3] = {0, 0, 0};

    Cajete::expanded_cartesian_complex_sort_stl(bucketsND, complementND, geoplex2D, system_graph);
    
    std::cout << "Computing bucket sums\n";

    std::size_t sumsND[3] = {0, 0, 0};
    
    for(int i = 0; i < 3; i++)
    {
        for(auto item : bucketsND[i])
        {
            sumsND[i] += item.second.size();
        }
    }
    for(auto i = 0; i < 3; i++)
    {   
        std::cout << "Dimension " << 2 - i << " will search a space of: " << sumsND[i] << " system objects\n";
        REQUIRE(sumsND[i] + complementND[i] == total_objects);
    }
}
