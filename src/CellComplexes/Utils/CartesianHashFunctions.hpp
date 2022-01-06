#ifndef CARTESIAN_HASH_FUNCTIONS_HPP
#define CARTESIAN_HASH_FUNCTIONS_HPP

#include <cstdint>

namespace Cajete
{

template <typename NodeType, typename CplexType>
int64_t cartesian_complex_expanded_hash2D(NodeType& node, CplexType& geoplex2D)
{
    auto& fine_grid = geoplex2D.getFineGrid();
    auto& coarse_grid = geoplex2D.getCoarseGrid();
    auto& geoplex_graph = geoplex2D.getGraph();
    auto&  node_data = node.getData();
        
    double xp = node_data.position[0];
    double yp = node_data.position[1];
    int ic, jc;

    coarse_grid.locatePoint(xp, yp, ic, jc); 
    if(ic > coarse_grid._nx || ic < 0) 
    {
        if(jc > coarse_grid._ny || jc < 0)
        {
            for(auto i = 0; i < 3; i++) node_data.tagND[i] = -1;
            //out of bounds 
            return -1;
        }
    }

    geoplex2D.coarse_cell_to_fine_lattice(ic, jc);
    
    auto cardinal = fine_grid.cardinalLatticeIndex(ic, jc);
    
    for(auto i = 0; i < 3; i++) node_data.tagND[i] = cardinal;

    return cardinal;
}

template <typename NodeType, typename CplexType>
int64_t cartesian_complex_expanded_hash1D(NodeType& node, CplexType& geoplex2D)
{
    auto& fine_grid = geoplex2D.getFineGrid();
    auto& coarse_grid = geoplex2D.getCoarseGrid();
    auto& geoplex_graph = geoplex2D.getGraph();
    auto&  node_data = node.getData();
        
    double xp = node_data.position[0];
    double yp = node_data.position[1];
    int ic, jc;

    coarse_grid.locatePoint(xp, yp, ic, jc);
    if(ic > coarse_grid._nx || ic < 0) 
    {
        if(jc > coarse_grid._ny || jc < 0)
        {
            for(auto i = 0; i < 3; i++) node_data.tagND[i] = -1;
            //out of bounds 
            return -1;
        }
    }
    geoplex2D.coarse_cell_to_fine_lattice(ic, jc);
    
    auto cardinal = fine_grid.cardinalLatticeIndex(ic, jc);
    
    for(auto jter = geoplex_graph.out_neighbors_begin(cardinal); jter != geoplex_graph.out_neighbors_end(cardinal); jter++)
    {
        auto jd = *jter;
        auto jnode = geoplex_graph.findNode(jd)->second;
        auto jnode_data = jnode.getData();
        if(jnode_data.type == 1 && jnode_data.interior) //check to see which 1D zone, if any the point belongs
        {
            auto& corners = jnode_data.corners;
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
                    return jd; //returns the expanded 1D cell is belongs to
                }
            }
        }
    }
    
    return cardinal; //returns the 2D cell it belongs to
}

template <typename NodeType, typename CplexType>
int64_t cartesian_complex_expanded_hash0D(NodeType& node, CplexType& geoplex2D)
{
    auto& fine_grid = geoplex2D.getFineGrid();
    auto& coarse_grid = geoplex2D.getCoarseGrid();
    auto& geoplex_graph = geoplex2D.getGraph();
    auto&  node_data = node.getData();
        
    double xp = node_data.position[0];
    double yp = node_data.position[1];
    int ic, jc;

    coarse_grid.locatePoint(xp, yp, ic, jc);
    if(ic >= coarse_grid._nx || ic < 0) 
    {
        if(jc >= coarse_grid._ny || jc < 0)
        {
            for(auto i = 0; i < 3; i++) node_data.tagND[i] = -1;
            //out of bounds 
            return -1;
        }
    }
    geoplex2D.coarse_cell_to_fine_lattice(ic, jc);
    
    auto cardinal = fine_grid.cardinalLatticeIndex(ic, jc);
    
    for(auto jter = geoplex_graph.out_neighbors_begin(cardinal); jter != geoplex_graph.out_neighbors_end(cardinal); jter++)
    {
        auto jd = *jter;
        auto jnode = geoplex_graph.findNode(jd)->second;
        auto jnode_data = jnode.getData();
        if(jnode_data.type == 1) //search the vertices connected to the edges
        {
            for(auto kter = geoplex_graph.out_neighbors_begin(jd); kter != geoplex_graph.out_neighbors_end(jd); kter++)
            {
                auto kd = *kter;
                auto knode = geoplex_graph.findNode(kd)->second;
                auto knode_data = knode.getData();
                if(knode_data.type == 2 && knode_data.interior) 
                {
                    auto& corners = knode_data.corners;
                    auto xmin = corners[0][0]; 
                    auto xmax = corners[0][0]; 
                    auto ymin = corners[0][1]; 
                    auto ymax = corners[0][1];
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
                            return kd; //returns the expanded 0D cell is belongs to
                        }
                    }
                }
            }
        }
    }
    
    return cardinal; //returns the 2D cell it belongs to
}

template <typename BucketType, typename ComplementType, typename GeoplexType, typename SysGraphType>
void expanded_cartesian_complex_sort_stl(BucketType& bucketsND, ComplementType& complementND, GeoplexType& geoplex2D, SysGraphType& system_graph)
{
    auto& geoplex_graph = geoplex2D.getGraph();
    
    for(auto iter = geoplex_graph.node_list_begin(); iter != geoplex_graph.node_list_end(); iter++)
    {
        auto id = iter->first;
        auto itype = iter->second.getData().type;
        if(itype == 0) { 
            bucketsND[0].insert({id, {}});
        }
        if(itype == 1)
        {
            bucketsND[1].insert({id, {}});
        }
        if(itype == 2)
        {
            bucketsND[2].insert({id, {}});
        }
    }
    
    for(auto iter = system_graph.node_list_begin(); iter != system_graph.node_list_end(); iter++)
    {
        auto id = iter->first;
        auto& node = iter->second;
        
        for(auto i = 0; i < 3; i++) 
        {
            std::size_t cardinal;
            if(i == 0) cardinal = cartesian_complex_expanded_hash2D(node, geoplex2D);
            if(i == 1) cardinal = cartesian_complex_expanded_hash1D(node, geoplex2D);
            if(i == 2) cardinal = cartesian_complex_expanded_hash0D(node, geoplex2D);

            auto search = bucketsND[i].find(cardinal);
            if(search != bucketsND[i].end())
                bucketsND[i].find(cardinal)->second.push_back(id);
            else
                complementND[i]++;
        }
    }
}

}

#endif
