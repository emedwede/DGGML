#ifndef CELL_LIST_HPP
#define CELL_LIST_HPP

#include <vector>

#include "CartesianGrid2D.hpp"
#include "ComponentMap.hpp"
#include "PlantTypes.hpp"
#include "YAGL_Graph.hpp"

namespace DGGML
{

//type for holding the CellList Data in CSR format
//could be customized for various key types
using ViewType = std::vector<std::size_t>;

//This is the CSR layout version
struct CellListDataCSR
{
    //number of objects per cell
    ViewType counts; 

    //offsets for cells 
    ViewType offsets;

    //keys associated with a cell, indexed by offsets
    ViewType permutation;

};

//CellList Data in a noncontiguous binned format
struct CellListDataVec
{

};

template <typename GraphType>
struct CellList
{

    //CellListData data;
    using ComponentMapType = ComponentMap<typename GraphType::key_type>;
    using ComponentMapIter = typename ComponentMapType::iterator;
    using CellListObjectType = typename ComponentMapType::key_type;
    using CellListDataType = std::vector<std::vector<CellListObjectType>>;
    //using GraphType = YAGL::Graph<Plant::mt_key_type, Plant::MT_NodeData>;
    CellListDataType data;

    CartesianGrid2D grid;
    GraphType& graph;
    ComponentMapType& component_matches;
    int cell_range;

    CellList(CartesianGrid2D& _grid, GraphType& _graph, ComponentMapType& _component_matches)
        : grid(_grid), graph(_graph), component_matches(_component_matches), data(_grid.totalNumCells())
    {
        //TODO make the cell size ratio a paramater 
        double cell_ratio = 1;
        cell_range = std::ceil(1/cell_ratio);

        //build the cell list from the system graph 
        std::cout << "Building cell list\n";
        for(auto& item : component_matches) insert_match(item);
        std::cout << "Cell list has been built\n";
        std::cout << "The cell list has " << totalNumCells() << " reaction cells\n"; 
        std::cout << "Components mapped: " << getTotalSize() << "\n";
    }
    
    std::size_t getCellSize(const int cell)
    {
        return data[cell].size();
    }
    
    std::size_t getTotalSize()
    {
        auto sum = 0;
        for(const auto& cell : data) sum += cell.size();
        return sum;
    }
    template<typename T> //Temporary fix since ComponentMapIter throws an error
    std::size_t locate_cell(T& comp_iter)
    {
        auto& match = comp_iter.second;
        auto& anchor_data = graph.findNode(match.anchor)->second.getData();
        double xp = anchor_data.position[0];
        double yp = anchor_data.position[1];
        auto key = comp_iter.first;
        int ic, jc;
        grid.locatePoint(xp, yp, ic, jc);
        return grid.cardinalCellIndex(ic, jc);
    }

    template<typename T> //Temporary fix since ComponentMapIter throws an error
    void insert_match(T& comp_iter)
    {
        auto cardinal = locate_cell(comp_iter);
        data[cardinal].push_back(comp_iter.first);
    }
    
    std::size_t totalNumCells() { return grid.totalNumCells(); }
    
    std::size_t cardinalCellIndex(int i, int j) { return grid.cardinalCellIndex(i, j); }

    // give a cell, get it's stencil
    void getCells(const int cell, int& imin, int& imax, int& jmin, int& jmax) 
    {
        int i, j;
        grid.ijCellIndex(cell, i, j);
        //TODO: there may be issues because this formulation does not ignore 
        //ghosted reaction cells
        imin = ( i - cell_range > 0 ) ? i - cell_range : 0;
        imax = ( i + cell_range + 1 ) < grid._nx ? i + cell_range + 1 : grid._nx;
        
        jmin = ( j - cell_range > 0 ) ? j - cell_range : 0;
        jmax = ( j + cell_range + 1 ) < grid._ny ? j + cell_range + 1 : grid._ny;
    }
};

}
#endif 
