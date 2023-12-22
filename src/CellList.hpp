#ifndef CELL_LIST_HPP
#define CELL_LIST_HPP

#include <vector>

#include "CartesianGrid2D.hpp"
#include "ComponentMatchMap.hpp"
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

    //TODO: maybe we need a multidimensional iterator, to iterate over the outer and inner
    //CellListData data;
    using ComponentMapType = ComponentMatchMap<typename GraphType::key_type>;
    using ComponentMapIter = typename ComponentMapType::iterator;
    using CellListObjectType = typename ComponentMapType::key_type;
    using CellListDataType = std::vector<std::vector<CellListObjectType>>;
    //using GraphType = YAGL::Graph<Plant::mt_key_type, Plant::MT_NodeData>;

    using outer_iterator = typename CellListDataType::iterator;
    using inner_iterator = typename std::vector<CellListObjectType>::iterator;

    struct IteratorPair
    {
        typename CellListDataType::iterator outer_iterator;
        typename std::vector<CellListObjectType>::iterator inner_iterator;
    };

    using iterator_pair = IteratorPair;

    CellListDataType data;

    CartesianGrid2D grid;
    GraphType* graph;
    ComponentMapType* component_matches;
    int cell_range;

    CellList() = default;

    CellList(CartesianGrid2D& _grid, GraphType* _graph, ComponentMapType* _component_matches)
        : grid(_grid), graph(_graph), component_matches(_component_matches), data(_grid.totalNumCells())
    {
        //TODO make the cell size ratio a paramater 
        double cell_ratio = 1;
        cell_range = std::ceil(1/cell_ratio);

        //build the cell list from the system graph 
        std::cout << "Building cell list\n";
        for(auto it = component_matches->begin(); it != component_matches->end(); it++) insert_match(it);
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

    std::size_t locate_cell(ComponentMapIter& comp_iter)
    {
        auto& match = comp_iter->second;
        auto& anchor_data = graph->findNode(match.anchor)->second.getData();
        double xp = anchor_data.position[0];
        double yp = anchor_data.position[1];
        auto key = comp_iter->first;
        int ic, jc;
        grid.locatePoint(xp, yp, ic, jc);
        return grid.cardinalCellIndex(ic, jc);
    }

    void insert_match(ComponentMapIter& comp_iter)
    {
        auto cardinal = locate_cell(comp_iter);
        data[cardinal].push_back(comp_iter->first);
    }


    //The find function must search the current cell and potentially it's neigbhors to see where a component is
    iterator_pair find(ComponentMapIter& comp_iter)
    {
        //first see if the component is in the cell it would be currently mapped to
        auto c = locate_cell(comp_iter);
        bool found = false;
        for(auto i = 0; i < data[c].size(); i++)
        {
            if(comp_iter->first == data[c][i]) {
                found = true;
                std::cout << "Found component " << comp_iter->first << " in cell " << c << " at index " << i << "\n";
                auto outIter = data.begin()+c;
                auto inIter = outIter->begin()+i;
                return iterator_pair{outIter, inIter};
            }
        }

        if(!found) //expand the search, since currently mapped cell, may not be the original (components move)
        {
            int imin, imax, jmin, jmax;
            getCells(c, imin, imax, jmin, jmax);
            for(auto i = imin; i < imax; i++)
            {
                for(auto j = jmin; j < jmax; j++)
                {
                    auto d = grid.cardinalCellIndex(i, j);
                    //if(d == c) continue;
                    for(auto k = 0; k < data[d].size(); k++)
                    {
                        if(comp_iter->first == data[c][k]) {
                            found = true;
                            std::cout << "Found component " << comp_iter->first << " in cell " << c << " at index " << k << "\n";
                            auto outIter = data.begin()+c;
                            auto inIter = outIter->begin()+i;
                            return iterator_pair{outIter, inIter};
                        }
                    }
                }
            }
        }
        //else we return that we didn't find anything
        auto outIter = data.end();
        auto inIter = outIter->end();
        return iterator_pair {outIter, inIter};
    }

    void erase(ComponentMapIter& comp_iter)
    {
        auto r = find(comp_iter);
        if(r.outer_iterator != data.end())
            r.outer_iterator->erase(r.inner_iterator);
    }

    void rebuild()
    {
        std::cout << "Rebuilding cell list\n";
        for(auto& cell : data)
            cell.clear();
        std::cout << "Cleared cell list mapped components: " << getTotalSize() << "\n";
        for(auto it = component_matches->begin(); it != component_matches->end(); it++) insert_match(it);
        std::cout << "Cell list has been built\n";
        std::cout << "The cell list has " << totalNumCells() << " reaction cells\n";
        std::cout << "Components mapped after reset: " << getTotalSize() << "\n";
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
