#ifndef CELL_LIST_HPP
#define CELL_LIST_HPP

#include <vector>

namespace Cajete
{

//type for holding the CellList Data in CSR format
//could be customized for various key types
using ViewType = std::vector<std::size_t>;

//This is the CSR layout version
struct CellListData
{
    //number of neighbors per object
    ViewType counts; 

    //offsets for the neighbor list 
    ViewType offsets;

    //neighbor list 
    ViewType neighbors;

    //TODO: maybe an addNeighbor function
};

struct CellList
{

    CellListData data;

    CellList() {}

};

}
#endif 
