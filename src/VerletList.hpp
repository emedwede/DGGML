#ifndef VERLET_LIST_HPP
#define VERLET_LIST_HPP

#include <vector>

namespace Cajete
{

//type for holding the VerletList Data in CSR format
//could be customized for various key types
using ViewType = std::vector<std::size_t>;

//This is the CSR layout version
struct VerletListData
{
    //number of neighbors per object
    ViewType counts; 

    //offsets for the neighbor list 
    ViewType offsets;

    //neighbor list 
    ViewType neighbors;

    //TODO: maybe an addNeighbor function
};

struct VerletList
{

    VerletListData data;

    VerletList() {}

};

}
#endif 
