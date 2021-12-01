#ifndef CAJETE_CELL_COMPLEX_TYPES_HPP
#define CAJETE_CELL_COMPLEX_TYPES_HPP

#include "YAGL_Graph.hpp"

namespace Cajete 
{
    struct cell_complex_2D_node 
    {
        enum Fields 
        {
            Type, //Type of cell complex node
            Center, //Spatially embedded cneter location
            Corners //The fattened up area of this node 
        };
    };

    struct cell_complex_2D_node_t
    {
        std::size_t type; // Type of the cell complex node
        double center[2]; // Spatially embedded center location
        double corners[4][2]; //The fattened up area of this node
    };
    
    using cell_complex_2D_key_t = std::size_t;

    using CplexGraph2D_t = YAGL::Graph<cell_complex_2D_key_t, cell_complex_2D_node_t>;    

}// end namespace Cajete

#endif
