#ifndef DGGML_CELL_COMPLEX_TYPES_HPP
#define DGGML_CELL_COMPLEX_TYPES_HPP

#include "YAGL_Graph.hpp"

namespace DGGML
{
    struct cell_complex_2D_node 
    {
        enum Fields 
        {
            Type, //Type of cell complex node
            Center, //Spatially embedded cneter location
            Corners, //The fattened up area of this node 
            Ghosted  //If the cell is ghosted out
        };
    };
    
    enum corner_types
    {
        lower_left,
        lower_right,
        upper_left,
        upper_right
    };

    struct cell_complex_2D_node_t
    {
        std::size_t type; // Type of the cell complex node
        double position[3]; // Spatially embedded center location
        double corners[4][2]; //The fattened up area of this node
        bool ghosted;
        bool interior;
    };
    
    using cell_complex_2D_key_t = std::size_t;

    using CplexGraph2D_t = YAGL::Graph<cell_complex_2D_key_t, cell_complex_2D_node_t>;    

}// end namespace DGGML

#endif
