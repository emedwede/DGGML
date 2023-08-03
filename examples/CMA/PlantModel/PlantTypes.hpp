#ifndef DGGML_PLANT_TYPES_HPP
#define DGGML_PLANT_TYPES_HPP

#include "YAGL_Graph.hpp"
#include "YAGL_Node.hpp"

namespace DGGML
{
    namespace Plant
    {
        const uint8_t DIM3D = 3;

        using mt_key_type = std::size_t;

        enum mt_type 
        {
            negative,
            intermediate,
            positive,
            junction,
            zipper
        };

        struct MT_NodeData
        {
            double position[DIM3D];
            double velocity[DIM3D];
            int type; 
            int64_t tagND[3]; 
            double unit_vec[DIM3D];
        };

        using key_type = DGGML::Plant::mt_key_type;
        using data_type = DGGML::Plant::MT_NodeData;

        using graph_type = YAGL::Graph<key_type, data_type>;
        using node_type = YAGL::Node<key_type, data_type>;
    } // end namespace Plant
}// end namespace DGGML

#endif
