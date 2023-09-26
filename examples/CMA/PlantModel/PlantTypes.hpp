#ifndef DGGML_PLANT_TYPES_HPP
#define DGGML_PLANT_TYPES_HPP

#include "YAGL_Graph.hpp"
#include "YAGL_Node.hpp"
#include "SpatialData3D.hpp"

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

        struct Negative
        {
            double velocity[DIM3D];
            double unit_vec[DIM3D];
        };

        struct Intermediate
        {
            double velocity[DIM3D];
            double unit_vec[DIM3D];
        };

        struct Positive
        {
            double velocity[DIM3D];
            double unit_vec[DIM3D];
        };

        struct Junction
        {
            double unit_vec[DIM3D];
        };

        struct Zipper
        {
            double unit_vec[DIM3D];
        };

        using key_type = DGGML::Plant::mt_key_type;
        using data_type = DGGML::Plant::MT_NodeData;

        using graph_type = YAGL::Graph<key_type, data_type>;
        using node_type = YAGL::Node<key_type, data_type>;

        using spatial_graph_type = YAGL::Graph<key_type, SpatialNode3D<Negative, Intermediate, Positive, Junction, Zipper>>;
    } // end namespace Plant
}// end namespace DGGML

namespace Plant
{

    const uint8_t DIM3D = 3;

    struct Negative
    {
        double velocity[DIM3D];
        double unit_vec[DIM3D];
    };

    struct Intermediate
    {
        double velocity[DIM3D];
        double unit_vec[DIM3D];
    };

    struct Positive
    {
        double velocity[DIM3D];
        double unit_vec[DIM3D];
    };

    struct Junction
    {
        double unit_vec[DIM3D];
    };

    struct Zipper
    {
        double unit_vec[DIM3D];
    };

    using key_type = std::size_t;

    using graph_type = YAGL::Graph<key_type, SpatialNode3D<Negative, Intermediate, Positive, Junction, Zipper>>;
}

#endif