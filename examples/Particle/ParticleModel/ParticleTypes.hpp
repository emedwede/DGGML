#ifndef DGGML_PARTICLE_TYPES_HPP
#define DGGML_PARTICLE_TYPES_HPP

#include "YAGL_Graph.hpp"
#include "YAGL_Node.hpp"
#include "SpatialData3D.hpp"

namespace MD
{

    const uint8_t DIM3D = 3;

    struct Molecule
    {
        double velocity[DIM3D];
        double unit_vec[DIM3D];
    };

    using key_type = std::size_t;

    using graph_type = YAGL::Graph<key_type, SpatialNode3D<Molecule>>;
}

#endif
