#ifndef DGGML_SPATIALDATAE3D_HPP
#define DGGML_SPATIALDATA3D_HPP

#include "YAGL_Variant_Data.hpp"

template <typename ... Ts>
struct SpatialNode3D : YAGL::VariantData<Ts ...>
{
    double position[3];
};

#endif //DGGML_SPATIALDATA3D_HPP
