#ifndef CAJETE_PLANT_TYPES_HPP
#define CAJETE_PLANT_TYPES_HPP

#include "YAGL_Node.hpp"

namespace Cajete 
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
            junction
        };

        struct MT_NodeData
        {
            double position[DIM3D];
            double velocity[DIM3D];
            int type; 
            int64_t tagND[3]; 
        };
    } // end namespace Plant
}// end namespace Cajete

#endif
