#ifndef PLANT_UTILS_HPP
#define PLANT_UTILS_HPP

#include <random>

#include "PlantTypes.hpp"

namespace Cajete
{
namespace Plant 
{
    template <typename GraphType>
    void microtubule_unit_scatter(GraphType& graph)
    {
        std::random_device random_device;
        std::mt19937 random_engine(random_device());
        std::uniform_real_distribution<double> distribution_global(1.0, 9.0);
        std::uniform_real_distribution<double> distribution_local(0.5, 1.0);
        
        using node_type = typename GraphType::node_type;

        std::size_t num_mt = 32;
        std::size_t segments = 3;
        for(auto i = 0; i < num_mt; i++) 
        {
            auto x_c = distribution_global(random_engine);
            auto y_c = distribution_global(random_engine);
            auto z_c = 0.0; //distribution_global(random_engine);

            auto x_r = x_c + distribution_local(random_engine);
            auto y_r = y_c + distribution_local(random_engine); 
            auto z_r = 0.0; //z_c + distribution_local(random_engine);
            
            auto x_l = x_c - distribution_local(random_engine);
            auto y_l = y_c - distribution_local(random_engine); 
            auto z_l = 0.0; //z_c - distribution_local(random_engine);
            
            node_type node_l(i*segments, {{x_l, y_l, z_l}, {0.0, 0.0, 0.0}, negative});
            node_type node_c(i*segments+1, {{x_c, y_c, z_c}, {0.0, 0.0, 0.0}, intermediate});
            node_type node_r(i*segments+2, {{x_r, y_r, z_r}, {0.0, 0.0, 0.0}, positive});

            graph.addNode(node_l);
            graph.addNode(node_c);
            graph.addNode(node_r);

            graph.addEdge(node_l, node_c);
            graph.addEdge(node_r, node_c);
        }

    }
} //end namespace Plant
} //end namespace Cajete

#endif
