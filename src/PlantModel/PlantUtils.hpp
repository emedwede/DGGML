#ifndef PLANT_UTILS_HPP
#define PLANT_UTILS_HPP

#include <random>

#include "PlantTypes.hpp"
#include "MathUtils.hpp"

namespace Cajete
{
namespace Plant 
{
    template <typename GraphType, typename CplexType, typename ParamType>
    void microtubule_unit_scatter(GraphType& graph, CplexType& cplex, ParamType& settings)
    {
        double epsilon_min = settings.MT_MIN_SEGMENT_INIT;
        double epsilon_max = settings.MT_MAX_SEGMENT_INIT;
        std::random_device random_device;
        std::mt19937 random_engine(random_device());
        auto& grid = cplex.getCoarseGrid();
        
        std::uniform_real_distribution<double> 
            distribution_global_x(cplex.min_x+epsilon_max, 
                    cplex.max_x-epsilon_max);

        std::uniform_real_distribution<double>
            distribution_global_y(cplex.min_y+epsilon_max, 
                    cplex.max_y-epsilon_max);

        std::uniform_real_distribution<double> 
            distribution_local(epsilon_min, epsilon_max);
        
        std::uniform_real_distribution<double>
            distribution_angle(0.0, 2.0*3.14);
        using node_type = typename GraphType::node_type;

        std::size_t segments = 3;
        for(auto i = 0; i < settings.NUM_MT; i++) 
        {
            auto x_c = distribution_global_x(random_engine);
            auto y_c = distribution_global_y(random_engine);
            auto z_c = 0.0; 

            auto theta = distribution_angle(random_engine);
            auto seg_len = distribution_local(random_engine);
            std::cout << x_c << " " << y_c << "\n"; 
            std::cout << "Theta: " << theta << "\n";
            auto x_s = 0.0;
            auto y_s = seg_len;
            auto x_r_t = x_s*cos(theta) + y_s*sin(theta);
            auto y_r_t = -x_s*sin(theta) +y_s*cos(theta);
            auto x_r = x_c + x_r_t;
            auto y_r = y_c + y_r_t;
            auto z_r = 0.0;
            std::cout << x_r_t << " " << y_r_t << "\n";
            auto x_l = x_c - (x_r - x_c);
            auto y_l = y_c - (y_r - y_c);
            auto z_l = 0.0;
            std::cout << x_l << " " << y_l << "\n";
            //compute dist and unit vector
            double p1[3] = {x_r - x_c, y_r - y_c, 0.0};
            double p2[3] = {0.0, 0.0, 0.0};
            double p3[3] = {x_c - x_l, y_c - y_l, 0.0};
            double u1[3], u2[3];

            set_unit_vector(p1, p2, u1);
            set_unit_vector(p3, p2, u2);
            node_type node_l(i*segments, 
                    {{x_l, y_l, z_l}, 
                    {0.0, 0.0, 0.0}, 
                    negative, 
                    {-1, -1, -1}, 
                    {u2[0], u2[1], u2[2]}});

            node_type node_c(i*segments+1, 
                    {{x_c, y_c, z_c}, 
                    {0.0, 0.0, 0.0}, 
                    intermediate, 
                    {-1, -1, -1}, 
                    {u1[0], u1[1], u1[2]}});
            
            node_type node_r(i*segments+2, 
                    {{x_r, y_r, z_r}, 
                    {0.0, 0.0, 0.0}, 
                    positive, 
                    {-1, -1, -1}, 
                    {u1[0], u1[1], u1[2]}});

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
