#ifndef PLANT_UTILS_HPP
#define PLANT_UTILS_HPP

#include <random>

#include "PlantTypes.hpp"
#include "MathUtils.hpp"

namespace Cajete
{
namespace Plant 
{
    template <typename GraphType, typename CplexType>
    void microtubule_unit_scatter(GraphType& graph, CplexType& cplex, std::size_t num_mt = 128)
    {
        double epsilon = 1.0;
        std::random_device random_device;
        std::mt19937 random_engine(random_device());
        auto grid = cplex.getCoarseGrid();
        
        std::cout << "Min/Max: " << cplex.min_x << " " << cplex.min_y << " " << cplex.max_x << " " << cplex.max_y << "\n";

        std::uniform_real_distribution<double> 
            distribution_global(cplex.min_x+epsilon, cplex.max_x-epsilon);
        std::uniform_real_distribution<double> 
            distribution_local(epsilon/2.0, epsilon);
        
        using node_type = typename GraphType::node_type;

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
            
            //compute dist and unit vector
            double p1[3] = {x_r - x_c, y_r - y_c, 0.0};
            double p2[3] = {0.0, 0.0, 0.0};
            auto len = calculate_distance(p1, p2);
            double u[3] = {p1[0]/len, p1[1]/len, 0.0};
            
            double p3[3] = {x_c - x_l, y_c - y_l, 0.0};
            len = calculate_distance(p3, p2);
            double u2[3] = {p3[0]/len, p3[1]/len, 0.0};

            node_type node_l(i*segments, 
                    {{x_l, y_l, z_l}, 
                    {0.0, 0.0, 0.0}, 
                    negative, 
                    {-1, -1, -1}, 
                    {u2[0], u2[1], u2[3]}});

            node_type node_c(i*segments+1, 
                    {{x_c, y_c, z_c}, 
                    {0.0, 0.0, 0.0}, 
                    intermediate, 
                    {-1, -1, -1}, 
                    {u[0], u[1], u[3]}});
            
            node_type node_r(i*segments+2, 
                    {{x_r, y_r, z_r}, 
                    {0.0, 0.0, 0.0}, 
                    positive, 
                    {-1, -1, -1}, 
                    {u[0], u[1], u[3]}});

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
