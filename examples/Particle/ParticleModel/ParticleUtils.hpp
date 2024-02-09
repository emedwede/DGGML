#ifndef PARTICLE_UTILS_HPP
#define PARTICLE_UTILS_HPP

#include <random>
#include <algorithm>
#include <vector>

#include "ParticleTypes.hpp"
#include "Utlities/MathUtils.hpp"

#include "CartesianGrid2D.hpp"

namespace MD
{
    template <typename GraphType, typename CplexType, typename ParamType, typename GenType>
    void particle_uniform_scatter(GraphType& graph, CplexType& cplex, ParamType& settings, GenType& gen) {
        std::random_device random_device;
        std::mt19937 random_engine(random_device());

        auto epsilon = settings.MAXIMAL_REACTION_RADIUS;
        auto num_particles = settings.NUM_PARTICLES;

        std::uniform_real_distribution<double> distribution_global_x(cplex.min_x + epsilon, cplex.max_x - epsilon);
        std::uniform_real_distribution<double> distribution_global_y(cplex.min_y + epsilon, cplex.max_y - epsilon);
        std::uniform_real_distribution<double> distribution_angle(0.0, 2.0 * 3.14);

        using molecule_type = typename graph_type::node_type;

        for (auto i = 0; i < num_particles; i++) {
            auto p_x = distribution_global_x(random_engine);
            auto p_y = distribution_global_x(random_engine);
            auto p_z = 0.0;

            //rotate terminal unit vector by theta
            auto theta = distribution_angle(random_engine);
            auto u_x = 1.0 * cos(theta) + 0.0 * sin(theta);
            auto u_y = -1.0 * sin(theta) + 0.0 * cos(theta);
            auto u_z = 0.0;

            graph_type g;
            molecule_type molecule = {
                    gen.get_key(),
                    {
                        Molecule{0.0,0.0,0.0, u_x, u_y, u_z},
                        p_x, p_y, p_z
                    }
            };

            graph.addNode(molecule);
        }
    }
}

#endif
