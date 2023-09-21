#ifndef DGGML_APPROXIMATESSA_HPP
#define DGGML_APPROXIMATESSA_HPP

#include "PlantTypes.hpp"
#include "PlantGrammar.hpp"

#include "ExpandedComplex2D.hpp"
#include "YAGL_Graph.hpp"
#include "YAGL_Algorithms.hpp"

#include "Solver.h"

#include "RuleSystem.hpp"

#include <chrono>
#include <random>
#include <algorithm>
#include <vector>
#include <iostream>
#include <math.h>
#include <map>

namespace DGGML
{

//TODO: Move to it's own random class or header
auto RandomRealsBetween = [](double low, double high)
{
    auto randomFunc =
            [distribution_ = std::uniform_real_distribution<double>(low, high),
                    random_engine_ = std::mt19937{std::random_device{}() }]() mutable
            {
                return distribution_(random_engine_);
            };
    return randomFunc;
};

auto RandomIntsBetween = [](int low, int high)
{
    auto randomFunc =
            [distribution_ = std::uniform_int_distribution<int>(low, high),
                    random_engine_ = std::mt19937{std::random_device{}() }]() mutable
            {
                return distribution_(random_engine_);
            };
    return randomFunc;
};

template <typename T1, typename T2, typename T3, typename M1>
void approximate_ssa(T1& rule_system, T2& rule_map, T3& rule_instances, M1& model,
                                          std::size_t k, std::pair<double, double>& geocell_progress)
{
    auto& geoplex2D = model->geoplex2D;
    auto& system_graph = model->system_graph;
    auto& settings = model->settings;
    auto& gamma = model->gamma;

    std::cout << "Cell " << k << " has " << rule_map[k].size() << " rules\n";

    double delta_t, geocell_propensity;
    std::size_t events, steps;
    delta_t = 0.0; events = 0; steps = 0, geocell_propensity = 0.0;
    double tau = geocell_progress.first;
    double exp_sample = geocell_progress.second;

    if(exp_sample <= 0.0)
    {
        double uniform_sample = RandomRealsBetween(0.0, 1.0)();
        exp_sample = -log(1-uniform_sample);
        std::cout << "Warped wating time: " << exp_sample << "\n";
    }

    while(delta_t < settings.DELTA)
    {
        //the propensities calculated for a particular rule
        std::vector<std::pair<std::size_t, double>> propensity_space;

        if(tau < exp_sample)
        {
            // STEP(1) : sum all of the rule propensities

            for(const auto& r : rule_map[k])
            {
                //TODO: need to replace with key type, and need more efficient method for turning
                // keys into compatible inputs into propensity functions
                std::vector<std::size_t> inducers;
                auto& inst = rule_instances[r];
                if(inst.category != "stochastic")
                    continue;
                for(auto& c : inst.components)
                {
                    for(auto& id : rule_system[c])
                        inducers.push_back(id);
                }
                auto induced_graph = YAGL::induced_subgraph(system_graph, inducers);
                std::map<std::size_t, std::size_t> placeholder;
                //using at, because the [] requires they map_type to be default constructible
                auto rho = gamma.stochastic_rules.at(inst.name).propensity(induced_graph, placeholder);
                propensity_space.push_back({r,rho});
            }

            geocell_propensity = 0.0;
            for(auto& [name, data] : gamma.stochastic_rules)
            {
                auto sum = 0.0;
                for(auto& item : propensity_space)
                {
                    if(rule_instances[item.first].name == name) sum += item.second;
                }
                std::cout << "Rule " << name << " has total propensity of " << sum << "\n";
                geocell_propensity += sum;
            }
            std::cout << "Geocell " << k << " has total propensity of " << geocell_propensity << "\n";

            //the step adapts based on propensity or systems fastest dynamic
            //settings.DELTA_T_MIN = std::min(1.0/(10.0*geocell_propensity), settings.DELTA_DELTA_T);

            // STEP(2) : solve the system of ODES
            //ode_system.step();
            //ode_system.copy_back();
            // Alternate STEP(2) : use forward euler to solve the TAU ODE
            tau += geocell_propensity*settings.DELTA_T_MIN;
            std::cout << "Current tau: " << tau << "\n";
            // STEP(3) : advance the loop timer
            //delta_t = ode_system.t;
            delta_t += settings.DELTA_T_MIN;

            //std::cout << "tout: " << ode_system.t << " , DT: " << delta_t << "\n";
        }

        // if we get over our threshold an event can occur
        if(tau >= exp_sample) {
            //determine which rule to fire and fire it
            auto lam = [](double& s, auto& rp) { return s+rp.second; };
            auto sum = std::accumulate(propensity_space.begin(), propensity_space.end(), 0.0, lam);
            std::cout << "Size of sample space: " << sum << "\n";

            auto local_sample = RandomRealsBetween(0.0, 1.0)();

            auto eventFired = 0;

            auto local_progress = local_sample*geocell_propensity - propensity_space[0].second;
            while(local_progress > 0.0)
            {
                eventFired++;
                local_progress -= propensity_space[eventFired].second;
            }

            auto fired_id = propensity_space[eventFired].first;
            std::cout << "Selected rule id " << fired_id << " of type " << rule_instances[fired_id].name << "\n";

            //TODO: perform the rewrite

            //TODO: incrementally update the in memory set of matches

            //zero out tau since a rule has fired
            tau = 0.0;
            //reset the exp waiting_time sample
            double uniform_sample = RandomRealsBetween(0.0, 1.0)();
            //sample the exponential variable
            exp_sample = -log(1-uniform_sample);
        }
    }
    std::cout << "Total steps taken: " << steps << "\n";
    geocell_progress.first = tau;
    geocell_progress.second = exp_sample;
}

}
#endif //DGGML_APPROXIMATESSA_HPP
