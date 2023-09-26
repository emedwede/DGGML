#ifndef DGGML_APPROXIMATESSA_HPP
#define DGGML_APPROXIMATESSA_HPP

#include "PlantTypes.hpp"
#include "PlantGrammar.hpp"

#include "ExpandedComplex2D.hpp"
#include "YAGL_Graph.hpp"
#include "YAGL_Algorithms.hpp"

#include "Solver.h"

#include "RuleSystem.hpp"
#include "AnalyzedGrammar.hpp"

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


template<typename T1, typename T2, typename T3>
auto induce_from_set(T1& inst, T2& rule_system, T3& system_graph)
{
    std::vector<std::size_t> inducers;
    for(auto& c : inst.components)
    {
        for(auto& id : rule_system[c])
            inducers.push_back(id);
    }
    return YAGL::induced_subgraph(system_graph, inducers);
}

template <typename T1, typename T2, typename T3, typename T4, typename M1>
void approximate_ssa(RuleSystem<T1>& rule_system, AnalyzedGrammar<T2>& grammar_analysis, T3& rule_map, T4& rule_instances, M1& model,
                                          std::size_t k, std::pair<double, double>& geocell_progress)
{
    auto& geoplex2D = model->geoplex2D;
    auto& system_graph = model->system_graph;
    auto& settings = model->settings;

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
                auto induced_graph = induce_from_set(rule_instances[r], rule_system, system_graph);
                std::map<std::size_t, std::size_t> placeholder;
                //using at, because the [] requires they map_type to be default constructible
                auto rho = grammar_analysis.with_rules.at(inst.name).propensity(induced_graph, placeholder);
                propensity_space.push_back({r,rho});
            }

            geocell_propensity = 0.0;
            for(auto& [name, data] : grammar_analysis.with_rules)
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
            auto inst = rule_instances[fired_id];
            auto fired_name = inst.name;
            std::cout << "Selected rule id " << fired_id << " of type " << fired_name << "\n";

            //TODO: perform the rewrite
            auto& rewrite = grammar_analysis.with_rewrites[fired_name];
            rewrite.print_node_sets(fired_name);
            rewrite.print_edge_sets(fired_name);
            //create a copy of the old lhs
            auto lhs_graph_copy = induce_from_set(inst, rule_system, system_graph);
            auto rhs_graph = grammar_analysis.with_rules.at(fired_name).rhs_graph;
           for(auto& n : lhs_graph_copy.getNodeSetRef())
            {
                std::cout << n.first << "\n";
//                auto res = system_graph.findNode(n.first);
//                std::cout << res->second.getData().position[0] << "\n";
//                std::cout << n.second.getData().position[0] << "\n";
//                n.second.getData().position[0] = 2.2;
//                std::cout << res->second.getData().position[0] << "\n";
//                std::cout << n.second.getData().position[0] << "\n";
            }

           for(auto& item : rewrite.node_set_create)
           {
               auto n = rhs_graph.findNode(item)->second;
               std::cout << n.getData().type << "\n";
               lhs_graph_copy.addNode(n);
           }
            //------Idea------//
            //Step 1: Create two copies of the lhs instance, one for the prev state and the other
            //        for the transformation (Note we could alternatively just compute this ahead of time)
            //Step 2: Transform one of the lhs copies into a rhs through rewrites
            //Step 3: pass a lhs copy into an update with the rhs we just created, plus the f(g(lhs)) vertex mapping
            //Step 4: perform the rewrite on the system graph and copy in the values from the new rhs

            //TODO: incrementally update the in memory set of matches
            //Once we have updated the graph, we need to invalidated any old matches and then search for the
            //new ones!

            //When it comes to match invalidation, removing a node, changing node type, or removing and edge
            //are the actions that lead to invalidation. Why does removal not trigger a re-search? Because our
            // pattern matcher would've already found any other sub pattern if it existed. Changing type, on the
            // other hand, does trigger a research. Adding new a node or edge only opens up a new search path and also does.
            //There are different apporaches to handle finding newly created matches. A potentially memory intensive
            //method would be to remember all completed and partial search paths. My method is a bit more brute force,
            //but only locally.

            //------Idea------//
            //adding/

            return;
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
