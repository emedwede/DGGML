#ifndef DGGML_APPROXIMATESSA_HPP
#define DGGML_APPROXIMATESSA_HPP

//#include "PlantTypes.hpp"
//#include "PlantGrammar.hpp"

#include "ExpandedComplex2D.hpp"
#include "YAGL_Graph.hpp"
#include "YAGL_Algorithms.hpp"

#include "Solver.h"

#include "ComponentMatchMap.hpp"
#include "IncrementalUpdate.hpp"
#include "AnalyzedGrammar.hpp"
#include "HelperFunctions.hpp"
#include "RandomFunctions.hpp"
#include "HelperStructs.hpp"

#include <chrono>
#include <algorithm>
#include <vector>
#include <iostream>
#include <math.h>
#include <map>

namespace DGGML
{

template <typename T1, typename T2, typename T3, typename T4, typename M1>
void approximate_ssa(ComponentMatchMap<T1>& component_matches, AnalyzedGrammar<T2>& grammar_analysis, T3& rule_map, T4& rule_matches, M1& model,
                     std::size_t k, GeocellProperties<T1>& geocell_properties, CellList<T2>& cell_list)
{
    auto& geoplex2D = model->geoplex2D;
    auto& system_graph = model->system_graph;
    auto& settings = model->settings;
    auto& gen = model->gen;
    auto reaction_radius = model->settings.MAXIMAL_REACTION_RADIUS;

//    std::cout << "Cell " << k << " has " << rule_map[k].size() << " rules\n";

    double delta_t, geocell_propensity;
    std::size_t events, steps;
    delta_t = 0.0; events = 0; steps = 0, geocell_propensity = 0.0;
    double tau = geocell_properties.tau;
    double exp_sample = geocell_properties.exp_sample;

    if(exp_sample <= 0.0)
    {
        double uniform_sample = RandomRealsBetween(0.0, 1.0)();
        exp_sample = -log(1-uniform_sample);
//        std::cout << "Warped wating time: " << exp_sample << "\n";
    }

    std::vector<std::size_t> solving_rules;
    for(auto& item : rule_map[k])
        if(rule_matches[item].category == "deterministic")
            solving_rules.push_back(item);
//    std::cout << "There are " << solving_rules.size() << " solving rules\n";

    // the number of equations is defined by the parameters of the nodes
    // involved in all the ode rule types, it's up to the user to get the number per
    // rule right
    auto num_eq = 0;
    for(auto& item : solving_rules)
    {
        auto& name = rule_matches[item].name;
        num_eq += grammar_analysis.solving_rules.find(name)->second.num_eq;
    }
    num_eq++; //one extra equation for tau

    Solver ode_system(model, grammar_analysis, solving_rules, rule_matches, component_matches,
                      k, geocell_propensity, tau, exp_sample, num_eq);

    //return;
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
                auto& inst = rule_matches[r];
                if(inst.category != "stochastic")
                    continue;
                auto induced_graph = induce_from_set(rule_matches[r], component_matches, system_graph);
                std::map<std::size_t, std::size_t> lhs_vertex_map;
                construct_grammar_match_map(inst, grammar_analysis, lhs_vertex_map, component_matches);
                //using at, because the [] requires they map_type to be default constructible
                auto rho = grammar_analysis.with_rules.at(inst.name).propensity(induced_graph, lhs_vertex_map);
                propensity_space.push_back({r,rho});
            }

            geocell_propensity = 0.0;
            for(auto& [name, data] : grammar_analysis.with_rules)
            {
                auto sum = 0.0;
                for(auto& item : propensity_space)
                {
                    if(rule_matches[item.first].name == name) sum += item.second;
                }
//                std::cout << "Rule " << name << " has total propensity of " << sum << "\n";
                geocell_propensity += sum;
            }
            //std::cout << "Geocell " << k << " has total propensity of " << geocell_propensity << "\n";

            //the step adapts based on propensity or systems fastest dynamic
            //settings.DELTA_T_MIN = std::min(1.0/(10.0*geocell_propensity), settings.DELTA_DELTA_T);

            // STEP(2) : solve the system of ODES
            ode_system.step();
            //TODO: make sure the step is adaptively solving to the event, meaning many events can happen between a
            // delta_delta_t
//            std::cout << "finished stepping\n";
            ode_system.copy_back();
//            std::cout << "finished copying\n";
            // STEP(3) : advance the loop timer
            delta_t = ode_system.t;
            steps++;
            // Alternate STEP(2) : use forward euler to solve the TAU ODE
//            tau += geocell_propensity*settings.DELTA_T_MIN;
//            std::cout << "Current tau: " << tau << "\n";
            //Alternate STEP(3)
            //delta_t += settings.DELTA_T_MIN;

            //std::cout << "tout: " << ode_system.t << " , DT: " << delta_t << "\n";
        }

        // if we get over our threshold an event can occur
        if(tau >= exp_sample) {
            //determine which rule to fire and fire it
            auto lam = [](double& s, auto& rp) { return s+rp.second; };
            auto sum = std::accumulate(propensity_space.begin(), propensity_space.end(), 0.0, lam);
            //std::cout << "Size of sample space: " << sum << "\n";

            auto local_sample = RandomRealsBetween(0.0, 1.0)();

            auto eventFired = 0;

            auto local_progress = local_sample*geocell_propensity - propensity_space[0].second;
            while(local_progress > 0.0)
            {
                eventFired++;
                local_progress -= propensity_space[eventFired].second;
            }

            auto fired_id = propensity_space[eventFired].first;
            auto inst = rule_matches[fired_id];
            auto fired_name = inst.name;
//            std::cout << "Selected rule id " << fired_id << " of type " << fired_name  << " with components: { ";
//            for(auto& c : inst.components) std::cout << c << " ";
//            std::cout << "}\n";
//            if(fired_name == "retraction")
//                std::cin.get();

            //First step, I need to plug in the experimental rewrite code and clean the rest up
            //Problem: rewrite wants a component_match_set, not the components stored in the rule system like I had
            //done originally

            //TODO: a big thing to decide is how phi partitions the match and component set
            auto [changes, removals] = perform_invalidations_and_rewrite(inst, component_matches, gen, grammar_analysis,
                                                             system_graph, rule_matches, rule_map[k], cell_list, geocell_properties);
//            changes.print();
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
            // hierarchy: (geocells) -> rule instances -> components instances -> nodes
            //slow, but easy case to code - assume we know nothing and must search everywhere

            //should invalidate components and rule instances containing invalid components
//            using graph_t = typename std::remove_reference<decltype(system_graph)>::type;
//            auto removals = perform_invalidations<graph_t>(changes, component_matches,
//                                                           grammar_analysis, rule_matches,
//                                                           rule_map[k], cell_list);
//            removals.print();

            find_new_matches(changes, system_graph, component_matches,grammar_analysis,
                             rule_matches, rule_map[k], cell_list, geoplex2D,k, reaction_radius, geocell_properties);


            //update the count since a rule has fired
            geocell_properties.current_rules_fired++;
            geocell_properties.total_rules_fired++;

            //zero out tau since a rule has fired
            tau = 0.0;
            //reset the exp waiting_time sample
            double uniform_sample = RandomRealsBetween(0.0, 1.0)();
            //sample the exponential variable
            exp_sample = -log(1-uniform_sample);
            ode_system.user_data.waiting_time = exp_sample;
            //update the solving rules
            solving_rules.clear();
            for(auto& item : rule_map[k])
                if(rule_matches[item].category == "deterministic")
                    solving_rules.push_back(item);
//            std::cout << "Number of deterministic rules: " << solving_rules.size() << "\n";
            ode_system.reinit(solving_rules);
        }
        //break; //makes only one reaction occur break point
    }
//    std::cout << "Total steps taken: " << steps << "\n";
//    ode_system.print_stats();
    //std::cin.get();
    geocell_properties.tau = tau;
    geocell_properties.exp_sample = exp_sample;
}

}
#endif //DGGML_APPROXIMATESSA_HPP
