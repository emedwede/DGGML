#ifndef CAJETE_PLANT_SSA_HPP
#define CAJETE_PLANT_SSA_HPP 

#include "PlantTypes.hpp"
#include "PlantGrammar.hpp"

#include "ExpandedComplex2D.hpp"
#include "YAGL_Graph.hpp"

#include <chrono>
#include <random>
#include <algorithm>
#include <vector>
#include <iostream>
#include <math.h>
#include <map>

namespace Cajete 
{
namespace Plant 
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
        [distribution_ = std::uniform_real_distribution<double>(low, high),
         random_engine_ = std::mt19937{std::random_device{}() }]() mutable
       {
           return distribution_(random_engine_);
       };
    return randomFunc;
};

//template <typename BucketType>
//void plant_model_ssa(
//        BucketType& bucket,
//        ExpandedComplex2D<>& geoplex2D, 
//        YAGL::Graph<mt_key_type, MT_NodeData>& system_graph) 
template <typename BucketType, typename GeoplexType, typename GraphType, typename ParamType>
void plant_model_ssa(BucketType& bucket, GeoplexType& geoplex2D, GraphType& system_graph, ParamType& settings)
{
    double delta_t, exp_sample, tau, geocell_propensity;
    std::size_t events;
    
    delta_t = 0.0; events = 0;
    while(delta_t < settings.DELTA) 
    {
        //reset tau
        tau = 0.0;
        
        auto k = bucket.first;
        std::size_t num_patterns = 2;
        std::vector<std::vector<mt_key_type>> rule_matches[num_patterns];
        rule_matches[0] = microtubule_growing_end_matcher(system_graph, bucket.second); 
        rule_matches[1] = microtubule_retraction_end_matcher(system_graph, bucket.second); 
        
        double uniform_sample = RandomRealsBetween(0.0, 1.0)();
        
        //sample the exponential variable
        exp_sample = -log(1-uniform_sample);
        
        
        while(delta_t < settings.DELTA && tau < exp_sample)
        {
            // STEP(0) : store old state and find all the matches 
            //do a deep copy? TODO: only need to copy old state parameters
            auto system_graph_old = system_graph;
                   
            // STEP(1) : sum all of the rule propensities
            //zero the geocell_propensity
            geocell_propensity = 5.0;
            
            // STEP(2) : solve the system of ODES 
            microtubule_ode_solver(rule_matches, system_graph, system_graph_old, k, settings); 
                        
            // STEP(3) : use forward euler to solve the TAU ODE
            tau += geocell_propensity*settings.DELTA_DELTA_T; //TODO: we need to be careful not to oversolve
            
            // STEP(4) : advance the loop timer
            delta_t += settings.DELTA_DELTA_T; //TODO: make delta_t adaptive
        }
    

        // if we get over our threshold an event can occur
        if(tau > exp_sample) 
        {
            //determine which rule to file and fire it
            auto k = bucket.first;
            microtubule_rule_firing(rule_matches, system_graph, k);
        }
    }
}

//be careful with solving, could lead to a segmentation fault 
//in the sorting phase if a parameter is solved out to out of bounds
template <typename MatchSetType, typename GraphType, typename KeyType, typename ParamType>
void microtubule_ode_solver(MatchSetType* all_matches, 
        GraphType& system_graph,
        GraphType& system_graph_old,
        KeyType& k,
        ParamType& settings)
{
    for(auto i = 0; i < 2; i++) 
    {
        for(auto match : all_matches[i])
        {
           bool bad_match = false;
           //check match integrity
           for(auto key : match)
           {
                auto dtag = system_graph.findNode(key)->second.getData().tagND[0];
                if(dtag != k)
                {
                    bad_match = true;
                    //std::cout << "Bad match found, it'll be skipped!\n";
                    break;
                }
           }

           if(!bad_match)
           {
               if(i == 0)
                    microtubule_growing_end_polymerize_solve(system_graph, system_graph_old, match, settings);
               if(i == 1)
                   microtubule_retraction_end_depolymerize_solve(system_graph, system_graph_old, match, settings);
           }
        }
    }
}

template <typename MatchType, typename GraphType, typename KeyType>
void microtubule_rule_firing(MatchType* all_matches, GraphType& system_graph, KeyType& k)
{
    std::vector<std::vector<mt_key_type>> good_matches;
    for(auto i = 0; i < 2; i++) 
    {
        //check match integrity and only fire the good one
        for(auto match : all_matches[i])
        {
           bool bad_match = false;
           //check match integrity
           for(auto key : match)
           {
                auto dtag = system_graph.findNode(key)->second.getData().tagND[0];
                if(dtag != k)
                {
                    bad_match = true;
                    break;
                }
           }
           if(!bad_match)
           {
               good_matches.push_back(match);
           }
        }
        
        if(good_matches.size() == 0) return; //we cant proced without good matches
        auto rule_fired = RandomIntsBetween(0, good_matches.size())();
        if(i == 0)
            microtubule_growing_end_polymerize_rewrite(system_graph, good_matches[rule_fired]); 
        //if(i == 1)
            //microtubule_retraction_end_depolymerize_rewrite();
    }
}
    
} //end namespace plant 
} //end namespace cajete

#endif
