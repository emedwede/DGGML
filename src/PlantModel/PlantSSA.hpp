#ifndef CAJETE_PLANT_SSA_HPP
#define CAJETE_PLANT_SSA_HPP 

#include "PlantTypes.hpp"
#include "PlantGrammar.hpp"

#include "ExpandedComplex2D.hpp"
#include "YAGL_Graph.hpp"

#include <random>
#include <map>

namespace Cajete 
{
namespace Plant 
{
template <typename BucketType>
void plant_model_ssa(
        BucketType& bucket,
        ExpandedComplex2D<>& geoplex2D, 
        YAGL::Graph<mt_key_type, MT_NodeData>& system_graph) 
{
        
    auto k = bucket.first;
    
    std::cout << "Running ssa and matcher on geocell " << k << "\n";
    auto matches = microtubule_growing_end_matcher(system_graph, bucket.second); 
    std::cout << "Found " << matches.size() << " candidates\n";
    for(auto match : matches)
    {
       bool bad_match = false;
       //check match integrity
       for(auto key : match)
       {
            auto dtag = system_graph.findNode(key)->second.getData().tagND[0];
            if(dtag != k)
            {
                bad_match = true;
                std::cout << "Bad match found, it'll be skipped!\n";
                break;
                }
       }
       if(!bad_match)
       {
            microtubule_growing_end_polymerize_solve(system_graph, match);
            microtubule_growing_end_polymerize_rewrite(system_graph, match); 
       }
    }
}

} //end namespace plant 
} //end namespace cajete

#endif
