#ifndef PLANT_GRAMMAR_HPP
#define PLANT_GRAMMAR_HPP

#include "PlantTypes.hpp"

#include "YAGL_Graph.hpp"

#include <vector>

namespace Cajete
{

namespace Plant 
{

// search for growing ends 
template <typename GraphType>
std::vector<std::vector<mt_key_type>> microtubule_growing_end_matcher(GraphType& graph)
{
    //YAGL::Graph<mt_key_type, MT_NodeData> graph;    
    std::vector<std::vector<mt_key_type>> matches;
    //iterate the whole graph 
    for(auto iter = graph.node_list_begin(); iter != graph.node_list_end(); iter++)
    {
        auto i = iter->first;
        auto itype = iter->second.getData().type;

        if(itype != positive) continue;
        
        for(auto jter = graph.out_neighbors_begin(i); jter != graph.out_neighbors_end(i); jter++)
        {
            auto j = *jter;
            auto jtype = graph.findNode(j)->second.getData().type;
            
            if(jtype != intermediate) continue;
            std::vector<mt_key_type> temp;
            temp.push_back(i);
            temp.push_back(j);
            matches.push_back(temp);
        }
    }

    return matches;
}

// search for growing ends in a dimensional partition 
template <typename GraphType, typename BucketType>
std::vector<std::vector<mt_key_type>> microtubule_growing_end_matcher(GraphType& graph, BucketType& bucket)
{
    //YAGL::Graph<mt_key_type, MT_NodeData> graph;    
    std::vector<std::vector<mt_key_type>> matches;
    //iterate the whole graph 
    for(auto i : bucket)
    {
        auto itype = graph.findNode(i)->second.getData().type;

        if(itype != positive) continue;
        
        for(auto jter = graph.out_neighbors_begin(i); jter != graph.out_neighbors_end(i); jter++)
        {
            auto j = *jter;
            auto jtype = graph.findNode(j)->second.getData().type;
            
            if(jtype != intermediate) continue;
            std::vector<mt_key_type> temp;
            temp.push_back(i);
            temp.push_back(j);
            matches.push_back(temp);
        }
    }

    return matches;
}


// search for retracting ends 
template <typename GraphType>
std::vector<std::vector<mt_key_type>> microtubule_retraction_end_matcher(GraphType& graph)
{
    //YAGL::Graph<mt_key_type, MT_NodeData> graph;    
    std::vector<std::vector<mt_key_type>> matches;
    //iterate the whole graph 
    for(auto iter = graph.node_list_begin(); iter != graph.node_list_end(); iter++)
    {
        auto i = iter->first;
        auto itype = iter->second.getData().type;

        if(itype != negative) continue;
        
        for(auto jter = graph.out_neighbors_begin(i); jter != graph.out_neighbors_end(i); jter++)
        {
            auto j = *jter;
            auto jtype = graph.findNode(j)->second.getData().type;
            
            if(jtype != intermediate) continue;
            std::vector<mt_key_type> temp;
            temp.push_back(i);
            temp.push_back(j);
            matches.push_back(temp);
        }
    }

    return matches;
}

//Simple first attempt a polymerizing
template <typename GraphType>
void microtubule_growing_end_polymerize_rewrite(GraphType& graph, std::vector<mt_key_type>& match)
{
    if(match.size() != 2) return;
    auto i = match[0]; auto j = match[1];

    //TODO: need a unique key generator
    typename GraphType::key_type key = graph.numNodes()+1;

    double x3[3];

    auto x1 = graph.findNode(i)->second.getData().position;
    auto x2 = graph.findNode(j)->second.getData().position;
    
    auto gamma = 0.75;
    for(auto iter = 0; iter < 3; iter++)
    {
        x3[iter] = x2[iter] - ((x2[iter]-x1[iter]) * gamma); 
    }
    
    graph.addNode({key, {{x3[0], x3[1], x3[2]}, {0, 0, 0}, intermediate}});
    graph.removeEdge(i, j);
    graph.addEdge(i, key);
    graph.addEdge(j, key);

}

template <typename GraphType>
void microtubule_growing_end_polymerize_solve(GraphType& graph, std::vector<mt_key_type>& match)
{
    if(match.size() != 2) return;
    auto i = match[0]; auto j = match[1];

    auto& x1 = graph.findNode(i)->second.getData().position;
    auto& x2 = graph.findNode(j)->second.getData().position;
    
    auto dx = 0.2;
    for(auto iter = 0; iter < 3; iter++)
    {
        x1[iter] += (x1[iter] - x2[iter])*dx;  
    }
}


} // end namespace Plant 
} //end namespace Cajete

#endif 
