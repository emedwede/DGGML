#ifndef PLANT_GRAMMAR_HPP
#define PLANT_GRAMMAR_HPP

#include "PlantTypes.hpp"

#include "YAGL_Graph.hpp"

#include "MathUtils.hpp"

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

// search for retracting ends in dimensional partitions
template <typename GraphType, typename BucketType>
std::vector<std::vector<mt_key_type>> microtubule_retraction_end_matcher(GraphType& graph, BucketType& bucket)
{
    //YAGL::Graph<mt_key_type, MT_NodeData> graph;    
    std::vector<std::vector<mt_key_type>> matches;
    //iterate the whole graph 
    for(auto i : bucket) 
    {
        auto itype = graph.findNode(i)->second.getData().type;

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

template <typename GraphType, typename MatchType>
void microtubule_growing_end_polymerize_solve(GraphType& graph, GraphType& graph_old, MatchType& match)
{
    if(match.size() != 2) return;
    auto i = match[0]; auto j = match[1];
    auto& node_i_data = graph.findNode(i)->second.getData();
    auto& node_j_data = graph.findNode(j)->second.getData();
    auto& node_i_data_old = graph.findNode(i)->second.getData();
    auto& node_j_data_old = graph.findNode(j)->second.getData();
    double DELTA_DELTA_T = 0.1;
    double LENGTH_DIV_FACTOR = 1.2;
    double DIV_LENGTH = 2.0; // make sure this is bigger than any initialized length
    double V_PLUS = 1.0;
    double length_limiter = (1.0 - (calculate_distance(node_i_data_old.position, node_j_data_old.position)/DIV_LENGTH));
    
    for(auto iter = 0; iter < 3; iter++)
    {
        node_i_data.velocity[iter] = V_PLUS*node_i_data_old.unit_vec[iter]*length_limiter;
        node_i_data.position[iter] += node_i_data_old.velocity[iter]*DELTA_DELTA_T; 
    } 
}

template <typename GraphType, typename MatchType>
void microtubule_retraction_end_depolymerize_solve(GraphType& graph, GraphType& graph_old, MatchType& match)
{
    if(match.size() != 2) return;
    auto i = match[0]; auto j = match[1];
    auto& node_i_data = graph.findNode(i)->second.getData();
    auto& node_j_data = graph.findNode(j)->second.getData();
    auto& node_i_data_old = graph.findNode(i)->second.getData();
    auto& node_j_data_old = graph.findNode(j)->second.getData();
    double DELTA_DELTA_T = 0.1;
    double LENGTH_DIV_FACTOR = 1.2;
    double DIV_LENGTH = 2.0; // make sure this is bigger than any initialized length
    double V_PLUS = 1.0;
    double V_MINUS = V_PLUS / 2.0;
    double length_limiter = ((calculate_distance(node_i_data_old.position, node_j_data_old.position)/DIV_LENGTH));
    if(length_limiter <= 0.000001) length_limiter = 0.0; //absolutely needed 
    for(auto iter = 0; iter < 3; iter++)
    {
        node_i_data.velocity[iter] = V_MINUS*node_i_data_old.unit_vec[iter]*length_limiter;
        node_i_data.position[iter] += node_i_data_old.velocity[iter]*DELTA_DELTA_T; 
    } 
}



} // end namespace Plant 
} //end namespace Cajete

#endif 
