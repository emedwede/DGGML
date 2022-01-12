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

// search for retracting ends in dimensional partitions
template <typename GraphType, typename BucketType>
std::vector<std::vector<mt_key_type>> microtubule_retraction_end_two_intermediate_matcher(GraphType& graph, BucketType& bucket)
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
            
            for(auto kter =graph.out_neighbors_begin(j); kter != graph.out_neighbors_end(j); kter++)
            {
                auto k = *kter;
                auto ktype = graph.findNode(k)->second.getData().type;

                if(ktype != intermediate) continue;
                std::vector<mt_key_type> temp;
                temp.push_back(i);
                temp.push_back(j);
                temp.push_back(k);
                matches.push_back(temp); 
            }
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
template <typename GraphType, typename BucketType>
void microtubule_growing_end_polymerize_rewrite(GraphType& graph, std::vector<mt_key_type>& match, BucketType& bucket)
{
    if(match.size() != 2) return;
    auto i = match[0]; auto j = match[1];

    //TODO: need a unique key generator
    typename GraphType::key_type key = graph.numNodes()+1;
    while(graph.findNode(key) != graph.node_list_end()) key++; //TODO: fix, very greedy
    double x3[3];

    auto& x1 = graph.findNode(i)->second.getData().position;
    auto& x2 = graph.findNode(j)->second.getData().position;
    auto& u1 = graph.findNode(i)->second.getData().unit_vec;
    auto gamma = 0.75;
    for(auto iter = 0; iter < 3; iter++)
    {
        x3[iter] = x2[iter] - ((x2[iter]-x1[iter]) * gamma); 
    }
    
    int64_t k = bucket.first;
    graph.addNode({key, 
            {{x3[0], x3[1], x3[2]}, 
            {0, 0, 0}, 
            intermediate, 
            {k, k, k}, 
            {u1[0], u1[1], u1[2]}}});

    bucket.second.push_back(key);
    graph.removeEdge(i, j);
    graph.addEdge(i, key);
    graph.addEdge(j, key);

}

template <typename GraphType, typename ParamType>
double microtubule_growing_end_polymerize_propensity(GraphType& graph, std::vector<mt_key_type>& match, ParamType& settings)
{
    auto& node_i_data = graph.findNode(match[0])->second.getData();
    auto& node_j_data = graph.findNode(match[1])->second.getData();
    
    auto len = calculate_distance(node_i_data.position, node_j_data.position);
    //double propensity = heaviside(len, settings.DIV_LENGTH);
    double propensity = sigmoid((len/settings.DIV_LENGTH) - 1.0, settings.SIGMOID_K);
    return propensity;
}

template <typename GraphType, typename ParamType>
double microtubule_retraction_end_depolymerize_propensity(GraphType& graph, std::vector<mt_key_type>& match, ParamType& settings)
{
    auto& node_i_data = graph.findNode(match[0])->second.getData();
    auto& node_j_data = graph.findNode(match[1])->second.getData();
    
    auto len = calculate_distance(node_i_data.position, node_j_data.position);
    //double propensity = heaviside(-len, settings.DIV_LENGTH_RETRACT);
    double propensity = sigmoid(-(len/settings.DIV_LENGTH), settings.SIGMOID_K);
    return propensity;
}

template <typename GraphType, typename MatchType, typename ParamType>
void microtubule_growing_end_polymerize_solve(GraphType& graph, GraphType& graph_old, MatchType& match, ParamType& settings)
{
    if(match.size() != 2) return;

    auto dtdt = settings.DELTA_DELTA_T;
    auto l_d_f = settings.LENGTH_DIV_FACTOR;
    auto d_l = settings.DIV_LENGTH;
    auto v_plus = settings.V_PLUS;

    auto i = match[0]; auto j = match[1];
    auto& node_i_data = graph.findNode(i)->second.getData();
    auto& node_j_data = graph.findNode(j)->second.getData();
    auto& node_i_data_old = graph.findNode(i)->second.getData();
    auto& node_j_data_old = graph.findNode(j)->second.getData();

    double length_limiter = 
        (1.0 - (calculate_distance(node_i_data_old.position, node_j_data_old.position)/d_l));
    
    for(auto iter = 0; iter < 3; iter++)
    {
        node_i_data.velocity[iter] = v_plus*node_i_data_old.unit_vec[iter]*length_limiter;
        node_i_data.position[iter] += node_i_data_old.velocity[iter]*dtdt; 
    } 
}

template <typename GraphType, typename MatchType, typename ParamType>
void microtubule_retraction_end_depolymerize_solve(GraphType& graph, GraphType& graph_old, MatchType& match, ParamType& settings)
{
    if(match.size() != 2) return;
    
    auto dtdt = settings.DELTA_DELTA_T;
    auto l_d_f = settings.LENGTH_DIV_FACTOR;
    auto d_l = settings.DIV_LENGTH;
    auto v_minus = settings.V_MINUS;


    auto i = match[0]; auto j = match[1];
    auto& node_i_data = graph.findNode(i)->second.getData();
    auto& node_j_data = graph.findNode(j)->second.getData();
    auto& node_i_data_old = graph.findNode(i)->second.getData();
    auto& node_j_data_old = graph.findNode(j)->second.getData();

    double length_limiter = 
        ((calculate_distance(node_i_data_old.position, node_j_data_old.position)/d_l));

    if(length_limiter <= 0.000001) length_limiter = 0.0; //absolutely needed 

    for(auto iter = 0; iter < 3; iter++)
    {
        node_i_data.velocity[iter] = v_minus*node_i_data_old.unit_vec[iter]*length_limiter;
        node_i_data.position[iter] += node_i_data_old.velocity[iter]*dtdt; 
    } 
}

//Simple first attempt at depolymerizing
template <typename GraphType, typename BucketType>
void microtubule_retraction_end_depolymerize_rewrite(GraphType& graph, std::vector<mt_key_type>& match, BucketType& bucket)
{
    if(match.size() != 3) return;
    auto i = match[0]; auto j = match[1]; auto k = match[2];

    auto& x1 = graph.findNode(i)->second.getData().position;
    auto& x3 = graph.findNode(k)->second.getData().position;
    auto& u1 = graph.findNode(i)->second.getData().unit_vec;
    auto& u2 = graph.findNode(j)->second.getData().unit_vec;
    auto len = calculate_distance(x1, x3);

    //calculate the new unit vector
    for(auto iter = 0; iter < 3; iter++) u1[iter] = (x3[iter] - x1[iter])/len;
    
    auto node_j = graph.findNode(j)->second;
    graph.removeNode(node_j);
    graph.addEdge(i, k);
    std::size_t found;
    for(auto iter = 0; iter < bucket.second.size(); iter++)
    {
        if(bucket.second[iter] == j)
        {
            found = iter;
            break;    
        }
    }
    //TODO: improve this O(N) search
    bucket.second.erase(bucket.second.begin()+found); //remove from bucket

}



} // end namespace Plant 
} //end namespace Cajete

#endif 
