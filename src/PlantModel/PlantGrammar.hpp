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

template <typename GraphType, typename BucketType>
std::vector<std::vector<mt_key_type>> three_intermediate_matcher(GraphType& graph, BucketType& bucket)
{
    std::vector<std::vector<mt_key_type>> matches;
    for(auto& i : bucket)
    {
        auto& itype = graph.findNode(i)->second.getData().type;

        if(itype != intermediate) continue;

        for(auto jter = graph.out_neighbors_begin(i); jter != graph.out_neighbors_end(i); jter++)
        {
            auto j = *jter;
            
            auto& jtype = graph.findNode(j)->second.getData().type;
            if(jtype != intermediate) continue;

            //since we are doing a wildcard, the type of j does not matter
            for(auto kter = graph.out_neighbors_begin(i); kter != graph.out_neighbors_end(i); kter++)
            {
                auto k = *kter;
                
                auto& ktype = graph.findNode(k)->second.getData().type;
                if(ktype != intermediate) continue;

                if(j != k) //as long as j and k are not the same, we have a match 
                {
                    matches.push_back({{i, j, k}});
                }
            }
        }
    }

    return matches;
}

//search for wild carded match types
template <typename GraphType, typename BucketType>
std::vector<std::vector<mt_key_type>> wildcard_intermediate_wildcard_matcher(GraphType& graph, BucketType& bucket)
{
    std::vector<std::vector<mt_key_type>> matches;
    for(auto& i : bucket)
    {
        auto& itype = graph.findNode(i)->second.getData().type;

        if(itype != intermediate) continue;

        for(auto jter = graph.out_neighbors_begin(i); jter != graph.out_neighbors_end(i); jter++)
        {
            auto j = *jter;
            //since we are doing a wildcard, the type of j does not matter
            for(auto kter = graph.out_neighbors_begin(i); kter != graph.out_neighbors_end(i); kter++)
            {
                auto k = *kter;
                if(j != k) //as long as j and k are not the same, we have a match 
                {
                    matches.push_back({{i, j, k}});
                }
            }
        }
    }

    return matches;
}


template <typename GraphType, typename MatchType>
void collision_match_refinement(GraphType& graph, double cutoff, MatchType& growing_matches, MatchType& wildcard_matches, MatchType& collision_matches)
{
    std::vector<mt_key_type> temp; temp.reserve(6);

    //loop over all the matches and ensure that they
    //do not share any keys
    for (auto& match_g : growing_matches)
    {
        for(auto& match_w : wildcard_matches)
        {
            bool valid = true;
            for(auto& i : match_g)
            {
                for(auto& j : match_w)
                {
                    if(i == j) valid = false;
                }
            }
            //as long as the match is valid, return it only if it passes
            //the distance check
            if(valid) 
            {
                //we choose to anchor the subgraph match at a point in order to compute the 
                //nearness
                auto& anchor_pos_g = graph.findNode(match_g[0])->second.getData().position;
                auto& anchor_pos_w = graph.findNode(match_w[0])->second.getData().position;

                auto distance = calculate_distance(anchor_pos_g, anchor_pos_w);
                if(distance <= cutoff)
                {
                    for(auto& i : match_g) temp.push_back(i);
                    for(auto& j : match_w) temp.push_back(j);
                    collision_matches.push_back(temp);
                    temp.clear();    
                }
            }
        }
    }
}

//Simple first attempt a polymerizing
template <typename GraphType>
void microtubule_growing_end_polymerize_rewrite(GraphType& graph, std::vector<mt_key_type>& match)
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
    
    graph.addNode({key, 
            {{x3[0], x3[1], x3[2]}, 
            {0, 0, 0}, 
            intermediate, 
            {-1, -1, -1}, 
            {u1[0], u1[1], u1[2]}}});

    graph.removeEdge(i, j);
    graph.addEdge(i, key);
    graph.addEdge(j, key);
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

template <typename GraphType, typename BucketType>
void microtubule_positive_to_negative_rewrite(GraphType& graph, std::vector<mt_key_type>& match, BucketType& bucket)
{
    auto& data = graph.findNode(match[0])->second.getData();
    //change type and switch directions
    data.type = negative;
    for(int i = 0; i < 3; i++) data.unit_vec[i] *= -1;
}

template <typename GraphType, typename BucketType>
void microtubule_negative_to_positive_rewrite(GraphType& graph, std::vector<mt_key_type>& match, BucketType& bucket)
{
    auto& data = graph.findNode(match[0])->second.getData();
    //change type and switch directions
    data.type = positive;
    for(int i = 0; i < 3; i++) data.unit_vec[i] *= -1;
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

template <typename GraphType, typename ParamType>
double microtubule_positive_to_negative_propensity(GraphType& graph, std::vector<mt_key_type>& match, ParamType& settings)
{
    return 0.2;
}

template <typename GraphType, typename ParamType>
double microtubule_negative_to_positive_propensity(GraphType& graph, std::vector<mt_key_type>& match, ParamType& settings)
{
    return 0.2;
}

template <typename GraphType, typename BucketType, typename ParamType>
void microtubule_crossover_rewrite(GraphType& graph, std::vector<mt_key_type>& match, BucketType& bucket, ParamType& settings)
{
    //TODO: need a unique key generator
    typename GraphType::key_type key = graph.numNodes()+1;
    
    while(graph.findNode(key) != graph.node_list_end()) key++; //TODO: fix, very greedy
    
    auto x1 = match[2]; //intermediate wildcard rule
    auto x2 = match[3]; //wildcard left relative to intermediate
    auto x3 = match[4]; //wildcard right relative to intermediate
    auto x4 = match[1]; //intermediate growing end
    auto x5 = match[0]; //positive growing end
   
    //find all the node data
    auto& dat1 = graph.findNode(x1)->second.getData();
    auto& dat2 = graph.findNode(x2)->second.getData();
    auto& dat3 = graph.findNode(x3)->second.getData();
    auto& dat4 = graph.findNode(x4)->second.getData();
    auto& dat5 = graph.findNode(x5)->second.getData();
   
    //get references to position vector
    auto& pos1 = dat1.position;
    auto& pos2 = dat2.position;
    auto& pos3 = dat3.position;
    auto& pos4 = dat4.position;
    auto& pos5 = dat5.position;

    //get references to unit vectors
    auto& u1 = dat1.unit_vec;
    auto& u2 = dat2.unit_vec;
    auto& u3 = dat3.unit_vec;
    auto& u4 = dat4.unit_vec;
    auto& u5 = dat5.unit_vec;

    set_unit_vector(pos1, pos4, u5);

    //place x5 past the point of the intersecting node, keep unit vec the same
    pos5[0] = pos1[0] + settings.MT_MIN_SEGMENT_INIT*u5[0];
    pos5[1] = pos1[1] + settings.MT_MIN_SEGMENT_INIT*u5[1];
     
    //create an intermediate node between the new x5 and x1
    int64_t k = bucket.first;
    graph.addNode({key, 
            {{pos1[0]+(pos5[0]-pos1[0])/2.0, pos1[1]+(pos5[1]-pos1[1])/2.0, pos1[2]+(pos5[2]-pos1[2])/2.0}, 
            {0, 0, 0}, 
            intermediate, 
            {k, k, k}, 
            {u5[0], u5[1], u5[2]}}});
    bucket.second.push_back(key);

    //rewrite x5s connections and connect up the new node, key
    graph.removeEdge(x5, x4);
    graph.addEdge(x5, key);
    graph.addEdge(key, x1);
    graph.addEdge(x1, x4);
    
    set_unit_vector(pos1, pos4, u4);
    //chage the type of x1
    dat1.type = junction;
}

template <typename GraphType, typename BucketType, typename ParamType>
void microtubule_zippering_rewrite(GraphType& graph, std::vector<mt_key_type>& match, BucketType& bucket, ParamType& settings)
{
    auto x1 = match[2]; //intermediate wildcard rule
    auto x4 = match[1]; //intermediate growing end
    auto x5 = match[0]; //positive growing end
   
    //find all the node data
    auto& dat1 = graph.findNode(x1)->second.getData();
    auto& dat4 = graph.findNode(x4)->second.getData();
   
    //get references to position vector
    auto& pos1 = dat1.position;
    auto& pos4 = dat4.position;

    //get references to unit vectors
    auto& u4 = dat4.unit_vec;
    
    //remove the positive node and connect it's intermediate to
    //the threeway juntion
    graph.removeNode(graph.findNode(x5)->second);
    graph.addEdge(x1, x4);
    set_unit_vector(pos1, pos4, u4);
    
    dat1.type = zipper;
    
    //remove deleted node from the bucket
    std::size_t found;
    for(auto iter = 0; iter < bucket.second.size(); iter++)
    {
        if(bucket.second[iter] == x5)
        {
            found = iter;
            break;    
        }
    }
    //TODO: improve this O(N) search
    bucket.second.erase(bucket.second.begin()+found); //remove from bucket

}


template <typename GraphType, typename BucketType, typename ParamType>
void microtubule_catastrophe_rewrite(GraphType& graph, std::vector<mt_key_type>& match, BucketType& bucket, ParamType& settings)
{
    auto x5 = match[0]; //positive growing end
    auto& dat5 = graph.findNode(x5)->second.getData();
    auto& pos5 = dat5.position;
    auto& u5 = dat5.unit_vec;
    
    //change the positive type to negative and reverse the unit vector
    dat5.type = negative;
    u5[0] = -u5[0]; u5[1] = -u5[1]; u5[2] = -u5[2];
}


template <typename GraphType, typename ParamType>
double microtubule_collision_crossover_propensity(GraphType& graph, std::vector<mt_key_type>& match, ParamType& settings)
{
    auto x1 = match[2]; //intermediate wildcard rule
    auto x2 = match[3]; //wildcard left relative to intermediate
    auto x3 = match[4]; //wildcard right relative to intermediate
    auto x4 = match[1]; //intermediate growing end
    auto x5 = match[0]; //positive growing end
   
    //find all the node data
    auto& dat1 = graph.findNode(x1)->second.getData();
    auto& dat2 = graph.findNode(x2)->second.getData();
    auto& dat3 = graph.findNode(x3)->second.getData();
    auto& dat4 = graph.findNode(x4)->second.getData();
    auto& dat5 = graph.findNode(x5)->second.getData();
   
    //get references to position vector
    auto& pos1 = dat1.position;
    auto& pos2 = dat2.position;
    auto& pos3 = dat3.position;
    auto& pos4 = dat4.position;
    auto& pos5 = dat5.position;

    //get references to unit vectors
    auto& u1 = dat1.unit_vec;
    auto& u2 = dat2.unit_vec;
    auto& u3 = dat3.unit_vec;
    auto& u4 = dat4.unit_vec;
    auto& u5 = dat5.unit_vec;
    
    double sol_l[2];
    double sol_r[2];
    double theta_l = unit_dot_product(u5, u2);
    double theta_r = unit_dot_product(u5, u3);

    double propensity = 0.0;

    if(theta_l != 1.0 && theta_l != -1.0)
    {
        paramaterized_intersection(pos5, pos2, pos1, u5, sol_l);
        
        if(sol_l[0] > 0.0 && sol_l[1] >= 0.0 && sol_l[1] <= 1.0) 
        {
            std::cout << "Sol l:" << sol_l[0] << " " << sol_l[1] << "\n";
            propensity = 
                exp(-pow(calculate_distance(pos1, pos5), 2.0) / pow(0.5*settings.DIV_LENGTH, 2.0)); 
            return propensity;
        } 
    } 
    else if(theta_r != 1.0 && theta_r != -1.0)
    {
        paramaterized_intersection(pos5, pos3, pos1, u5, sol_r);
        
        if(sol_r[0] > 0.0 && sol_r[1] >= 0.0 && sol_r[1] <= 1.0) 
        {
            std::cout << "Sol r:" << sol_r[0] << " " << sol_r[1] << "\n";
            propensity = 
                exp(-pow(calculate_distance(pos1, pos5), 2.0) / pow(0.5*settings.DIV_LENGTH, 2.0));
            return propensity;
        }
    } 
    else 
    {
        propensity = 0.0;
        return propensity;
    }



    //propensity = exp(-pow(calculate_distance(pos1, pos5), 2.0) / pow(0.5*settings.DIV_LENGTH, 2.0));
    
    return propensity;
    /*double gamma_l = 
        cross_product((pos1[0]-pos2[0]), (pos1[1]-pos2[1]), (pos2[0]-pos5[0]), (pos2[1]-pos5[1]))
        /
        cross_product((pos1[0]-pos2[0]), (pos1[1]-pos2[1]), u5[0], u5[1]);

    double gamma_r = 
        cross_product((pos3[0]-pos2[0]), (pos3[1]-pos2[1]), (pos2[0]-pos5[0]), (pos2[1]-pos5[1]))
        /
        cross_product((pos3[0]-pos2[0]), (pos3[1]-pos2[1]), u5[0], u5[1]);

    double alpha_l = 
        -cross_product((pos2[0]-pos5[0]), (pos2[1]-pos5[1]), u5[0], u5[1])
        /
        cross_product((pos1[0]-pos2[0]), (pos1[1]-pos2[1]), u5[0], u5[1]);
    
    double alpha_r = 
        -cross_product((pos2[0]-pos5[0]), (pos2[1]-pos5[1]), u5[0], u5[1])
        /
        cross_product((pos3[0]-pos2[0]), (pos3[1]-pos2[1]), u5[0], u5[1]);
    
    double alpha_1 = calculate_alpha(pos2, pos1, pos5);
    double alpha_2 = calculate_alpha(pos2, pos3, pos5);

    double beta_1 = calculate_beta(pos2, pos1, pos5, alpha_1);
    double beta_2 = calculate_beta(pos2, pos3, pos5, alpha_2);

    double beta_max = std::min(beta_1, beta_2);
    
    //if(beta_max >= CYLINDRICAL_CUTOFF*V_PLUS) return 0.0;
    
    if(gamma_l <= 0.0 || gamma_r <= 0.0) return 0.0;

    double incoming_angle = unit_dot_product(u4, u1);
    
    if(incoming_angle < 0.0) incoming_angle = -incoming_angle;

    double THETA_CRIT = 0.5;

    double critical_angle = cos(THETA_CRIT);

    int crit_out_range = 0;
    
    if(incoming_angle < 1.0 && incoming_angle >= critical_angle) 
    {
        crit_out_range += 1;
    }

    double e = exp(-1.0*beta_max*beta_max/2.0);
    
    int in_range = 0;

    double EPSILON = 0.2;
    if((alpha_r >= EPSILON) && (alpha_r <= 1.0 - EPSILON)) 
    {
        in_range = 1;
    }

    propensity = e*in_range;//crit_out_range 
    
    return propensity;*/
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

    if(length_limiter <= 0.01) length_limiter = 0.0; //absolutely needed 

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
