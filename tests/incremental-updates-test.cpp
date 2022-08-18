#include <iostream>

#include "catch.hpp"

#include "YAGL_Graph.hpp"
#include "YAGL_Algorithms.hpp"
#include "PlantTypes.hpp"
#include "PlantGrammar.hpp"
#include "PlantUtils.hpp"
#include <vector> 
#include <random> 
#include <unordered_map> 
#include <set>
#include <algorithm> 
#include <utility> 

using key_type = Cajete::Plant::mt_key_type; 
using node_type = Cajete::Plant::MT_NodeData;

using graph_type = YAGL::Graph<key_type, node_type>;

void print_matches(std::vector<std::vector<key_type>>& matches)
{
    std::cout << "Found " << matches.size() << " Matches: \n";
    for(auto& m : matches)
    {
        std::cout << "\t{ ";
        for(auto& v : m)
        {
            std::cout << v << " ";
        } std::cout << "}\n";
    }
}

void print_matches(std::unordered_map<key_type, std::vector<key_type>>& matches)
{
    std::cout << "Found " << matches.size() << " Matches: \n";
    for(auto& [key, m] : matches)
    {
        std::cout << "\t{ ";
        for(auto& v : m)
        {
            std::cout << v << " ";
        } std::cout << "}\n";
    }
}

template <typename GraphType>
void microtubule_scatter(GraphType& graph, std::size_t num_mt)
{
    double epsilon_min = 0.5;
    double epsilon_max = 1.0;
    std::random_device random_device;
    std::mt19937 random_engine(random_device());
    
    std::uniform_real_distribution<double> distribution_global_x(0, 10);

    std::uniform_real_distribution<double> distribution_global_y(0, 10);

    std::uniform_real_distribution<double> distribution_local(epsilon_min, epsilon_max);
    
    std::uniform_real_distribution<double> distribution_angle(0.0, 2.0*3.14);
    using node_type = typename GraphType::node_type;

    std::size_t segments = 3;
    for(auto i = 0; i < num_mt ; i++) 
    {
        double x_c, y_c;
        double z_c = 0;
        x_c = distribution_global_x(random_engine);
        y_c = distribution_global_y(random_engine);
        auto theta = distribution_angle(random_engine);
        auto seg_len = distribution_local(random_engine);

        auto x_s = 0.0;
        auto y_s = seg_len;
        auto x_r_t = x_s*cos(theta) + y_s*sin(theta);
        auto y_r_t = -x_s*sin(theta) +y_s*cos(theta);
        auto x_r = x_c + x_r_t;
        auto y_r = y_c + y_r_t;
        auto z_r = 0.0;

        auto x_l = x_c - (x_r - x_c);
        auto y_l = y_c - (y_r - y_c);
        auto z_l = 0.0;
        
        //compute dist and unit vector
        double p1[3] = {x_r - x_c, y_r - y_c, 0.0};
        double p2[3] = {0.0, 0.0, 0.0};
        double p3[3] = {x_c - x_l, y_c - y_l, 0.0};
        double u1[3], u2[3];

        Cajete::set_unit_vector(p1, p2, u1);
        Cajete::set_unit_vector(p3, p2, u2);

        node_type node_l(i*segments, 
                {{x_l, y_l, z_l}, 
                {0.0, 0.0, 0.0}, 
                Cajete::Plant::mt_type::negative, 
                {-1, -1, -1}, 
                {u2[0], u2[1], u2[2]}});

        node_type node_c(i*segments+1, 
                {{x_c, y_c, z_c}, 
                {0.0, 0.0, 0.0}, 
                Cajete::Plant::mt_type::intermediate, 
                {-1, -1, -1}, 
                {u1[0], u1[1], u1[2]}});
        
        node_type node_r(i*segments+2, 
                {{x_r, y_r, z_r}, 
                {0.0, 0.0, 0.0}, 
                Cajete::Plant::mt_type::positive, 
                {-1, -1, -1}, 
                {u1[0], u1[1], u1[2]}});

        graph.addNode(node_l);
        graph.addNode(node_c);
        graph.addNode(node_r);

        graph.addEdge(node_l, node_c);
        graph.addEdge(node_r, node_c);
    }
}


//Simple first attempt a polymerizing
template <typename GraphType>
std::pair<std::set<key_type>, std::set<key_type>> test_rewrite(GraphType& graph, std::vector<key_type>& match)
{
    //if(match.size() != 2) return;
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
    std::random_device random_device; std::mt19937 random_engine(random_device());
    std::uniform_real_distribution<double> distribution_angle(-3.14/8.0, 3.14/8.0); 
    double theta = distribution_angle(random_engine);
    auto u10_rot = u1[0]*cos(theta) + u1[1]*sin(theta);
    auto u11_rot = -u1[0]*sin(theta) +u1[1]*cos(theta);
    u1[0] = u10_rot; u1[1] = u11_rot;

    graph.addNode({key, 
            {{x3[0], x3[1], x3[2]}, 
            {0, 0, 0}, 
            Cajete::Plant::mt_type::intermediate, 
            {0, 0, 0}, 
            {u1[0], u1[1], u1[2]}}});

    graph.removeEdge(i, j);
    graph.addEdge(i, key);
    graph.addEdge(j, key);
    
    return std::pair<std::set<key_type>, std::set<key_type>>{{i, j}, {i, j, key}}; //matches to invalidate and induce
}



TEST_CASE("Incremental Update Test", "[update-test]")
{
    graph_type graph;
   
    std::size_t n = 10; // num mt
    microtubule_scatter(graph, n);

    std::vector<Cajete::Plant::mt_key_type> bucket;
    
    for(const auto& [key, value] : graph.getNodeSetRef())
        bucket.push_back(key);
    
    // compute all the matches
    auto matches = Cajete::Plant::microtubule_growing_end_matcher(graph, bucket);
    //print_matches(matches);
    REQUIRE(matches.size() == n);

    //randomly select a rewrite to occur 
    std::random_device random_device;
    std::mt19937 random_engine(random_device());
    
    std::uniform_int_distribution<std::size_t> rand_rewrite(0, n-1);

    auto selected = rand_rewrite(random_engine);
    std::cout << "Rule selected to rewrite: " << selected << "\n";
    
    std::unordered_map<std::size_t, std::vector<key_type>> gi_map;
    for(std::size_t i = 0; i < matches.size(); i++)
    {
        gi_map.insert({i, matches[i]});
    }
    print_matches(gi_map); 
    //each node needs to know what match it belongs to 
    std::unordered_map<key_type, std::vector<std::size_t>> gi_inverse;
    for(const auto& key : bucket)
       gi_inverse.insert({key, {}}); 
    
    for(const auto& [key, match] : gi_map)
    {
        for(const auto& item : match)
        {
            auto search = gi_inverse.find(item);
            if(search != gi_inverse.end())
            {
                search->second.push_back(key);
            }
        }
    }

    REQUIRE(gi_map.size() == n);
    
    auto lambda_print = [](std::unordered_map<key_type, std::vector<std::size_t>>& gi_inverse)
    {
        for(const auto& [key, value] : gi_inverse)
        {
            std::cout << "Node " << key << ": {";
            for(const auto& node : value)
                std::cout << " " << node << " ";
            std::cout << "}\n";
        }
    };
    lambda_print(gi_inverse);
    // rewrite the graph and gather the nodes the rule invalidates
    auto [invalidations, inducers] = test_rewrite(graph, matches[selected]);
    
    //first test the graph was rewritten
    REQUIRE(YAGL::connected_components(graph) == n);
    // only a single new unique node is added 
    REQUIRE(graph.numNodes() == n*3+1);
    
    REQUIRE(invalidations.size() == 2);

    //invalidate old matches
    for(const auto& i : invalidations)
    {
        const auto& rules = gi_inverse.find(i);
        if(rules != gi_inverse.end())
        {
            //TODO: need to remove matches from inverse map in a less destructive way 
            for(const auto& j : rules->second)
            {
                gi_map.erase(j);
            }
            gi_inverse.erase(i);
        }
    }
    
    REQUIRE(gi_map.size() == n - 1);
    print_matches(gi_map);

    auto temp = inducers;
    // get all nodes one hop away and add those to inducers
    for(const auto& i : temp)
    {
        auto res = YAGL::recursive_dfs(graph, i, 2);
        for(auto& j : res)
            inducers.insert(j);
    }
    REQUIRE(inducers.size() == 4);
    
    //induce a subgraph to search for new matches
    auto subgraph = YAGL::induced_subgraph(graph, inducers);

    REQUIRE(subgraph.numNodes() == 4);
    
    std::vector<Cajete::Plant::mt_key_type> sub_bucket;
    
    for(const auto& [key, value] : subgraph.getNodeSetRef())
        sub_bucket.push_back(key);
    
    // compute all the matches
    auto incremental_matches = Cajete::Plant::microtubule_growing_end_matcher(graph, sub_bucket);
    //print_matches(matches);
    REQUIRE(incremental_matches.size() == 1);
    for(auto& key : sub_bucket)
    {
        if(gi_inverse.find(key) == gi_inverse.end())
            gi_inverse.insert({key, {}});
    }
    int iii = 1;
    for(const auto& match : incremental_matches)
    {
        auto key = n + iii;
        gi_map.insert({key, match}); iii++;
        for(const auto& item : match)
        {
            auto search = gi_inverse.find(item);
            if(search != gi_inverse.end())
            {
                search->second.push_back(key);
            }
        }
    }
    REQUIRE(gi_map.size() == n);
    print_matches(gi_map);
    lambda_print(gi_inverse);
}
