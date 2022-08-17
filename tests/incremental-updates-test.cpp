#include <iostream>

#include "catch.hpp"

#include "YAGL_Graph.hpp"
#include "YAGL_Algorithms.hpp"
#include "PlantTypes.hpp"
#include "PlantGrammar.hpp"
#include "PlantUtils.hpp"
#include <vector> 
#include <random> 

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


TEST_CASE("Incremental Update Test", "[update-test]")
{
    graph_type graph;
   
    std::size_t n = 10; // num mt
    microtubule_scatter(graph, n);

    std::vector<Cajete::Plant::mt_key_type> bucket;
    
    for(const auto& [key, value] : graph.getNodeSetRef())
        bucket.push_back(key);
    
    auto matches = Cajete::Plant::microtubule_growing_end_matcher(graph, bucket);
    print_matches(matches);
    REQUIRE(matches.size() == n);
    
}
