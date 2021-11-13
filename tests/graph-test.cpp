#include <iostream>
#include <type_traits>

#include "catch.hpp"
#include "../src/SimpleGraph.hpp"

TEST_CASE("Base Graph Host Test", "[graph host test]")
{    
    std::cout << "Running the Simple Graph Test...\n";
    
    using GraphType = Cajete::SimpleGraph<std::size_t, std::size_t>;
    
    GraphType::index_type node_capacity = 10;
    GraphType::index_type avg_degree = 4; //sets the edge block default size
    GraphType test_graph(node_capacity, avg_degree);
   
    //Initializing Test
    REQUIRE(test_graph.node_capacity == node_capacity);
    REQUIRE(test_graph.edge_block_size == avg_degree);
    REQUIRE(test_graph.edge_capacity == node_capacity*avg_degree);
    REQUIRE(test_graph.num_nodes == 0);
    
    std::cout << "Type size: " << sizeof(GraphType::mapping_type) << "\n";
    GraphType::node_type keys[] = {0, 3, 2, 1, 4};
    
    for(auto i = 0; i < 5; i++) {
        test_graph.insert_node(keys[i]);
    }
   
    //Node Insertion Test
    REQUIRE(test_graph.num_nodes == 5);
    //REQUIRE(test_graph.unique_id == 5);

    for(auto i = 2; i < 5; i++) {
        test_graph.insert_edge(0, i);
    }
    
    for(auto i = 1; i < 3; i++) {
        test_graph.insert_edge(4, i);
    }
    
    for(auto i = 2; i < 5; i++) {
        test_graph.insert_edge(0, i);
    }
    
    auto inserted_edges = 8;
    REQUIRE(test_graph.num_edges == inserted_edges);

    for(auto i = 0; i < test_graph.edge_capacity; i++) {
        auto src = test_graph.edge_set[i].out_edge;
        auto loc = test_graph.edge_set[i].location;
        auto dst = test_graph.node_set[loc];
        if( !(i % 5) )
            std::cout << "\n";
        std::cout << "[Key:" << src << ", Loc: " << loc << "] ";

    } std::cout << std::endl;
    
    std::cout << test_graph;

    REQUIRE(test_graph.out_degree(4) == 2);
    REQUIRE(test_graph.isolated(4) == false);
    REQUIRE(test_graph.isolated(1) == true);
    
    for(auto i = 0; i <5; i++) {
        std::cout << test_graph.node_set[i] << " " << test_graph.out_degree(i) << " " << test_graph.isolated(i) << "\n";
    }
    std::cout << test_graph.out_neighbors_begin(0) << " " << test_graph.out_neighbors_end(0) << "\n";
}
