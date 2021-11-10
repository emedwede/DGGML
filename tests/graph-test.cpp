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

    for(auto i = 0; i < test_graph.edge_capacity; i++) {
        std::cout << test_graph.edge_set[i] << " ";
    } std::cout << std::endl;


}
