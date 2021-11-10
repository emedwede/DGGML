#ifndef CAJETE_SIMPLE_GRAPH_HPP
#define CAJETE_SIMPLE_GRAPH_HPP

#include "MemoryManager.hpp"

namespace Cajete 
{

template <typename NodeType, typename CountingType>
struct SimpleGraph {
    using index_type = CountingType;
    using node_type = NodeType;
    
    // Note: The sizeof for a struct is not always equal to the sum of 
    // sizeof of each individual member. This is because of the padding 
    // added by the compiler to avoid alignment issues. Padding is only 
    // added when a structure member is followed by a member with a 
    // larger size or at the end of the structure
    using mapping_type = struct Mapping {
        bool has_edges;
        std::size_t block_id;
        std::size_t block_size;
        std::size_t spots_used;
    }; // turns out this is 32bytes on my system :D

    SimpleGraph() 
        : num_nodes(0)
        , node_capacity(0)
        , edge_capacity(0)
        , edge_block_size(0)
        , next_free_block(0)
        , unique_id(0)
        , node_set(nullptr)
        , edge_set(nullptr)
        , mapping_set(nullptr) {
        std::cout << "Default constructor\n";
    }
    
    SimpleGraph(index_type c, index_type d) {
        num_nodes = 0;
        node_capacity = c;
        edge_block_size = d; //avg degree
        edge_capacity = c*d;
        next_free_block = 0;

        unique_id = 0;
        node_set = Cajete::MemoryManager::allocate_std<node_type>(c);
        edge_set = Cajete::MemoryManager::allocate_std<node_type>(c*d);
        mapping_set = Cajete::MemoryManager::allocate_std<mapping_type>(c);

        
        //set the mapping set and other defaults
        for(auto i = 0; i < c; i++) {
            node_set[i] = 0;
            mapping_set[i].has_edges = false;
            mapping_set[i].block_size = 0;
            mapping_set[i].spots_used = 0;
            //no need to set the block_id since it will be set 
            //when a block is requested
        }

        //set edge_set defaults
        for(auto i = 0; i < edge_capacity; i++) {
            edge_set[i] = 0;
        }
        std::cout << "Allocated node_set, edge_set, mapping_set\n";
    }

   ~SimpleGraph() {
        if(node_set) {
            Cajete::MemoryManager::deallocate_std(node_set);
           std::cout << "deleted node_set\n";
        } else {
            std::cout << "node_set is a nullptr, cannot delete\n";
        }
        if(edge_set) {
            Cajete::MemoryManager::deallocate_std(edge_set);
            std::cout << "deleted edge_set\n";
        } else {
            std::cout << "edge_set is a nullptr, cannot delete\n";
        }

        if(mapping_set) { 
           Cajete::MemoryManager::deallocate_std(mapping_set);
           std::cout << "deleted mapping_set\n";
        } else {
            std::cout << "mapping_set is a nullptr, cannot delete\n";
        }
    }

   //Leave it to the user to guarantee key uniqueness
   void insert_node(const node_type &node) {
       //In a standard interface, we should try to find this node
       //first, and if it exists, we silently fail, but for our 
       //purpose we don't search. No duplication checking can lead
       //to duplicate keys if the user is not careful

       if(num_nodes < node_capacity) {
           node_set[num_nodes] = node;
           num_nodes++; 
           std::cout << "Inserted [key: " << node << "] at [index: " << (num_nodes - 1) << "]\n";
           //unique_id++;
       } else {
           std::cout << "Failure, out of capacity\n";
       }
   }

   //standard interface, we won't use this though
   //void insert_edge(node_type src, node_type dst) {
        //Try to find the src, if not add it 

       //Try to find the dst, if not add it
   //}

   //assume a graph grammar rule already knows the correct pointers
   //to the keys, i.e. insertion by iterators of sorts
   void insert_edge(index_type src, index_type dst) {
        //in one world we check to see if the nodes already exists
        //and if not, we create them. For simplicity, we leave it to
        //the user to ensure nodes already exists.
        //
        //we also could check to see if this edge pair already exists
        //but, again, we assume the user is smarter, so we put no 
        //restrictions in place
        
       //First check to see if the src node has an edge block reserved
       //by checking the mapping set
       if(!mapping_set[src].has_edges) {
           //request a block
           mapping_set[src].block_id = next_free_block;
           next_free_block += edge_block_size;
           mapping_set[src].block_size = edge_block_size;
           mapping_set[src].has_edges = true;
       }
       
       //Next check to see if the edge block has enough space to add
       //a new edge, if not, we need to expand the block size, and we 
       //should periodically perform garbage collection manually
       if(mapping_set[src].spots_used == mapping_set[src].block_size) {
           //we're going to do a full copy, so we need to see how many 
           //blocks we need
           index_type blocks_needed = (mapping_set[src].block_size / edge_block_size) + 1;
           index_type start_range = next_free_block;
           next_free_block += edge_block_size*blocks_needed;
           //TODO: finish this
       }

       //Insert the edge
       auto id = mapping_set[src].block_id + mapping_set[src].spots_used;
       edge_set[id] = dst; //TODO: need a node_type as well, i.e. data and ptr
       mapping_set[src].spots_used++;
   } 

    index_type num_nodes;
    index_type node_capacity;
    index_type edge_capacity;
    index_type edge_block_size;
    index_type next_free_block;
    index_type unique_id;

    node_type* node_set; //TODO: needs a bit to indicate if it's empty
    node_type* edge_set;
    mapping_type* mapping_set;
};

template<typename T>
void print_node(T &data, typename T::index_type idx) {
    std::cout << "\n------------------------------------\n";
    std::cout << "Key at index " << idx << ": "
        << data.node_set[idx] << "\n";
    std::cout << "------------------------------------\n";

}

} //end namespace Cajete
#endif
