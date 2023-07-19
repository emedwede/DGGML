#ifndef DGGML_SIMPLE_GRAPH_HPP
#define DGGML_SIMPLE_GRAPH_HPP

#include "MemoryManager.hpp"

namespace DGGML 
{

template <typename NodeType, typename CountingType>
struct SimpleGraph {
    using index_type = CountingType;
    using node_type = NodeType;
    
    using edge_type = struct EdgeList {
        index_type out_edge;
        index_type location;
    };
    
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
        , num_edges(0)
        , node_capacity(0)
        , edge_capacity(0)
        , edge_block_size(0)
        , next_free_block(0)
        , node_set(nullptr)
        , edge_set(nullptr)
        , mapping_set(nullptr)
        , undirected(true) //undirected by default
    {
        std::cout << "Default constructor\n";
    }
    
    SimpleGraph(index_type c, index_type d) {
        num_nodes = 0;
        num_edges = 0;
        node_capacity = c;
        edge_block_size = d; //avg degree
        edge_capacity = c*d;
        next_free_block = 0;
        undirected = true;
        node_set = DGGML::MemoryManager::allocate_std<node_type>(c);
        edge_set = DGGML::MemoryManager::allocate_std<edge_type>(c*d);
        mapping_set = DGGML::MemoryManager::allocate_std<mapping_type>(c);

        
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
            edge_set[i].out_edge = 0;
            edge_set[i].location = 0;
        }
        std::cout << "Allocated node_set, edge_set, mapping_set\n";
    }

   ~SimpleGraph() {
        if(node_set) {
            DGGML::MemoryManager::deallocate_std(node_set);
           std::cout << "deleted node_set\n";
        } else {
            std::cout << "node_set is a nullptr, cannot delete\n";
        }
        if(edge_set) {
            DGGML::MemoryManager::deallocate_std(edge_set);
            std::cout << "deleted edge_set\n";
        } else {
            std::cout << "edge_set is a nullptr, cannot delete\n";
        }

        if(mapping_set) { 
           DGGML::MemoryManager::deallocate_std(mapping_set);
           std::cout << "deleted mapping_set\n";
        } else {
            std::cout << "mapping_set is a nullptr, cannot delete\n";
        }
    }
  
   // Note: if the block is used but no edges 
   // The index returned is simply 0 for end 
   // and begin

   //Returns the start id of the edge list
   index_type out_neighbors_begin(index_type src) {
        return mapping_set[src].block_id;
   }

   //Returns the end id of the edge list 
   index_type out_neighbors_end(index_type src) {
       return ( mapping_set[src].block_id + mapping_set[src].spots_used - 1);
   }

   index_type out_degree(index_type src) {
       //we're being risky, so we don't find first
       return mapping_set[src].spots_used;
   }

   //We could have an in degrre as well
   index_type degree(index_type src) {
       //we're being risky so we don't find first 
       return out_degree(src);
   }

   //returns true if node has no edges
   bool isolated(index_type src) {
      return ( out_degree(src) == 0 ); 
   }
   bool is_undirected() const {
       return undirected;
   }

   bool is_directed() const {
       return !undirected;
   }

   void set_undirected() {
       undirected = true;
   }

   //We could implement a find vertex id function
   // iterator find(node_type &src) {}

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
   
   //TODO: implement node removal
   bool remove_node(index_type src) {
       node_set[src] = 0;
       return true;
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
           index_type next_start = next_free_block;
           next_free_block += edge_block_size*blocks_needed;
           //TODO: finish this
           auto old_used = mapping_set[src].spots_used;
           auto old_start = mapping_set[src].block_id;
        
           //Copy over and don't clean up
           for(auto i = 0; i < old_used; i++) {
               edge_set[next_start+i].out_edge = edge_set[old_start+i].out_edge;
               edge_set[next_start+i].location = edge_set[old_start+i].location;
           }
           //set the mapping set new start
           mapping_set[src].block_id = next_start;
           mapping_set[src].block_size = edge_block_size*blocks_needed;
       }

       //Insert the edge
       auto id = mapping_set[src].block_id + mapping_set[src].spots_used;
       edge_set[id].out_edge = node_set[dst];
       edge_set[id].location = dst; 
       mapping_set[src].spots_used++;
       num_edges++;
   }

   //TODO: implement edge removal
   bool remove_edge(index_type src, index_type dst) {
       auto start_range = out_neighbors_begin();
       auto end_range = out_neighbors_end();
       
       bool found = false;

       while(start_range != end_range || found) {

       }
       return true;
   }

   //Determines the number of connected components
   void connectedComponents() 
   {
       //TODO: implement connectedComponents 
   }
    
    //TODO: implement a use for this
    bool undirected; //is the graph undirected

    index_type num_nodes;
    index_type num_edges;

    index_type node_capacity;
    index_type edge_capacity;
    
    index_type edge_block_size;
    index_type next_free_block;
    
    node_type* node_set; //TODO: needs a bit to indicate if it's empty
    edge_type* edge_set;
    mapping_type* mapping_set;
};

template <typename T, typename U>
std::ostream &operator<<(std::ostream &os, const SimpleGraph<T, U> &graph)
{
    std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    os << "The graph has " << graph.num_nodes << " nodes and " << graph.num_edges << " edges";
    std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    
    return os;
}

template<typename T>
void print_node(T &data, typename T::index_type idx) {
    std::cout << "\n------------------------------------\n";
    std::cout << "Key at index " << idx << ": "
        << data.node_set[idx] << "\n";
    std::cout << "------------------------------------\n";

}

} //end namespace DGGML
#endif
