#ifndef RULE_SYSTEM_HPP
#define RULE_SYSTEM_HPP

#include <vector> 
#include <algorithm> 
#include <utility>
#include <iostream>
#include <map>

//TODO seperate rule system from plant grammar later
#include "YAGL_Algorithms.hpp"
#include "YAGL_Graph.hpp"
#include "AnalyzedGrammar.hpp"

namespace DGGML
{
    //TODO: Please remove after the rework
    enum class Rule {G, R, I};

    std::ostream& operator << (std::ostream& os, const Rule& r)
    {
        if(Rule::G == r)
            os << "Growing";
        else if(Rule::R == r)
            os << "Retraction";
        else
            os << "Intermediate";
        return os;
    }

    //We need matching rewrite functions for every enum class member

    template<typename KeyType>
    struct KeyGenerator 
    {
        using key_type = KeyType;
        
        key_type current_key;

        KeyGenerator(KeyType key) : current_key(key) {}

        KeyGenerator() : current_key(0) {}

        key_type get_key() { return current_key++; }
    };

    template <typename KeyType>
    struct RuleInstType {
        std::string category;
        std::string name;
        std::vector<std::size_t> components;
        KeyType anchor;
    };

    template <typename KeyType>
    struct Instance 
    {
        using key_type = KeyType;
        using instance_type = std::vector<key_type>;
        instance_type match;
        
        std::size_t type;
        key_type anchor;

        using iterator = typename std::vector<key_type>::iterator;

        Instance() {}

        Instance(instance_type _match, std::size_t _type)
            : match(_match), type(_type)
        {
            anchor = match[0];    
        }

        iterator begin() { return match.begin(); }
        iterator end() { return match.end(); }

        key_type& operator[](unsigned int i) { return match[i]; }

        const auto size() { return match.size(); }
    };
    
    template<typename T>
    std::ostream& operator<<(std::ostream& os, Instance<T>& match)
    {
        os << "{ ";
        for(const auto& item : match)
            os << item << " ";
        os << "}\n";
        return os;
    }

    template <typename KeyType>
    struct RuleSystem 
    {
        using key_type = KeyType;
        using data_type = Instance<KeyType>;
        using pair_type = std::pair<data_type, key_type>;
        using gen_type = KeyGenerator<KeyType>;
        //TODO: in the upgraded class, just use a map for container type
        using container_type = std::vector<std::pair<data_type, key_type>>;
        using index_type = std::unordered_map<key_type, std::size_t>;  
        using inverse_type = std::unordered_map<key_type, std::vector<key_type>>;
        using cell_list_type = std::unordered_map<key_type, key_type>;
        using iterator = typename container_type::iterator;    

        //containes the matches and rule id bundled as a pair
        container_type matches;
        
        //indexes the container data 
        index_type index;
        
        //inversly maps the rules a node belongs too
        inverse_type inverse_index; 

        //generates unique keys for the new rules 
        gen_type key_gen;
        
        //default constructor
        
        void push_back(data_type match)
        {
            //generate a new key and add the match 
            auto k = key_gen.get_key();
            matches.push_back({match, k});
            
            //build the index
            index.insert({k, matches.size()-1});
            
            //build the reverse index 
            for(const auto& item : match)
            {
                auto search = inverse_index.find(item);

                if(search != inverse_index.end())
                    search->second.push_back(k);
                else 
                    inverse_index[item] = {k};
            }
        }

        //TODO: Incorrect since I think it invalidates all matches a node participates in, but it doesn't
        // remove those rule keys from the inverse index set for the other nodes that participate
        void invalidate(key_type k)
        {
            const auto& rules = inverse_index.find(k);
            if(rules != inverse_index.end())
            {
                for(const auto& item : rules->second)
                {
                    if(index.find(item) != index.end())
                    {
                        //currently uses a back swap method for the matches,
                        //but could be upgraded to allow key reclamation?
                        auto a = index[item];
                        auto b = matches.size()-1;
                        std::swap(matches[a], matches[b]);
                        std::swap(index[matches[a].second], index[matches[b].second]);
                        matches.pop_back();
                        index.erase(item);
                    }
                }
                inverse_index.erase(k);
            }
        }
        
        //TODO: also may not be working correctly
        void invalidate_rule(key_type k)
        {
            auto match = matches[index[k]].first; 
            for(auto& item : match)
            {
                auto& node_rules = inverse_index.find(item)->second;
                for(auto i = 0; i < node_rules.size(); i++)
                {
                    if(node_rules[i] == k)
                    {
                        std::swap(node_rules[i], node_rules[node_rules.size()-1]);
                        node_rules.pop_back();
                        break;
                    }
                }
            }
            auto a = index[k];
            auto b = matches.size()-1;
            std::swap(matches[a], matches[b]);
            std::swap(index[matches[a].second], index[matches[b].second]);
            matches.pop_back();
            index.erase(k);

            //TODO: first problem is we are looking up an inverse index of a rule not a node key
            const auto& rules = inverse_index.find(k);
            std::cout << "here\n";
            //TODO: see above, this check is bugged and only works if rule key is the same a node key
            if(rules != inverse_index.end())
            {
                for(const auto& item : rules->second)
                {
                    if(index.find(item) != index.end())
                    {
                        //currently uses a back swap method for the matches,
                        //but could be upgraded to allow key reclamation?
                        auto a = index[item];
                        auto b = matches.size()-1;
                        std::swap(matches[a], matches[b]);  
                        std::swap(index[matches[a].second], index[matches[b].second]);
                        matches.pop_back();
                        index.erase(item);
                    }
                }
                std::cout << "here\n";
                inverse_index.erase(k);
            }
        }

        data_type& operator[](key_type k)
        {
            return matches[index[k]].first;
        }
        
        //perhaps these should be over the index not the match container?
        iterator begin() { return matches.begin(); }
        iterator end() { return matches.end(); }

        const auto size() const { return matches.size(); }
        
        //TODO could be a template?
        auto count(std::size_t r)
        {
            auto total = 0;
            for(auto it = matches.begin(); it != matches.end(); it++)
            {
                if(it->first.type == r)
                    total++;
            }
            return total;
        }

        void print_index()
        {
            for(const auto& [k, v] : index)
            {
                std::cout << "{ Key: " << k << ", Value: " 
                    << v << " }\n";
                    //<< matches[v].first << " }\n";               
            }
        };
    };

    struct RewriteUpdates
    {
        std::set<std::size_t> node_removals, node_updates;
        using cmp = struct PairComparator {
            using pair_t = std::pair<std::size_t,std::size_t>;
            bool operator()(const pair_t& a, const pair_t& b) const {
                int min_a = std::min(a.first, a.second);
                int max_a = std::max(a.first, a.second);

                int min_b = std::min(b.first, b.second);
                int max_b = std::max(b.first, b.second);

                if (min_a < min_b) {
                    return true;
                } else if (min_a > min_b) {
                    return false;
                } else {
                    return max_a < max_b;
                }
            }
        };
        std::set<std::pair<std::size_t, std::size_t>, PairComparator> edge_removals, edge_updates;

        void print()
        {
            std::cout << "Nodes to invalidate: { ";
            for(auto& item : node_removals) std::cout << item << " ";
            std::cout << "}\n";
            std::cout << "Edges to invalidate: { ";
            for(auto& item : edge_removals) std::cout << "( " << item.first << ", " << item.second  << " ) ";
            std::cout << "}\n";
            std::cout << "Nodes to validate: { ";
            for(auto& item : node_updates) std::cout << item << " ";
            std::cout << "}\n";
            std::cout << "Edges to validate: { ";
            for(auto& item : edge_updates) std::cout << "( " << item.first << ", " << item.second  << " ) ";
            std::cout << "}\n";
        }
    };

    template <typename GraphType>
    RewriteUpdates perform_rewrite(DGGML::RuleInstType<std::size_t>& inst,
                         std::map<std::size_t, DGGML::Instance<std::size_t>>& component_match_set,
                         KeyGenerator<std::size_t>& gen,
                         DGGML::AnalyzedGrammar<GraphType>& gamma_analysis, GraphType& system_graph)
    {

        RewriteUpdates changes;
        std::string rname = inst.name;

        //construct a vertex map for the lhs to the rule instance
        std::map<std::size_t, std::size_t> lhs_vertex_map;
        std::vector<std::size_t> left, mid, right;
        for(auto& m : gamma_analysis.ccuv_mappings[rname])
        {
            for (auto &[k, v]: m)
            {
                left.push_back(k);
                mid.push_back(v);
            }
        }

        for(auto& c : inst.components)
            for(auto& k : component_match_set[c].match)
                right.push_back(k);

        //print out the mapping info, so we know how a lhs numbering maps to an instance numbering
        std::cout << "\nMappings: { LHS Key -> Minimal Component Key -> Rule Instance Key }\n";
        for(auto i = 0; i < left.size(); i++)
        {
            lhs_vertex_map[left[i]] = right[i];
            std::cout << "{ " << left[i] << " -> " << mid[i] << " -> " << right[i] << " }\n";
        }

        auto& rewrite = gamma_analysis.with_rewrites.at(rname);
        rewrite.print_node_sets(rname);
        std::cout << "\n";
        rewrite.print_edge_sets(rname);

        auto lhs_match = YAGL::induced_subgraph(system_graph, right);
        auto lhs_match_copy = lhs_match;
        auto rhs_rule_copy = gamma_analysis.with_rules.at(rname).rhs_graph;

        //we can build rhs of the vertex map by making a copy of the lhs and deleting
        auto rhs_vertex_map = lhs_vertex_map;

        for(auto& k : rewrite.node_set_destroy)
        {
            rhs_vertex_map.erase(k);
            lhs_match_copy.removeNode(lhs_match_copy.findNode(lhs_vertex_map[k])->second);
            changes.node_removals.insert(lhs_vertex_map[k]);
        }

        for(auto& k : rewrite.node_set_create)
        {
            //doing this way helps cheese my way into creating the correct types
            auto n = rhs_rule_copy.findNode(k)->second;
            decltype(n) node(gen.get_key(), n.getData());
            lhs_match_copy.addNode(node);
            rhs_vertex_map[k] = node.getKey();
            changes.node_updates.insert(rhs_vertex_map[k]);
        }

        for(auto& [u, v] : rewrite.edge_set_destroy)
        {
            //need to use the removal list edges and map those to the correct edges for the match
            lhs_match_copy.removeEdge(lhs_vertex_map[u], lhs_vertex_map[v]);
            changes.edge_removals.insert({lhs_vertex_map[u], lhs_vertex_map[v]});
        }

        // I think I need to have a vertex mapping for the rhs, which is incomplete until
        // all new nodes are created since their keys are uniquely generated
        for(auto& [u, v] : rewrite.edge_set_create) {
            lhs_match_copy.addEdge(rhs_vertex_map[u], rhs_vertex_map[v]);
            changes.edge_updates.insert({rhs_vertex_map[u], rhs_vertex_map[v]});
        }
        //TODO: I also need to do some work to make sure any nodes that change only type
        // are changed. This is because the node rewrite set currently finds set difference
        // on keys not types, below is a first attempt, but there may be a more efficient way
        for(auto& [k, _] : rhs_rule_copy.getNodeSetRef())
        {
            auto& v1 = lhs_match_copy[rhs_vertex_map[k]];

            auto& v2 = rhs_rule_copy[k];

            if(v1.type != v2.type)
            {
                v1.setData(v2.data);
                //need to add both since rules containing it are invalid, and it's also a candidate since it's new
                changes.node_removals.insert(rhs_vertex_map[k]);
                changes.node_updates.insert(rhs_vertex_map[k]);
                //copy assignment doesn't work since it won't invoke updating type
                //v1.data = v2.data;
            }
        }

        //I think we actually need a map for the lhs, and the rhs
        gamma_analysis.with_rules.at(rname).update(lhs_match, lhs_match_copy, lhs_vertex_map, rhs_vertex_map); //h

        //update the system graph
        for(auto& k : rewrite.node_set_destroy)
            system_graph.removeNode(system_graph.findNode(lhs_vertex_map[k])->second);

        for(auto& [u, v] : rewrite.edge_set_destroy)
            system_graph.removeEdge(lhs_vertex_map[u], lhs_vertex_map[v]);

        for(auto& [k, v] : lhs_match_copy.getNodeSetRef())
        {
            auto res = system_graph.findNode(k);
            if(res != system_graph.node_list_end())
                res->second = v; //just update the data
            else
                system_graph.addNode(v);
        }

        for(auto& [u, v] : rewrite.edge_set_create)
            system_graph.addEdge(rhs_vertex_map[u], rhs_vertex_map[v]);

        return changes;
    }
} //end namespace DGGML

#endif 


