#ifndef RULE_SYSTEM_HPP
#define RULE_SYSTEM_HPP

#include <vector> 
#include <algorithm> 
#include <utility>
#include <iostream>
#include <map>

//TODO seperate rule system from plant grammar later 
#include "PlantTypes.hpp"
#include "YAGL_Algorithms.hpp"
#include "YAGL_Graph.hpp"

namespace Cajete 
{
    using key_type = Cajete::Plant::mt_key_type; 
    using data_type = Cajete::Plant::MT_NodeData;

    using graph_type = YAGL::Graph<key_type, data_type>;
    using node_type = YAGL::Node<key_type, data_type>;

    struct LHS
    {
        graph_type graph;
        std::vector<graph_type> components;
        
        LHS() {}

        LHS(graph_type& graph) : graph(graph) 
        {
            //build the list of components 
            std::unordered_set<key_type> visited;
            std::size_t count = 0; //no connected components found to start
            for(auto i = graph.node_list_begin(); i != graph.node_list_end(); i++)
            {
                auto v = i->first;
                //node hasn't been visited so it must be the start of a new connected component 
                if(visited.find(v) == visited.end())
                {
                    std::vector<key_type> path;
                    //we could use whatever search method we feel like 
                    YAGL::impl_iterative_bfs(graph, v, visited, path);
                    components.push_back(YAGL::induced_subgraph(graph, path));
                    count++;
                }
            }
        }
    };

    //should LHS and RHS be the same object?
    struct RHS
    {
        graph_type graph;
        std::vector<graph_type> components;
        
        RHS() {}

        RHS(graph_type& graph) : graph(graph) {}
    };

    struct RuleType
    {
        LHS lhs;
        RHS rhs;
        std::string name;
        
        RuleType() {}

        RuleType(std::string name, graph_type& lhs_graph, graph_type& rhs_graph) : name(name), lhs(lhs_graph), rhs(rhs_graph) {}
    };

    struct Grammar 
    {
        std::map<std::string, RuleType> rule_set;
        std::map<std::string, std::vector<std::size_t>> rule_component;
        std::map<std::size_t, std::vector<std::string>> component_rule;
        std::vector<graph_type> minimal_set;

        void addRule(RuleType& r)
        {
            if(rule_set.find(r.name) == rule_set.end())
            {
                rule_set[r.name] = r;
                incremental_build_min_set(r);
            }
            else
                std::cout << "name already exists in rule_set\n";
        }

        void incremental_build_min_set(RuleType& r)
        {
            
            for(auto& c : r.lhs.components)
            {
                auto idx = 0;
                bool found = false;
                for(auto i = 0; i < minimal_set.size(); i++)
                {
                    auto matches = YAGL::subgraph_isomorphism2(c, minimal_set[i]);
                    if(!matches.empty())
                    {
                        idx = i;
                        found = true;
                        break;
                    }
                }
                if(!found)
                {
                    //add forward and backward references
                    rule_component[r.name].push_back(minimal_set.size());
                    component_rule[minimal_set.size()].push_back(r.name);
                    minimal_set.push_back(c);
                }
                else //found
                {
                    rule_component[r.name].push_back(idx);
                    component_rule[idx].push_back(r.name);
                }
            }
        }

        void print() 
        {
            std::cout << "Rules: ";
            for(auto& rule : rule_set)
                std::cout << "\t" << rule.first << "\n";
        }
    };

    void print_mapping(Grammar& gamma)
    {
        std::cout << "\nRule Component Mapppings:\n";
        for(auto& [key, value] : gamma.rule_set)
        {
            std::cout << "\t" << key << ": ";
            for(auto& item : gamma.rule_component[key])
            {
                std::cout << item << " ";
            } std::cout << "\n";
        }

        std::cout << "\nComponent Rule Mapppings:\n";
        for(auto key = 0; key < gamma.minimal_set.size(); key++)
        {
            std::cout << "\t" << key << ": ";
            for(auto& item : gamma.component_rule[key])
            {
                std::cout << item << " ";
            } std::cout << "\n";
        }

    }


    //TODO: Perhaps make this class a member of the rule system 
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

        KeyGenerator() : current_key(0) {}

        key_type get_key() { return current_key++; }
    };

    template <typename KeyType>
    struct Instance 
    {
        using key_type = KeyType;
        using instance_type = std::vector<key_type>;
        instance_type match;
        
        Rule type;
        key_type anchor;

        using iterator = typename std::vector<key_type>::iterator;

        Instance() {}

        Instance(instance_type _match, Rule _type) 
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
            return;
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

        data_type& operator[](key_type k)
        {
            return matches[index[k]].first;
        }
        
        //perhaps these should be over the index not the match container?
        iterator begin() { return matches.begin(); }
        iterator end() { return matches.end(); }

        const auto size() const { return matches.size(); }
        
        //TODO could be a template?
        auto count(enum Rule r)
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
} //end namespace cajete 

#endif 


