#ifndef COMPONENT_MAP_HPP
#define COMPONENT_MAP_HPP

#include <vector> 
#include <algorithm> 
#include <utility>
#include <iostream>
#include <map>

#include "YAGL_Algorithms.hpp"
#include "YAGL_Graph.hpp"
#include "AnalyzedGrammar.hpp"
#include "KeyGenerator.hpp"

namespace DGGML
{
    //forward declaration
    template <typename GraphType>
    struct CellList;

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
    struct ComponentMap
    {
        using key_type = KeyType;
        using data_type = Instance<KeyType>;
        using gen_type = KeyGenerator<KeyType>;
        using container_type = std::unordered_map<key_type, data_type>;
        using iterator = typename container_type::iterator;    

        //the component matches
        container_type matches;

        //generates unique keys for the new rules 
        gen_type key_gen;
        
        //default constructor
        
        auto insert(data_type match)
        {
            //generate a new key and add the match 
            auto k = key_gen.get_key();
            return matches.insert({k, match});
        }

        void erase(key_type k)
        {
            matches.erase(k);
        }

        iterator find(key_type k)
        {
            return matches.find(k);
        }

        data_type& operator[](key_type k)
        {
            return matches[k];
        }
        
        //perhaps these should be over the index not the match container?
        iterator begin() { return matches.begin(); }
        iterator end() { return matches.end(); }

        auto size() const { return matches.size(); }
        
        //TODO could be a template?
        auto count(std::size_t r)
        {
            auto total = 0;
            for(const auto& m : matches)
            {
                if(m.second.type == r)
                    total++;
            }
            return total;
        }

        void print_index()
        {
            for(const auto& [k, v] : matches)
            {
                std::cout << "{ Key: " << k << ", Value: " 
                    << v << " }\n";
                    //<< matches[v].first << " }\n";               
            }
        };
    };
} //end namespace DGGML

#endif 


