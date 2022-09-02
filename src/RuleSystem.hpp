#ifndef RULE_SYSTEM_HPP

#include <vector> 
#include <algorithm> 
#include <utility>
#include <iostream> 

namespace Cajete 
{
    //TODO: Perhaps make this class a member of the rule system 
    enum class Rule {G, R, I};

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
        
        using iterator = typename std::vector<key_type>::iterator;

        Instance() {}

        Instance(instance_type _match, Rule _type) 
            : match(_match), type(_type) 
        {
            
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
        using gen_type = KeyGenerator<KeyType>;
        using container_type = std::vector<std::pair<data_type, key_type>>;
        using index_type = std::unordered_map<key_type, std::size_t>;  
        using inverse_type = std::unordered_map<key_type, std::vector<key_type>>;
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

        data_type& operator[](key_type k)
        {
            return matches[index[k]].first;
        }
        
        //perhaps these should be over the index not the match container?
        iterator begin() { return matches.begin(); }
        iterator end() { return matches.end(); }

        auto size() { return matches.size(); }

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


