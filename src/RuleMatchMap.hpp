#ifndef DGGML_RULEMATCHMAP_HPP
#define DGGML_RULEMATCHMAP_HPP

#include <vector>
#include <unordered_map>
#include "KeyGenerator.hpp"

namespace DGGML
{

    template <typename KeyType>
    struct RuleMatch {
        std::string category;
        std::string name;
        std::vector<std::size_t> components;
        KeyType anchor;
    };

    //Could combine ComponentMatchMap and RuleMatchMap into a KeyGeneratingMap class
    //so long as they don't end up diverging too much
    template <typename KeyType>
    struct RuleMatchMap
    {
        using key_type = KeyType;
        using data_type = RuleMatch<KeyType>;
        using gen_type = KeyGenerator<KeyType>;
        using container_type = std::unordered_map<key_type, data_type>;
        using iterator = typename container_type::iterator;

        //the Rule matches
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

        auto empty() const { return matches.empty(); }

        //TODO: should this also reset the key generator?
        void clear() const { matches.clear(); }

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

    //TODO add a function that can count the number of rule instance types
//    auto count(std::size_t r)
//    {
//        auto total = 0;
//        for(const auto& m : matches)
//        {
//            if(m.second.type == r)
//                total++;
//        }
//        return total;
//    }
}
#endif //DGGML_RULEMATCHMAP_HPP
