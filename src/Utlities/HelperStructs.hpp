#ifndef DGGML_HELPERSTRUCTS_HPP
#define DGGML_HELPERSTRUCTS_HPP

#include <map>

#include "RuleMatchMap.hpp"

namespace DGGML {

    template <typename NodeKeyType>
    struct CellMatchPair
    {
        std::size_t cell_id;
        RuleMatch<NodeKeyType> rule_match;
    };
    template <typename NodeKeyType>
    struct GeocellProperties {
        double tau;
        double exp_sample;
        std::size_t current_rules_fired;
        std::size_t total_rules_fired;
        std::vector<CellMatchPair<NodeKeyType>> rejected_rule_matches;
        std::vector<std::size_t> invalidated_components;

        GeocellProperties() : tau(0.0), exp_sample(0.0), current_rules_fired(0), total_rules_fired(0) {}

        void print()
        {
            std::cout << "Current tau: " << tau << "\n";
            std::cout << "Current exp sample: " << exp_sample << "\n";
            std::cout << "Rules fired this time step: " << current_rules_fired << "\n";
            std::cout << "Total rules fired so far: " << total_rules_fired << "\n";
            std::cout << "Rule matches rejected: " << rejected_rule_matches.size() << "\n";
            std::cout << "Components invalidated: " << invalidated_components.size() << "\n";
        }
    };

    template<typename GeoplexKeyType, typename NodeKeyType>
    using GeocellPropertiesList = std::map<GeoplexKeyType, GeocellProperties<NodeKeyType>>;

}
#endif //DGGML_HELPERSTRUCTS_HPP
