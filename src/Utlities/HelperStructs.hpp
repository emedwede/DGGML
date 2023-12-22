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
        std::vector<CellMatchPair<NodeKeyType>> rejected_rule_matches;
    };

    template<typename GeoplexKeyType, typename NodeKeyType>
    using GeocellPropertiesList = std::map<GeoplexKeyType, GeocellProperties<NodeKeyType>>;

}
#endif //DGGML_HELPERSTRUCTS_HPP
