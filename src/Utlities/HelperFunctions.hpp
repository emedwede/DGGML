#ifndef DGGML_HELPERFUNCTIONS_HPP
#define DGGML_HELPERFUNCTIONS_HPP

#include <string>

template <typename T1, typename T2, typename T3>
auto induce_from_set(T1& inst, T2& component_matches, T3& system_graph)
{
    std::vector<std::size_t> inducers;
    for(auto& c : inst.components)
    {
        for(auto& id : component_matches[c])
            inducers.push_back(id);
    }
    return YAGL::induced_subgraph(system_graph, inducers);
}

template <typename InstType, typename GrammarAnalysisType, typename MapType, typename ComponentSetType>
void construct_grammar_match_map(InstType& inst, GrammarAnalysisType& grammar_analysis, MapType& lhs_vertex_map,
                                 ComponentSetType& component_match_set) {
    //I think this is the fix
    for(auto i = 0; i < grammar_analysis.ccuv_mappings[inst.name].size(); i++)
    {
        auto& m = grammar_analysis.ccuv_mappings[inst.name][i];
        auto& c = inst.components[i];
        auto cid = component_match_set[c].type;
        auto& o = grammar_analysis.orderings[cid];
        std::map<int, int> ro;
        for(int j = 0; j < o.size(); j++)
            ro[o[j]] = j;
        for(auto& [k, v] : m)
            lhs_vertex_map[k] = component_match_set[c].match[ro[v]];
    }
    //original code
//    //construct a vertex map for the lhs to the rule instance
//    std::vector<std::size_t> left, right;
//    for (auto &m: grammar_analysis.ccuv_mappings[inst.name])
//        for (auto &[k, v]: m)
//            left.push_back(k);
//    //og rule keys -> matched component keys -> ordering -> instance
//    for (auto &c: inst.components)
//        for (auto& k : component_match_set[c].match)
//            right.push_back(k);
//
//    for (auto i = 0; i < left.size(); i++)
//        lhs_vertex_map[left[i]] = right[i];
}

#endif //DGGML_HELPERFUNCTIONS_HPP
