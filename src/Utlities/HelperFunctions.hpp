#ifndef DGGML_HELPERFUNCTIONS_HPP
#define DGGML_HELPERFUNCTIONS_HPP

template<typename T1, typename T2, typename T3>
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

#endif //DGGML_HELPERFUNCTIONS_HPP
