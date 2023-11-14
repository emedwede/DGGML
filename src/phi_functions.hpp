#ifndef DGGML_PHI_FUNCTIONS_HPP
#define DGGML_PHI_FUNCTIONS_HPP

#include <cmath>

template<typename RuleInstance, typename GraphType, typename CellComplexType>
std::size_t anchored_phi(RuleInstance& inst, GraphType& system_graph, CellComplexType& geoplex2D)
{
    auto& node_data = system_graph.findNode(inst.anchor)->second.getData();

    double xp = node_data.position[0];
    double yp = node_data.position[1];
    int ic, jc;

    geoplex2D.reaction_grid.locatePoint(xp, yp, ic, jc);
    auto cardinal = geoplex2D.reaction_grid.cardinalCellIndex(ic, jc);
    auto max_cell = geoplex2D.cell_label[cardinal];

    return max_cell;
}

template<typename RuleInstance, typename GraphType, typename CellComplexType, typename ComponentType>
std::size_t min_dim_phi(RuleInstance& inst, GraphType& system_graph, CellComplexType& geoplex2D, ComponentType& component_matches)
{
    std::size_t min_cell;
    std::size_t min_dim;
    for(auto i = 0; i < inst.components.size(); i++)
    {
        auto& m = component_matches[inst.components[i]].match;
        for(auto j = 0; j < m.size(); j++)
        {
            auto& item = m[j];
            auto& node_data = system_graph.findNode(item)->second.getData();

            double xp = node_data.position[0];
            double yp = node_data.position[1];
            int ic, jc;
            geoplex2D.reaction_grid.locatePoint(xp, yp, ic, jc);
            auto cardinal = geoplex2D.reaction_grid.cardinalCellIndex(ic, jc);
            if(i == 0 && j == 0)
            {
                //TODO: dim label are incorrectly in reverse order in the ECC
                min_dim = geoplex2D.dim_label[cardinal];
                min_cell = geoplex2D.cell_label[cardinal];
            } else
            {
                auto temp_min_dim = geoplex2D.dim_label[cardinal];
                auto temp_min_cell = geoplex2D.cell_label[cardinal];
                //TODO: we need to fix the bug above to change the sign
                if(temp_min_dim > min_dim)
                {
                    min_dim = temp_min_dim;
                    min_cell = temp_min_cell;
                }
            }
        }
    }

    return min_cell;
}
#endif //DGGML_PHI_FUNCTIONS_HPP
