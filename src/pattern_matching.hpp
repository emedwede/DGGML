#ifndef DGGML_PATTERN_MATCHING_HPP
#define DGGML_PATTERN_MATCHING_HPP
#include <string>
#include <vector>

#include "ComponentMatchMap.hpp"
#include "RuleMatchMap.hpp"

//TODO: improve organization, both functions are just speacialzed versions of a combinatorial backtracking i.e.
//void backtrack(int a[], int k, data input) {
//    int c[MAXCANDIDATES]; //candidates for next position
//    int nc; //next position candidate count
//    int i; //counter
//
//    if (is_a_solution(a, k, input)) {
//        process_solution(a, k, input);
//    } else {
//        k = k + 1;
//        construct_candidates(a, k, input, c, &nc);
//        for (i = 0; i < nc; i++) {
//            a[k] = c[i];
//            make_move(a, k, input);
//            backtrack(a, k, input);
//            unmake_move(a, k, input);
//            if (finished) {
//                return;
//            }
//        }
//    }
//}
/* terminate early */
namespace DGGML {
    // This is a recursive backtracking function, and it is memory efficient because it does a DFS (inorder traversal)
    // so only one vector is need for finding the result
    // TODO: Check for bugs patterns > 2, currently some permutations may not be picked up correctly
    template<typename SimType>
    void reaction_instance_backtracker(SimType& sim, std::string name, std::vector<std::size_t> &result, int k,
                                       std::vector<std::size_t> &pattern, std::size_t c) {
        //if we reach the end of the pattern it's a solution
        if (k == pattern.size()) {
            //a fail condition if components share any nodes
            bool pass = true;
            std::set<std::size_t> all_diff_check;
            for(auto& c1 : result)
            {
                for(auto& c2 : sim.component_matches[c1].match)
                {
                    auto res = all_diff_check.insert(c2);
                    if(!res.second)
                    {
                        pass = false;
                        break;
                    }
                }
                if(!pass) break;
            }
            if(pass) {
                RuleMatch<std::size_t> inst;
                inst.name = name;
                inst.category = "stochastic";
                inst.components = result;
                sim.rule_matches.insert(inst);
            }
        } else {
            //we search in all nearby cells as candidates for the next element in the pattern
            int imin, imax, jmin, jmax;
            sim.cell_list.getCells(c, imin, imax, jmin, jmax);
            for (auto i = imin; i < imax; i++) {
                for (auto j = jmin; j < jmax; j++) {
                    auto nbr_idx = sim.cell_list.cardinalCellIndex(i, j);
                    for (const auto &m2: sim.cell_list.data[nbr_idx]) {
                        bool found = false;
                        for (auto v = 0; v < k; v++) {
                            if (result[v] == m2) found = true;
                        }
                        if (!found && sim.component_matches[m2].type == pattern[k]) {

                            auto a1 = sim.component_matches[result[0]].anchor;
                            auto a2 = sim.component_matches[m2].anchor;
                            auto &p1 = sim.model->system_graph.findNode(a1)->second.getData().position;
                            auto &p2 = sim.model->system_graph.findNode(a2)->second.getData().position;
                            auto d = calculate_distance(p1, p2);

                            // TODO: a constraint for nearness is needed, there are other ways to do it
                            //  this is just a placeholder
                            //TODO: we can actually just accept matches outside the boundary to keep our list valid longer
                            // they'll just be near zero with a proper propensity function
                            //if (d < sim.model->settings.MAXIMAL_REACTION_RADIUS) {
                                result[k] = m2;
                                reaction_instance_backtracker(sim, name, result, k + 1, pattern, c);
                           // }
                        }
                    }
                }
            }
        }
    };

    template<typename GraphType>
    void incremental_reaction_instance_backtracker(std::vector<RuleMatch<std::size_t>>& accepted_rule_matches,
                                                   std::vector<std::size_t>& validated_components,
                                                   GraphType& system_graph,
                                                   ComponentMatchMap<std::size_t>& component_matches,
                                                   DGGML::AnalyzedGrammar<GraphType>& grammar_analysis,
                                                   CellList<GraphType>& cell_list,
                                                   std::string name, std::vector<std::size_t> &result,
                                                   int k, std::vector<std::size_t> &pattern, std::size_t c,
                                                   double reaction_radius) {
        //if we reach the end of the pattern it's a candidate solution
        if (k == pattern.size()) {
            //check if the pattern contains any of the new components
            bool valid = false;
            for(auto& item : validated_components)
            {
                if(std::find(result.begin(), result.end(), item) != result.end())
                {
                    valid = true;
                    break;
                }
            }
            //a fail condition if components share any nodes
            bool pass = true;
            if(valid) {
                std::set<std::size_t> all_diff_check;
                for (auto &c1: result) {
                    for (auto &c2: component_matches[c1].match) {
                        auto res = all_diff_check.insert(c2);
                        if (!res.second) {
                            pass = false;
                            break;
                        }
                    }
                    if (!pass) break;
                }
            }
            if(valid && pass) {
                std::cout << "found: "; for(auto& item : result) std::cout << item << " "; std::cout << "\n";
                RuleMatch<std::size_t> inst;
                inst.name = name;
                inst.category = "stochastic";
                inst.components = result;
                inst.anchor = component_matches[result[0]].anchor;
                accepted_rule_matches.push_back(inst);
            }
        } else {
            //we search in all nearby cells as candidates for the next element in the pattern
            int imin, imax, jmin, jmax;
            cell_list.getCells(c, imin, imax, jmin, jmax);
            for (auto i = imin; i < imax; i++) {
                for (auto j = jmin; j < jmax; j++) {
                    auto nbr_idx = cell_list.cardinalCellIndex(i, j);
                    for (const auto &m2: cell_list.data[nbr_idx]) {
                        bool found = false;
                        for (auto v = 0; v < k; v++) {
                            if (result[v] == m2) found = true;
                        }
                        if (!found && component_matches[m2].type == pattern[k]) {

                            auto a1 = component_matches[result[0]].anchor;
                            auto a2 = component_matches[m2].anchor;
                            auto &p1 = system_graph.findNode(a1)->second.getData().position;
                            auto &p2 = system_graph.findNode(a2)->second.getData().position;
                            auto d = calculate_distance(p1, p2);

                            // TODO: a constraint for nearness is needed, there are other ways to do it
                            //  this is just a placeholder
                            if (d < reaction_radius) {
                                result[k] = m2;
                                incremental_reaction_instance_backtracker(accepted_rule_matches, validated_components,
                                                                          system_graph,component_matches,
                                                                          grammar_analysis, cell_list, name,
                                                                          result, k + 1, pattern, c, reaction_radius);
                            }
                        }
                    }
                }
            }
        }
    };

}
#endif //DGGML_PATTERN_MATCHING_HPP
