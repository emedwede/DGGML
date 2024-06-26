
#ifndef DGGML_INCREMENTALUPDATE_HPP
#define DGGML_INCREMENTALUPDATE_HPP

#include <vector>
#include <algorithm>
#include <utility>
#include <iostream>
#include <map>

#include "ComponentMatchMap.hpp"
#include "YAGL_Algorithms.hpp"
#include "YAGL_Graph.hpp"
#include "AnalyzedGrammar.hpp"
#include "CellList.hpp" //Would cause a circular dependency
#include "pattern_matching.hpp"
#include "phi_functions.hpp"
#include "HelperStructs.hpp"

namespace DGGML
{
    struct Invalidations
    {
        std::set<std::size_t> component_invalidations;
        std::set<std::size_t> rule_invalidations;

        void print()
        {
            std::cout << "Invalid components: { ";
            for(auto& item : component_invalidations) std::cout << item << " ";
            std::cout << "}\n";
            std::cout << "Invalid rules: { ";
            for(auto& item : rule_invalidations) std::cout << item << " ";
            std::cout << "}\n";
        }
    };

    struct RewriteUpdates
    {
        std::set<std::size_t> node_removals, node_updates;
        using cmp = struct PairComparator {
            using pair_t = std::pair<std::size_t,std::size_t>;
            bool operator()(const pair_t& a, const pair_t& b) const {
                int min_a = std::min(a.first, a.second);
                int max_a = std::max(a.first, a.second);

                int min_b = std::min(b.first, b.second);
                int max_b = std::max(b.first, b.second);

                if (min_a < min_b) {
                    return true;
                } else if (min_a > min_b) {
                    return false;
                } else {
                    return max_a < max_b;
                }
            }
        };
        std::set<std::pair<std::size_t, std::size_t>, PairComparator> edge_removals, edge_updates;

        void print()
        {
            std::cout << "Nodes to invalidate: { ";
            for(auto& item : node_removals) std::cout << item << " ";
            std::cout << "}\n";
            std::cout << "Edges to invalidate: { ";
            for(auto& item : edge_removals) std::cout << "( " << item.first << ", " << item.second  << " ) ";
            std::cout << "}\n";
            std::cout << "Nodes to validate: { ";
            for(auto& item : node_updates) std::cout << item << " ";
            std::cout << "}\n";
            std::cout << "Edges to validate: { ";
            for(auto& item : edge_updates) std::cout << "( " << item.first << ", " << item.second  << " ) ";
            std::cout << "}\n";
        }
    };

    template <typename GraphType, typename MatchType>
    std::pair<RewriteUpdates, Invalidations> perform_invalidations_and_rewrite(DGGML::RuleMatch<std::size_t>& inst,
                                                     MatchType& component_match_set,
                                                     KeyGenerator<std::size_t>& gen,
                                                     DGGML::AnalyzedGrammar<GraphType>& grammar_analysis,
                                                     GraphType& system_graph,
                                                     RuleMatchMap<std::size_t>& rule_matches,
                                                     std::vector<std::size_t>& rule_map,
                                                     CellList<GraphType>& cell_list,
                                                     GeocellProperties<std::size_t>& geocell_properties)
    {

        RewriteUpdates changes;
        std::string rname = inst.name;

        //TODO: double check this mapping code is actually correct, see changes added in the previous commit before TODO
        //construct a vertex map for the lhs to the rule instance
        std::map<std::size_t, std::size_t> lhs_vertex_map;
        std::vector<std::size_t> left, mid, right;
        for(auto& m : grammar_analysis.ccuv_mappings[rname])
        {
            for (auto &[k, v]: m)
            {
                left.push_back(k);
                mid.push_back(v);
            }
        }

        for(auto& c : inst.components)
            for(auto& k : component_match_set[c].match)
                right.push_back(k);

        //print out the mapping info, so we know how a lhs numbering maps to an instance numbering
//        std::cout << "\nMappings: { LHS Key -> Minimal Component Key -> Rule ComponentMatch Key }\n";
        for(auto i = 0; i < left.size(); i++)
        {
            lhs_vertex_map[left[i]] = right[i];
//            std::cout << "{ " << left[i] << " -> " << mid[i] << " -> " << right[i] << " }\n";
        }

        auto& rewrite = grammar_analysis.with_rewrites.at(rname);
//        rewrite.print_node_sets(rname);
//        std::cout << "\n";
//        rewrite.print_edge_sets(rname);

        //TODO: Induced subgraph may include edges we didn't expect and create a rewrite error
        auto lhs_match = YAGL::induced_subgraph(system_graph, right);
        //TODO: potentially fix by searching the induced graph and removing edges not in the match
        auto lhs_match_copy = lhs_match;
        auto rhs_rule_copy = grammar_analysis.with_rules.at(rname).rhs_graph;

        //we can build rhs of the vertex map by making a copy of the lhs and deleting
        auto rhs_vertex_map = lhs_vertex_map;

        for(auto& k : rewrite.node_set_destroy)
        {
            rhs_vertex_map.erase(k);
            lhs_match_copy.removeNode(lhs_match_copy.findNode(lhs_vertex_map[k])->second);
            changes.node_removals.insert(lhs_vertex_map[k]);
        }

        for(auto& k : rewrite.node_set_create)
        {
            //doing this way helps cheese my way into creating the correct types
            auto n = rhs_rule_copy.findNode(k)->second;
            decltype(n) node(gen.get_key(), n.getData());
            lhs_match_copy.addNode(node);
            rhs_vertex_map[k] = node.getKey();
            changes.node_updates.insert(rhs_vertex_map[k]);
        }

        for(auto& [u, v] : rewrite.edge_set_destroy)
        {
            //need to use the removal list edges and map those to the correct edges for the match
            lhs_match_copy.removeEdge(lhs_vertex_map[u], lhs_vertex_map[v]);
            changes.edge_removals.insert({lhs_vertex_map[u], lhs_vertex_map[v]});
        }

        // I think I need to have a vertex mapping for the rhs, which is incomplete until
        // all new nodes are created since their keys are uniquely generated
        for(auto& [u, v] : rewrite.edge_set_create) {
            lhs_match_copy.addEdge(rhs_vertex_map[u], rhs_vertex_map[v]);
            changes.edge_updates.insert({rhs_vertex_map[u], rhs_vertex_map[v]});
        }
        //TODO: I also need to do some work to make sure any nodes that change only type
        // are changed. This is because the node rewrite set currently finds set difference
        // on keys not types, below is a first attempt, but there may be a more efficient way
        for(auto& [k, _] : rhs_rule_copy.getNodeSetRef())
        {
            auto& v1 = lhs_match_copy[rhs_vertex_map[k]];

            auto& v2 = rhs_rule_copy[k];

            if(v1.type != v2.type)
            {
                v1.setData(v2.data);
                //need to add both since rules containing it are invalid, and it's also a candidate since it's new
                changes.node_removals.insert(rhs_vertex_map[k]);
                changes.node_updates.insert(rhs_vertex_map[k]);
                //copy assignment doesn't work since it won't invoke updating type
                //v1.data = v2.data;
            }
        }

        //I think we actually need a map for the lhs, and the rhs
        grammar_analysis.with_rules.at(rname).update(lhs_match, lhs_match_copy, lhs_vertex_map, rhs_vertex_map); //h

//        std::cout << "currently mapped rules: { ";
//        for(auto& key : rule_map) std::cout << key << " "; std::cout << "}\n";

        //TODO: need a better fix, currently here because the cell list uses the graph to find the anchor
        // for a component so that we can find its cell, but if the node has already been removed by a rewrite
        // there is nothing to find leading to an error, potential alternate fix, add an inverse mapping
        // of components to cells
        auto removals = perform_invalidations(changes, component_match_set, grammar_analysis,
                                              rule_matches, rule_map, cell_list, geocell_properties);
        //update the system graph
        for(auto& k : rewrite.node_set_destroy)
            system_graph.removeNode(system_graph.findNode(lhs_vertex_map[k])->second);

        for(auto& [u, v] : rewrite.edge_set_destroy)
            system_graph.removeEdge(lhs_vertex_map[u], lhs_vertex_map[v]);

        for(auto& [k, v] : lhs_match_copy.getNodeSetRef())
        {
            auto res = system_graph.findNode(k);
            if(res != system_graph.node_list_end())
                res->second = v; //just update the data
            else
                system_graph.addNode(v);
        }

        for(auto& [u, v] : rewrite.edge_set_create)
            system_graph.addEdge(rhs_vertex_map[u], rhs_vertex_map[v]);

        return std::make_pair(changes, removals);
    }

    //TODO: maybe this function needs to return even more info on the rules it invalidated
    //TODO: we could accelerate the invalidations by using the cell list to prune the search space
    // or we could create an inverse relation between components and rules, but pay the price of keeping
    // track of it
    template<typename GraphType>
    Invalidations perform_invalidations(RewriteUpdates& changes,
                                        ComponentMatchMap<std::size_t>& component_matches,
                                        DGGML::AnalyzedGrammar<GraphType>& grammar_analysis,
                                        RuleMatchMap<std::size_t>& rule_matches,
                                        std::vector<std::size_t>& rule_map,
                                        CellList<GraphType>& cell_list, GeocellProperties<std::size_t>& geocell_properties)
    {
        //std::cout << "This function invalidates\n";
        Invalidations removals;
        //we must first search for any components containing the nodes/edges that were invalidated
        //note: worst case we search the whole list, better case we use the cell list to speed up search
        auto& component_invalidations = removals.component_invalidations;

        //search components for nodes and edges
        for(auto& [k1, m] : component_matches)
        {
            for(auto i = 0; i < m.match.size(); i++)
            {
                auto k2 = m.match[i];
                // check for node and if found mark for invalidation
                if(changes.node_removals.find(k2) != changes.node_removals.end())
                {
                    component_invalidations.insert(k1);
                }
                //now check for edges and if found mark for invalidation
                for(auto j = 0; j < m.match.size(); j++)
                {
                    auto k3 = m.match[j];
                    auto& c = grammar_analysis.unique_components[m.type];
                    //match contains the nodes
                    std::vector<std::size_t> ordering = grammar_analysis.orderings[m.type];
                    if(changes.edge_removals.find({k2, k3}) != changes.edge_removals.end()
                       && c.adjacent(ordering[i], ordering[j]))
                    {
                        //std::cout << "{ " << ordering[i] << " -> " << k2 << " }\n";
                        //std::cout << "{ " << ordering[j] << " -> " << k3 << " }\n";
                        //std::cout << "component " << k1 << " contains the edge\n";
                        component_invalidations.insert(k1);
                    }
                }
            }
        }

        //add in the components we invalidated, so we can check lower dimensions for any invalid rules
        //TODO: we really only need to add in ones near the boundary
        for(auto& item : component_invalidations)
        {
            geocell_properties.invalidated_components.insert(item);
        }

        //std::cout << "CellList size before invalidations " << cell_list.getTotalSize() << "\n";
        //std::cout << "ComponentMatchMap size before invalidations " << component_matches.size() << "\n";
        //use list of removable components to delete items from component matches and cell_list
        for(auto& item : component_invalidations)
        {
            //Delete from cell list first otherwise, we would need a copy of the erased component
            //to locate its cell via it's position in the list
            auto iter = component_matches.find(item);
            cell_list.erase(iter);
            component_matches.erase(item);
            //std::cout << "erasing component " << item << "\n";
//            if(auto search = component_matches.find(item); search != component_matches.end())
//                std::cout << "somehow we still found " << item << "\n";
        }
        //std::cout << "CellList size after invalidations " << cell_list.getTotalSize() << "\n";
        //std::cout << "ComponentMatchMap size after invalidations " << component_matches.size() << "\n";

        //if a component is removed is in a boundary cell, it may participate in a rule instance of another dimension,
        //so we could mark it for that as well

        //TODO: instead of a direct search we should search the cell of which a rule instance resides and it's nbrs
        // for any rule instance containing the invalidated components
        auto& rule_invalidations = removals.rule_invalidations;
        for(auto& k : rule_map)
        {
            auto& inst = rule_matches[k];
            for(auto& c : inst.components)
            {
                if(component_invalidations.find(c) != component_invalidations.end())
                {
                    rule_invalidations.insert(k);
                }
            }
        }

        //we remove these now invalid rules
        //std::cout << "Rule matches before erase: " << rule_matches.size() << "\n";
        //std::cout << "Rule map before erase: " << rule_map.size() << "\n";
        {
            int count = 0;
            for(auto& item : rule_map)
                if(rule_matches[item].category == "deterministic")
                    count++;
            //std::cout << "Number of deterministic rules: " << count << "\n";
         }
        for(auto& item : rule_invalidations) {
            rule_matches.erase(item);
            //we also need to remove the invalid rules from the rule map
            rule_map.erase(std::find(rule_map.begin(), rule_map.end(), item));
        }
        //std::cout << "Rule matches after erase: " << rule_matches.size() << "\n";
        //std::cout << "Rule map after erase: " << rule_map.size() << "\n";
        //return the list of boundary components invalidated

        return removals;
    }

    template<typename GraphType, typename CellComplexType>
    void find_new_matches(RewriteUpdates& changes,
                          GraphType& system_graph,
                          ComponentMatchMap<std::size_t>& component_matches,
                          DGGML::AnalyzedGrammar<GraphType>& grammar_analysis,
                          RuleMatchMap<std::size_t>& rule_matches,
                          std::vector<std::size_t>& rule_map,
                          CellList<GraphType>& cell_list,
                          CellComplexType& geoplex2D,
                          std::size_t cell_k,
                          double reaction_radius,
                          GeocellProperties<std::size_t>& geocell_properties)
    {
        //std::cout << "This function finds new matches\n";
        //goal, take the candidate nodes, do a search the depth of the height of the tallest rooted spanning tree
        //for each candidate and create a set of nodes used to induce a graph that we will search for new components
        //What should the candidate nodes be?
        std::set<std::size_t> candidate_nodes = changes.node_updates;
        //candidate nodes are nodes updated/added or and nodes part of an edge pair
        for(auto& edge : changes.edge_updates)
        {
            candidate_nodes.insert(edge.first);
            candidate_nodes.insert(edge.second);
        }
        std::set<std::size_t> inducers;
        for(auto n : candidate_nodes)
        {
            //TODO: depth needs to be as long as the deepest rule pattern
            auto res = YAGL::recursive_dfs(system_graph, n, 5);
            for(auto& item : res)
                inducers.insert(item);
        }

        auto candidate_graph = YAGL::induced_subgraph(system_graph, inducers);
        //std::cout << candidate_graph << "\n";

        //TODO: fuse the acceptance functions
        //search the candidate graph for components, the code below finds and validates in one loop
        std::vector<std::size_t> validated_components;
        //TODO: I think we need to store the ordering or the rooted spanning tree
        for(auto& [k, pattern] : grammar_analysis.unique_components) {
            //need to actually get the ordering for the mapping
            auto matches = YAGL::subgraph_isomorphism2(pattern, candidate_graph);
            //std::cout << "Match sizes found for " << k << " " << matches.size() << "\n";
            if (!matches.empty()) {
                for (auto &m: matches) {
                    bool accepted = false;
                    for (auto& [v1, v2] : m) {
                        for (auto& n: changes.node_updates) {
                            if (v2 == n) {
                                accepted = true;
                                std::vector<typename GraphType::key_type> order;
                                std::vector<typename GraphType::key_type> match;
                                for (auto &[key, value]: m) {
                                    order.push_back(key);
                                    match.push_back(value);
                                }
                                ComponentMatch<typename GraphType::key_type> inst;
                                inst.match = match;
                                inst.order = order;
                                inst.type = k;
                                inst.anchor = match[0];
                                auto res = component_matches.insert(inst);
                                if (res.second)
                                    validated_components.push_back(res.first->first);
                                break;
                            }
                        }
                        if (accepted) break;
                    }
                    if (accepted) continue;
                    for (auto& [v1, v2] : m) {
                        for (auto& [u1, u2] : m) {
                            for (auto& [n1, n2]: changes.edge_updates) {
                                if (v2 == n1 && u2 == n2 && pattern.adjacent(v1, u1)) {
                                    //The accepted code is copy pasted and could be fused
                                    accepted = true;
                                    std::vector<typename GraphType::key_type> order;
                                    std::vector<typename GraphType::key_type> match;
                                    for (auto &[key, value]: m) {
                                        order.push_back(key);
                                        match.push_back(value);
                                    }
                                    ComponentMatch<typename GraphType::key_type> inst;
                                    inst.match = match;
                                    inst.order = order;
                                    inst.type = k;
                                    inst.anchor = match[0];
                                    auto res = component_matches.insert(inst);
                                    if (res.second)
                                        validated_components.push_back(res.first->first);
                                    break;
                                }
                                if (accepted) break;
                            }
                        }
                        if (accepted) break;
                    }
                    if (accepted) continue;
                }
            }
        }
//        std::cout << "Validated components: { ";
//        for(auto& item : validated_components)
//            std::cout << item << " ";
//        std::cout << "}\n";
//        for(auto& item : validated_components)
//        {
//            std::cout <<  component_matches[item].type << ": " << component_matches[item] << "\n";
//        }

        //new components need to be added back into their dynamic cell list
        //std::cout << "CellList size before insertion " << cell_list.getTotalSize() << "\n";
        for(auto& item : validated_components)
        {
            auto iter = component_matches.find(item);
            cell_list.insert_match(iter);
        }
        //std::cout << "CellList size after insertion " << cell_list.getTotalSize() << "\n";

        //after finding components and adding them to the cell list, we need find all new rule instances
        //we only need to search in the cells surrounding a newly inserted component
        //TODO: there is definitely a more efficient way of doing this other than searching the whole neighborhood of cells
        // for every pattern. For example, why search for a pattern if that type of validated component isn't even in it.

        std::vector<RuleMatch<std::size_t>> accepted_rule_matches;
        // We need all the components from the inserted cells and neighbors as starting components for our search
        std::set<std::size_t> candidate_cells;
        for(auto& item : validated_components)
        {
            auto m1 = component_matches.find(item);
            auto c = cell_list.locate_cell(m1);
            int imin, imax, jmin, jmax;
            cell_list.getCells(c, imin, imax, jmin, jmax);
            for (auto i = imin; i < imax; i++) {
                for (auto j = jmin; j < jmax; j++) {
                    auto nbr_idx = cell_list.cardinalCellIndex(i, j);
                    candidate_cells.insert(nbr_idx);
                }
            }
        }
        for(auto& cell : candidate_cells)
        {
            for(auto& item : cell_list.data[cell]) {
                auto m1 = component_matches.find(item);
//                if(m1 == component_matches.end())
//                    std::cout << "did not find " << item << "\n";
                for (auto &[name, pattern]: grammar_analysis.rule_component) {
                    int k = 0;
                    std::vector<std::size_t> result;
                    result.resize(pattern.size());
//                    std::cout << "searching for " << name << "\n";
                    if (pattern.size() && m1->second.type == pattern.front()) {
                        result[k] = m1->first;
                        k++;
                        //using this cell should work since the component was just added and can't drift
                        auto c = cell_list.locate_cell(m1);
//                        std::cout << "match " << m1->first << " is located in cell " << c << "\n";
                        incremental_reaction_instance_backtracker(accepted_rule_matches, validated_components,
                                                                  system_graph, component_matches,
                                                                  grammar_analysis, cell_list, name,
                                                                  result, k, pattern, c, reaction_radius);
                    }
                }
            }
        }
//        std::cout << "Number of accepted reaction instances: " << accepted_rule_matches.size() << "\n";

        //rule instances are only accepted if they contain the new components and phi maps them to the current cell
        //if phi maps them to a different cell, we may be able to add them to a separate future validations
        // list rather than do nothing with them
//        std::cout << "Size of rule instances before phi: " << rule_matches.size() << "\n";
//        std::cout << "Size of rules mapped to the geocell before phi: " << rule_map.size() << "\n";
        for(auto& inst : accepted_rule_matches)
        {
            auto max_cell = min_dim_phi(inst, system_graph, geoplex2D, component_matches);
//            std::cout << "we are in cell " << cell_k << " and accepted instance of type " << inst.name << " maps to cell " << max_cell << "\n";
            //this cell owns the instance so it can add it back in
            if(max_cell == cell_k)
            {
                auto res = rule_matches.insert(inst);
                //if successfully inserted, added to the rule_map
                if(res.second)
                {
                    rule_map.push_back(res.first->first);
                }
            } else
            {
//                std::cout << "####################REJECTED#################\n";
//                std::cout << "reject rule is of type: " << inst.name << "\n";
                geocell_properties.rejected_rule_matches.push_back({cell_k, inst});
                //std::cin.get();
            }
        }
//        std::cout << "Size of rule instances after phi: " << rule_matches.size() << "\n";
//        std::cout << "Size of rules mapped to the geocell after phi: " << rule_map.size() << "\n";
        {
            int count = 0;
            for(auto& item : rule_map)
                if(rule_matches[item].category == "deterministic")
                    count++;
//            std::cout << "Number of deterministic rules: " << count << "\n";
        }
    }
}
#endif //DGGML_INCREMENTALUPDATE_HPP
