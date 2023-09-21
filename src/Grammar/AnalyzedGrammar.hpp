#ifndef DGGML_ANALYZEDGRAMMAR_HPP
#define DGGML_ANALYZEDGRAMMAR_HPP

#include "Grammar.h"

namespace DGGML
{
    template<typename GraphType, typename MapType>
    void compute_and_store_connected_components(GraphType& graph, MapType& m) {
        // build the connected components for each side
        std::unordered_set<typename GraphType::key_type> visited;
        std::size_t count = 0; //no connected components found to start
        for (auto i = graph.node_list_begin(); i != graph.node_list_end(); i++) {
            auto v = i->first;
            //node hasn't been visited, so it must be the start of a new connected component
            if (visited.find(v) == visited.end()) {
                std::vector<typename GraphType::key_type> path;
                //we could use whatever search method we feel like
                YAGL::impl_iterative_bfs(graph, v, visited, path);
                auto c = YAGL::induced_subgraph(graph, path);
                m[count] = c;
                count++;
            }
        }
    }

    template <typename GraphType>
    struct Rewrite
    {
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

        using key_type = typename GraphType::key_type;
        std::set<key_type> node_set_left, node_set_right, node_set_create, node_set_destroy;

        std::set<std::pair<std::size_t, std::size_t>, PairComparator>
                edge_set_left, edge_set_right, edge_set_create, edge_set_destroy;

        //Rewrite() = default;

        Rewrite(GraphType& lhs_graph, GraphType& rhs_graph)
        {
            for(auto& node : lhs_graph.getNodeSetRef())
                node_set_left.insert(node.first);
            for(auto& node : rhs_graph.getNodeSetRef())
                node_set_right.insert(node.first);

            std::set_difference(node_set_left.begin(), node_set_left.end(),
                                node_set_right.begin(), node_set_right.end(),
                                std::inserter(node_set_destroy, node_set_destroy.begin()));
            std::set_difference(node_set_right.begin(), node_set_right.end(),
                                node_set_left.begin(), node_set_left.end(),
                                std::inserter(node_set_create, node_set_create.begin()));

            //builds edgeset without inverses for an undirected graph using a custom comparator
            auto build_edge_set = [&](auto& edge_set, auto& g)
            {
                for(auto& node : g.getNodeSetRef())
                {
                    auto& nbrs = g.out_neighbors(node.first);
                    for(auto& item : nbrs)
                    {
                        edge_set.insert({node.first, item});
                    }
                }
            };

            build_edge_set(edge_set_left, lhs_graph);
            build_edge_set(edge_set_right, rhs_graph);

            //finds set differences between edge sets, using the custom comparator
            std::set_difference(edge_set_left.begin(), edge_set_left.end(),
                                edge_set_right.begin(), edge_set_right.end(),
                                std::inserter(edge_set_destroy, edge_set_destroy.begin()));
            std::set_difference(edge_set_right.begin(), edge_set_right.end(),
                                edge_set_left.begin(), edge_set_left.end(),
                                std::inserter(edge_set_create, edge_set_create.begin()));

        }

        void print_node_sets(const std::string& name)
        {
            auto print = [](auto&& n, auto& s) {
                std::cout << n << ": ";
                for(auto& i : s) std::cout << i << " ";
                std::cout << "\n";
            };
            std::cout << name << ":\n";
            print("left", node_set_left);
            print("right", node_set_right);
            print("destroy", node_set_destroy);
            print("create", node_set_create);
        }

        void print_edge_sets(const std::string& name)
        {
            auto print_edge_set = [](std::string&& name, auto& edge_set)
            {
                std::cout << name << ": ";
                for(auto& e : edge_set)
                    std::cout << "( " << e.first << ", " << e.second << " ) ";
                std::cout << "\n";
            };

            print_edge_set("left", edge_set_left);
            print_edge_set("right", edge_set_right);
            print_edge_set("create", edge_set_create);
            print_edge_set("destroy", edge_set_destroy);
        }
    };

    /*
     * Grammar analysis currently consists of:
     *
     * 1. Looking at LHS and finding patterns to be fed into search code.
     *    The search is responsible for building a state of the system.
     *
     * 2. Pre-computing rewrite functions
     */
    template <typename GraphType>
    struct AnalyzedGrammar
    {

        // lhs[1].x1
        // std::get<NodeAData>(lhs[1].data).energy
        // 3 -> 1 -> 10 => lhs[m2[m1[1]]]
        std::map<std::string, WithRule<GraphType>> with_rules;
        std::map<std::string, SolvingRule<GraphType>> solving_rules;
        std::vector<std::string> rule_names;
        std::map<std::string, std::map<std::size_t, GraphType>> lhs_connected_components;
        std::map<std::string, Rewrite<GraphType>> with_rewrites;

        //many-to-many relationship of lhs components to minimum set
        std::map<std::string, std::vector<std::size_t>> rule_component;
        std::map<std::size_t, std::vector<std::string>> component_rule;
        std::map<std::size_t, GraphType> unique_components;

        //storing the vertex mapping (connected component unique vertex mappings)
        std::map<std::string, std::vector<std::map<std::size_t, std::size_t>>> ccuv_mappings;

        AnalyzedGrammar() = default;

        AnalyzedGrammar(Grammar<GraphType>& grammar) {
            std::cout << "Creating the grammar data structure used for analysis\n";
            with_rules = grammar.stochastic_rules;
            solving_rules =grammar.deterministic_rules;

            for(auto& item : with_rules)
            {
                rule_names.emplace_back(item.first);
                lhs_connected_components[item.first] = {};
                compute_and_store_connected_components(item.second.lhs_graph, lhs_connected_components[item.first]);
            }
            for(auto& item : solving_rules)
            {
                rule_names.emplace_back(item.first);
                lhs_connected_components[item.first] = {};
                compute_and_store_connected_components(item.second.lhs_graph, lhs_connected_components[item.first]);
            }

            for(auto& [n, c] : lhs_connected_components)
                std::cout << "Rule " << n << " has " << c.size() << " components\n";

            for(auto& [name, rule] : with_rules)
                with_rewrites.insert({name, Rewrite<GraphType>{rule.lhs_graph, rule.rhs_graph}});

            for(auto& item : with_rewrites)
            {
                item.second.print_node_sets(item.first);
            }

            // compute the unique set of connected components
            std::size_t unique_component_key = 0;
            for(auto& [name, components] : lhs_connected_components)
            {

                rule_component.insert({name, {}});
                ccuv_mappings.insert({name, {}});
                for(auto& [ki, ci] : components)
                {
                    bool found = false;
                    for(auto& [kj, cj] : unique_components)
                    {
                        auto matches = YAGL::graph_isomorphism(cj, ci);
                        if(!matches.empty()) found = true;
                        if(found)
                        {
                            std::cout << "Match for " << name << "\n";
                            rule_component[name].push_back(kj);
                            component_rule[kj].push_back(name);
                            ccuv_mappings[name].push_back(matches[0]);
                            break;
                        }
                    }
                    if(!found)
                    {
                        unique_components[unique_component_key] = ci;
                        rule_component[name].push_back(unique_component_key);
                        component_rule[unique_component_key].push_back(name);
                        //not found so we just find a mapping to itself
                        auto matches = YAGL::graph_isomorphism(ci, ci);
                        ccuv_mappings[name].push_back(matches[0]);
                        unique_component_key++;
                    }
                }
            }

            std::cout << "There are " << unique_components.size() << " unique component(s)\n";

            //could be computed during the unique connected component loop
            for(auto& [name, components] : lhs_connected_components)
            {
                for(auto i = 0; i < components.size(); i++)
                {
                    std::cout << "component " << i << " of " << name << " is isomorphic to unique component " << rule_component[name][i] << "\n";
                    for(auto& [k1, k2] : ccuv_mappings[name][i])
                        std::cout << "{ " << k1 << ", " << k2 << " }\n";
                    std::cout << "\n";
                }
            }

        }
    };
}

#endif //DGGML_ANALYZEDGRAMMAR_HPP
