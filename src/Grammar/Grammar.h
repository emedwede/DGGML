#ifndef DGGML_GRAMMAR_H
#define DGGML_GRAMMAR_H

namespace DGGML {

    using key_type = DGGML::Plant::mt_key_type;
    using data_type = DGGML::Plant::MT_NodeData;

    using graph_type = YAGL::Graph<key_type, data_type>;
    using node_type = YAGL::Node<key_type, data_type>;

    struct LHS
    {
        graph_type graph;
        std::vector<graph_type> components;

        LHS() {}

        LHS(graph_type& graph) : graph(graph)
        {
            //build the list of components
            std::unordered_set<key_type> visited;
            std::size_t count = 0; //no connected components found to start
            for(auto i = graph.node_list_begin(); i != graph.node_list_end(); i++)
            {
                auto v = i->first;
                //node hasn't been visited so it must be the start of a new connected component
                if(visited.find(v) == visited.end())
                {
                    std::vector<key_type> path;
                    //we could use whatever search method we feel like
                    YAGL::impl_iterative_bfs(graph, v, visited, path);
                    components.push_back(YAGL::induced_subgraph(graph, path));
                    count++;
                }
            }
        }
    };

    //should LHS and RHS be the same object?
    struct RHS
    {
        graph_type graph;
        std::vector<graph_type> components;

        RHS() {}

        RHS(graph_type& graph) : graph(graph) {}
    };

    struct RuleType
    {
        LHS lhs;
        RHS rhs;
        std::string name;

        RuleType() {}

        RuleType(std::string name, graph_type& lhs_graph, graph_type& rhs_graph) : name(name), lhs(lhs_graph), rhs(rhs_graph) {}
    };

    struct Grammar
    {
        std::map<std::string, RuleType> rule_set;
        std::map<std::string, std::vector<std::size_t>> rule_component;
        std::map<std::size_t, std::vector<std::string>> component_rule;
        std::vector<graph_type> minimal_set;

        void addRule(RuleType& r)
        {
            if(rule_set.find(r.name) == rule_set.end())
            {
                rule_set[r.name] = r;
                incremental_build_min_set(r);
            }
            else
                std::cout << "name already exists in rule_set\n";
        }

        void incremental_build_min_set(RuleType& r)
        {

            for(auto& c : r.lhs.components)
            {
                auto idx = 0;
                bool found = false;
                for(auto i = 0; i < minimal_set.size(); i++)
                {
                    auto matches = YAGL::subgraph_isomorphism2(c, minimal_set[i]);
                    if(!matches.empty())
                    {
                        idx = i;
                        found = true;
                        break;
                    }
                }
                if(!found)
                {
                    //add forward and backward references
                    rule_component[r.name].push_back(minimal_set.size());
                    component_rule[minimal_set.size()].push_back(r.name);
                    minimal_set.push_back(c);
                }
                else //found
                {
                    rule_component[r.name].push_back(idx);
                    component_rule[idx].push_back(r.name);
                }
            }
        }

        void print()
        {
            std::cout << "Rules: ";
            for(auto& rule : rule_set)
                std::cout << "\t" << rule.first << "\n";
        }
    };

    void print_mapping(Grammar& gamma)
    {
        std::cout << "\nRule Component Mapppings:\n";
        for(auto& [key, value] : gamma.rule_set)
        {
            std::cout << "\t" << key << ": ";
            for(auto& item : gamma.rule_component[key])
            {
                std::cout << item << " ";
            } std::cout << "\n";
        }

        std::cout << "\nComponent Rule Mapppings:\n";
        for(auto key = 0; key < gamma.minimal_set.size(); key++)
        {
            std::cout << "\t" << key << ": ";
            for(auto& item : gamma.component_rule[key])
            {
                std::cout << item << " ";
            } std::cout << "\n";
        }
    }

} // DGGML

#endif //DGGML_GRAMMAR_H
