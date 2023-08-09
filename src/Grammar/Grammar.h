#ifndef DGGML_GRAMMAR_H
#define DGGML_GRAMMAR_H

#include <functional>

namespace DGGML {

//    struct LHS
//    {
//        graph_type graph;
//        std::vector<graph_type> components;
//
//        LHS() {}
//
//        LHS(graph_type& graph) : graph(graph)
//        {
//            //build the list of components
//            std::unordered_set<key_type> visited;
//            std::size_t count = 0; //no connected components found to start
//            for(auto i = graph.node_list_begin(); i != graph.node_list_end(); i++)
//            {
//                auto v = i->first;
//                //node hasn't been visited so it must be the start of a new connected component
//                if(visited.find(v) == visited.end())
//                {
//                    std::vector<key_type> path;
//                    //we could use whatever search method we feel like
//                    YAGL::impl_iterative_bfs(graph, v, visited, path);
//                    components.push_back(YAGL::induced_subgraph(graph, path));
//                    count++;
//                }
//            }
//        }
//    };
//
//    //should LHS and RHS be the same object?
//    struct RHS
//    {
//        graph_type graph;
//        std::vector<graph_type> components;
//
//        RHS() {}
//
//        RHS(graph_type& graph) : graph(graph) {}
//    };
//
//    struct RuleType
//    {
//        LHS lhs;
//        RHS rhs;
//        std::string name;
//
//        RuleType() {}
//
//        RuleType(std::string name, graph_type& lhs_graph, graph_type& rhs_graph) : name(name), lhs(lhs_graph), rhs(rhs_graph) {}
//    };

    template <typename GraphType>
    struct WithRule
    {
        GraphType lhs_graph;
        GraphType rhs_graph;
        std::string rname;
        std::function<double(GraphType& lhs)> propensity;
        std::function<void(GraphType& lhs, GraphType& rhs)> update;

        WithRule() {}

        WithRule(std::string rname, GraphType& lhs_graph, GraphType& rhs_graph) :
        rname(rname), lhs_graph(lhs_graph), rhs_graph(rhs_graph) {}

        // Fluent interface
        WithRule<GraphType>& name(std::string&& n)
        {
            rname = n;
            return *this;
        }
        WithRule<GraphType>& lhs(GraphType& graph)
        {
            lhs_graph = graph;
            return *this;
        }

        WithRule<GraphType>& rhs(GraphType& graph)
        {
            rhs_graph = graph;
            return *this;
        }

        WithRule<GraphType>& with(std::function<double(GraphType& lhs)>&& p)
        {
            propensity = p;
            return *this;
        }

        WithRule<GraphType>& where(std::function<void(GraphType& lhs, GraphType& rhs)>&& u)
        {
            update = u;
            return *this;
        }

    };

    template <typename GraphType>
    struct SolvingRule
    {
        GraphType lhs;
        GraphType rhs;
        std::string name;
        std::function<double(GraphType& lhs)> ode;

        SolvingRule() {}
        SolvingRule(std::string name, GraphType& lhs_graph, GraphType& rhs_graph) : name(name), lhs(lhs_graph), rhs(rhs_graph) {}
    };

    //TODO: we need to make sure there is a numbered graph class for rule definitions
    // the numbered graphs has a numbering that is the maximum size of the left and right
    // hand side. All integer labels must differ from each other and they should cover the
    // numbering space and have no gaps. We also need some helper functions. One that
    // returns input only’s = destroyed, shared = maybe changed params not created or
    // destroyed, and the output helper for ones that are only created
    template<typename GraphType>
    struct Grammar
    {
        std::map<std::string, WithRule<GraphType>> stochastic_rules;
        std::map<std::string, SolvingRule<GraphType>> deterministic_rules;

        void addRule(WithRule<GraphType>& r)
        {
            if(stochastic_rules.find(r.rname) == stochastic_rules.end())
            {
                stochastic_rules[r.rname] = r;
            }
            else
                std::cout << "stochastic rule name already exists in rule_set\n";
        }

        void addRule(SolvingRule<GraphType>& r)
        {
            if(deterministic_rules.find(r.name) == deterministic_rules.end())
            {
                deterministic_rules[r.name] = r;
            }
            else
                std::cout << "deterministic rule name already exists in rule_set\n";
        }

        void print()
        {
            std::cout << "Stochastic Rules: \n";
            for(auto& rule : stochastic_rules)
                std::cout << "\t" << rule.first << "\n";

            std::cout << "Deterministic Rules: \n";
            for(auto& rule : deterministic_rules)
                std::cout << "\t" << rule.first << "\n";
        }
    };

} // DGGML

#endif //DGGML_GRAMMAR_H
