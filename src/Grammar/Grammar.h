#ifndef DGGML_GRAMMAR_H
#define DGGML_GRAMMAR_H

#include <functional>
#include <map>

namespace DGGML {
    template <typename GraphType>
    struct WithRule
    {
        GraphType lhs_graph;
        GraphType rhs_graph;
        std::string rname;

        using GraphMapType = std::map<typename GraphType::key_type, typename GraphType::key_type>;
        std::function<double(GraphType& lhs, GraphMapType& m)> propensity;
        std::function<void(GraphType& lhs, GraphType& rhs, GraphMapType& m)> update;

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

        WithRule<GraphType>& with(std::function<double(GraphType& lhs, GraphMapType& m)>&& p)
        {
            propensity = p;
            return *this;
        }

        WithRule<GraphType>& where(std::function<void(GraphType& lhs, GraphType& rhs, GraphMapType& m)>&& u)
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
        using graph_type = GraphType;
        using key_type = typename GraphType::key_type;

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
