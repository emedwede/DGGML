#ifndef DGGML_GRAMMAR_H
#define DGGML_GRAMMAR_H

#include <functional>
#include <map>
#include <set>
#include <nvector/nvector_serial.h>

namespace DGGML {

    template <typename GraphType>
    struct RuleBase
    {
        using key_type = typename GraphType::key_type;
        std::string rname;

        GraphType lhs_graph;
        GraphType rhs_graph;

        //TODO: figure out why I need the base constructor and can't delete
        //RuleBase() = delete;
        //RuleBase() = default;

        explicit RuleBase(std::string rname, GraphType& lhs_graph, GraphType& rhs_graph) :
        rname(rname), lhs_graph(lhs_graph), rhs_graph(rhs_graph) {}
    };

    template <typename GraphType>
    struct WithRule : RuleBase<GraphType>
    {
        using GraphMapType = std::map<typename GraphType::key_type, typename GraphType::key_type>;
        using propensity_t = std::function<double(GraphType& lhs, GraphMapType& m)>;
        using update_t =  std::function<void(GraphType& lhs, GraphType& rhs, GraphMapType& m1, GraphMapType& m2)>;

        propensity_t propensity;
        update_t update;

        //TODO: figure out why I need the base constructor and can't delete
        //WithRule() = delete;
        //WithRule() {};

        explicit WithRule(std::string rname, GraphType& lhs_graph, GraphType& rhs_graph, propensity_t&& p, update_t&& u)
            : RuleBase<GraphType>(rname, lhs_graph, rhs_graph), propensity(p), update(u) {}

    };

    template <typename GraphType>
    struct SolvingRule : RuleBase<GraphType>
    {
        using GraphMapType = std::map<typename GraphType::key_type, typename GraphType::key_type>;
        using initial_condition_t = std::function<void(GraphType& lhs, GraphMapType& m1, N_Vector ydot, std::size_t offset)>;
        using solving_t = std::function<void(GraphType& lhs, N_Vector ydot)>;

        initial_condition_t ic;
        solving_t ode;
        std::size_t num_eq;

        explicit SolvingRule(std::string rname, GraphType& lhs_graph, GraphType& rhs_graph,
                             std::size_t num_eq, initial_condition_t&& ic, solving_t&& ode)
                : RuleBase<GraphType>(rname, lhs_graph, rhs_graph), num_eq(num_eq), ic(ic), ode(ode) {}
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
                //using insert, vs [] = because the [] requires they map_type to be default constructible
                stochastic_rules.insert({r.rname, r});
            }
            else
                std::cout << "stochastic rule name already exists in rule_set\n";
        }

        void addRule(SolvingRule<GraphType>& r)
        {
            if(deterministic_rules.find(r.rname) == deterministic_rules.end())
            {
                //using insert, vs [] = because the [] requires they map_type to be default constructible
                deterministic_rules.insert({r.rname, r});
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
