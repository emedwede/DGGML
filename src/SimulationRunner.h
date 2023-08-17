#ifndef DGGML_SIMULATIONRUNNER_H
#define DGGML_SIMULATIONRUNNER_H

#include <iostream>
#include <memory>
#include <map>
#include <random>
#include <chrono>
#include <string>
#include <filesystem>

#include "DggModel.hpp"
#include "YAGL_Graph.hpp"
#include "YAGL_Node.hpp"
#include "YAGL_Algorithms.hpp"
#include "Utlities/VtkWriter.hpp"
#include "CartesianComplex2D.hpp"
#include "RuleSystem.hpp"
#include "Utlities/MathUtils.hpp"


namespace DGGML {
    template<typename ModelType>
    class SimulationRunner {
    public:
        using key_type = typename ModelType::key_type;
        using gplex_key_type = typename DGGML::ExpandedComplex2D<>::graph_type::key_type;
        //using data_type = typename ModelType::graph_type;
        //using graph_type = YAGL::Graph<key_type, data_type>;
        using graph_type = typename ModelType::graph_type;
        using node_type = typename graph_type::node_type;
        explicit SimulationRunner(const ModelType& model) : model(std::make_shared<ModelType>(model)) {}

        //May want to make this part of the constructor instead
        void initialize()
        {
            //Should the model initialize even happen here? or sooner
            std::cout << "Initializing " << model->name << "\n";
            model->initialize();
            std::cout << "Printing grammar rules...\n";
            model->gamma.print();

            //order matters here, which indicates maybe I should have a
            //file writer class which initializes with the save directory?
            create_save_directory();
            write_cell_complex();
            write_system_graph(0);

            analyze_grammar();
            compute_matches();

            set_geocell_propensities();

            model->collect();
            model->print_metrics();

        }

        void run() {

            for (auto &[n, r]: compTab.rule_component)
            {
                if(r.size() == 1)
                {
                    std::cout << "Running test for rule " << n << " \n";
                    for(auto& c : r)
                    {
                        auto result = 0.0;
                        for(auto& m : instances[c])
                        {
                            auto lhs_match = YAGL::induced_subgraph(model->system_graph, m);
                            typename WithRule<graph_type>::GraphMapType temp; //need actual ordering
                            result += model->gamma.stochastic_rules[n].propensity(lhs_match, temp);
                        }
                        std::cout << "Result: " << result << "\n";
                    }
                }
            }

            return;
        }
    private:

        /*
         * Grammar analysis currently consists of:
         *
         * 1. Looking at LHS and finding patterns to be fed into search code.
         *    The search is responsible for building a state of the system.
         *
         * 2. C Pre computing rewrite functions
         */
        void analyze_grammar()
        {
            std::cout << "Performing grammar analysis\n";
            Grammar<graph_type>& gamma = model->gamma;

            auto& components = compTab.components;
            auto& rule_component = compTab.rule_component;
            auto& component_rule = compTab.component_rule;
            std::map<std::size_t, graph_type> minimal_set;

            std::size_t ckey = 0;
            for(auto& [name, rule] : gamma.stochastic_rules)
            {
                graph_type& graph = rule.lhs_graph;
                rule_component.insert({name, {}});
                //build the list of components
                std::unordered_set<key_type> visited;
                std::size_t count = 0; //no connected components found to start
                for(auto i = graph.node_list_begin(); i != graph.node_list_end(); i++) {
                    auto v = i->first;
                    //node hasn't been visited, so it must be the start of a new connected component
                    if (visited.find(v) == visited.end()) {
                        std::vector<key_type> path;
                        //we could use whatever search method we feel like
                        YAGL::impl_iterative_bfs(graph, v, visited, path);
                        auto c = YAGL::induced_subgraph(graph, path);

                        auto is_isomorphic = [&](auto& a) {
                            auto matches = YAGL::graph_isomorphism(a.second, c);
                            return !matches.empty();
                        };

                        if(!components.empty())
                        {
                            auto it = std::find_if(components.begin(), components.end(), is_isomorphic);
                            if(it == components.end()) {
                                components[ckey] = c;
                                rule_component[name].push_back(ckey);
                                component_rule[ckey].push_back(name);
                                ckey++;
                            } else
                            {
                                rule_component[name].push_back(it->first);
                                component_rule[it->first].push_back(name);
                            }
                        }
                        else
                        {
                            components[ckey] = c;
                            rule_component[name].push_back(ckey);
                            component_rule[ckey].push_back(name);
                            ckey++;
                        }
                        count++;
                    }
                }
                std::cout << name << " has " << count << " components { ";
                for(auto& item : rule_component[name]) std::cout << item << " ";
                std::cout << "}\n";
            }
            std::cout << components.size() << " unique components found\n";

            for(auto& [k, v] : component_rule)
            {
                    std::cout << "Component " << k << ": { ";
                    for(auto& s : v)
                        std::cout << s << " ";
                    std::cout << "}\n";
            }

            //compute the graph rewrite operations, the numbering of the graphs matters here
            for(auto& [name, rule] : gamma.stochastic_rules)
            {
                graph_type& lhs_graph = rule.lhs_graph;
                graph_type& rhs_graph = rule.rhs_graph;

                std::set<key_type> left_keys, right_keys, create, destroy;
                for(auto& node : lhs_graph.getNodeSetRef())
                    left_keys.insert(node.first);
                for(auto& node : rhs_graph.getNodeSetRef())
                    right_keys.insert(node.first);

                std::set_difference(left_keys.begin(), left_keys.end(),
                                    right_keys.begin(), right_keys.end(),
                                    std::inserter(destroy, destroy.begin()));
                std::set_difference(right_keys.begin(), right_keys.end(),
                                    left_keys.begin(), left_keys.end(),
                                    std::inserter(create, create.begin()));

                auto print = [](auto&& n, auto& s) {
                    std::cout << n << ": ";
                    for(auto& i : s) std::cout << i << " ";
                    std::cout << "\n";
                };

                //print("left", left_keys);
                //print("right", right_keys);
                //print("destroy", destroy);
                //print("create", create);

                auto lhs_graph_copy = lhs_graph;
                for(auto& k : create) {
                    auto n = rhs_graph.findNode(k)->second;
                    std::cout << n.getData().type << "\n";
                    lhs_graph_copy.addNode(n);
                }

                //TODO: finish computing rewrites for edge set

            }

        }

        void create_save_directory()
        {
            std::cout << "Cleaning up old results folder if it exists and creating a new one\n";
            results_dir_name = model->name + "_results";
            std::filesystem::remove_all(results_dir_name);
            std::filesystem::create_directory(results_dir_name);
        }

        void set_geocell_propensities()
        {
            std::cout << "Setting intial cell propensities to zero\n";
            for(auto& [key, value] : model->geoplex2D.graph.getNodeSetRef())
                geocell_progress[key] = {0.0, 0.0};
        }

        void write_cell_complex()
        {
            //Save expanded cell complex graph
            DGGML::VtkFileWriter<typename DGGML::ExpandedComplex2D<>::types::graph_type> writer;
            writer.save(model->geoplex2D.getGraph(), results_dir_name+"/factory_geoplex");
            DGGML::GridFileWriter grid_writer;
            grid_writer.save({model->geoplex2D.reaction_grid,
                              model->geoplex2D.dim_label},
                             results_dir_name+"/expanded_cell_complex");
        }

        void write_system_graph(std::size_t step)
        {
            std::string title = results_dir_name+"/simulation_step_";
            std::cout << "Saving the initial state of the system graph\n";
            vtk_writer.save(model->system_graph, title+std::to_string(step));
        }

        void compute_matches()
        {
            //single component matches are space invariant, but multi-component matches
            //are not
            std::vector<std::vector<key_type>> match_set;
            //TODO: I think we need to store the ordering or the rooted spanning tree
            for(auto& [k, pattern] : compTab.components)
            {
                instances[k] = {};
                //need to actually get the ordering for the mapping
                auto matches = YAGL::subgraph_isomorphism2(pattern, model->system_graph);
                if(!matches.empty())
                for(auto& item : matches)
                {
                    std::vector<typename ModelType::key_type> match;
                    for(auto& [key, value] : item)
                    {
                        match.push_back(value);
                    }
                    instances[k].push_back(std::move(match));
                }
                std::cout << "Found " << matches.size() << " instances\n";
            }

            for(auto& [k, p] : instances)
            {
                std::cout << "Component " << k << " has " << p.size() << " instances\n";
            }
        }
        std::shared_ptr<ModelType> model;

        //actually a component table and a "join"
        struct RuleComponentTable{
                std::map<std::string, std::vector<std::size_t>> rule_component;
                std::map<std::size_t, std::vector<std::string>> component_rule;
                std::map<std::size_t, graph_type> components;
        } compTab;

       std::map<std::size_t, std::vector<std::vector<typename ModelType::key_type>>> instances;
       std::map<std::size_t, std::vector<typename ModelType::key_type>> ordering;

        //Grammar gamma;
        RuleSystem<Plant::mt_key_type> rule_system;
        std::map<gplex_key_type, std::pair<double, double>> geocell_progress;
        DGGML::VtkFileWriter<typename ModelType::graph_type> vtk_writer;
        std::string results_dir_name;
    };
}

#endif //DGGML_SIMULATIONRUNNER_H
