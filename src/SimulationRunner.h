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

            //find single components and compute the anchor
            compute_component_matches();

            //should multicomp matches be copies of single comps or references?
            set_geocell_propensities();

            model->collect();
            model->print_metrics();

            CellList test_cell_list(model->geoplex2D.reaction_grid, model->system_graph, rule_system);

            for(auto& [name, value] : rule_instances)
                std::cout << "So far we have found " << value.size() << " instances of rule " << name << "\n";

            std::function<void(std::string, std::vector<std::size_t>&, int, std::vector<std::size_t>&, decltype(test_cell_list)&, std::size_t)> reaction_instance_backtracker =
                    [&](std::string name, std::vector<std::size_t>& result, int k, std::vector<std::size_t>& pattern, decltype(test_cell_list)& cell_list, std::size_t c)
            {
                if( k == pattern.size())
                {
                    rule_instances[name].push_back(result);
                }
                else {
                    int imin, imax, jmin, jmax;
                    cell_list.getCells(c, imin, imax, jmin, jmax);
                    for (auto i = imin; i < imax; i++) {
                        for (auto j = jmin; j < jmax; j++) {
                            auto nbr_idx = cell_list.cardinalCellIndex(i, j);
                            for (const auto& m2: cell_list.data[nbr_idx]) {
                               bool found = false;
                               for(auto v = 0; v < k; v++)
                               {
                                   if(result[v] == m2.second) found = true;
                               }
                               if(!found && m2.first.type == pattern[k])
                               {
                                   result[k] = m2.second;
                                   reaction_instance_backtracker(name, result, k+1, pattern, cell_list, c);
                               }
                            }
                        }
                    }
                }
            };

            for(auto& m1 : rule_system)
            {
                for(auto& [name, pattern] : compTab.rule_component)
                {
                    int k = 0;
                    std::vector<std::size_t> result;
                    result.resize(pattern.size());

                    if(pattern.size() && m1.first.type == pattern.front())
                    {
                        result[k] = m1.second;
                        k++;
                        auto c = test_cell_list.locate_cell(m1);
                        reaction_instance_backtracker(name, result, k, pattern, test_cell_list, c);
                    }
                }
            }

            auto sum  = 0;
            for(auto& [name, instances] : rule_instances)
            {
                std::cout << "Rule " << name << " has " << instances.size() << " instances\n";
                sum += instances.size();
            }
            std::cout << "There are " << sum << " rule instances in total\n";
//            auto growing_matches = 0;
//            for(auto c = 0; c < test_cell_list.totalNumCells(); c++)
//            {
//                int imin, imax, jmin, jmax;
//                test_cell_list.getCells(c, imin, imax, jmin, jmax);
//                for(auto& [name, patterns] : compTab.components)
//                {
//                    for(auto& p : patterns)
//                    {
//
//                    }
//                }
//
//                /*
//                 * for each pattern p, first search inside the cell
//                 */
//                for(auto i = imin; i < imax; i++)
//                {
//                    for(auto j = jmin; j < jmax; j++)
//                    {
//                        auto nbr_idx = test_cell_list.cardinalCellIndex(i, j);
//                        for(const auto& match1 : test_cell_list.data[c])
//                        {
//                            for(const auto& match2 : test_cell_list.data[nbr_idx])
//                            {
//                                if(match1.first.type == 2 && match2.first.type == 2
//                                   && match1.second != match2.second)
//                                {
//                                    auto a1 = match1.first.anchor;
//                                    auto a2 = match2.first.anchor;
//                                    auto& p1 = model->system_graph.findNode(a1)->second.getData().position;
//                                    auto& p2 = model->system_graph.findNode(a2)->second.getData().position;
//                                    auto d = calculate_distance(p1, p2);
//                                    //std::cout << "Distance: " << d << "\n";
//                                    if(d < test_cell_list.grid._dx)
//                                        growing_matches++;
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//            std::cout << "Potential Growing End Collision Matches Found: " << growing_matches << "\n";

        }

        void run() {

            return;

            using rule_key_t = std::size_t;
            using cplex_key_t = typename DGGML::CartesianComplex2D<>::graph_type::key_type;
            std::map<cplex_key_t, std::vector<rule_key_t>> rule_map;

            for(auto& [key, value] : model->geoplex2D.graph.getNodeSetRef())
                rule_map.insert({key, {}});

            auto count2d = std::count_if(rule_map.begin(), rule_map.end(), [](auto& iter){ return iter.first == 0; });
            std::cout << "Count2D: " << count2d << "\n";
            for(auto i = 0; i <= model->settings.NUM_STEPS; i++)
            {
                std::cout << "Running step " << i << "\n";
            }
        }
    private:

        /*
         * Grammar analysis currently consists of:
         *
         * 1. Looking at LHS and finding patterns to be fed into search code.
         *    The search is responsible for building a state of the system.
         *
         * 2. Pre-computing rewrite functions
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

        void compute_component_matches()
        {
            //single component matches are space invariant, but multi-component matches
            //are not
            std::vector<std::vector<key_type>> match_set;
            //TODO: I think we need to store the ordering or the rooted spanning tree
            for(auto& [k, pattern] : compTab.components)
            {
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
                    Instance<typename ModelType::key_type> inst;
                    inst.match = match;
                    inst.type = k;
                    inst.anchor = match[0];
                    rule_system.push_back(inst);
                }
                std::cout << "Found " << matches.size() << " instances\n";
            }

            //for(auto& [k, p] : instances)
            for(auto& [k, pattern] : compTab.components)
            {
                std::cout << "Component " << k << " has " << rule_system.count(k) << " instances\n";
            }

//            for(auto& [name, rule] : model->gamma.stochastic_rules)
//            {
//                rule_instances[name] = {};
//                if(compTab.rule_component[name].size() == 1)
//                {
//                    std::cout << name << " is a " << compTab.rule_component[name].size() << " component rule\n";
//                    for(auto& ckey : compTab.rule_component[name])
//                    {
//                        std::cout << "There are " << rule_system.count(ckey) << " instances of component " << ckey << "\n";
//                        for(auto& cinst : rule_system)
//                        {
//                            if(cinst.first.type == ckey)
//                                rule_instances[name].push_back({cinst.second});
//                        }
//                    }
//                }
//
//            }
        }
        std::shared_ptr<ModelType> model;

        //actually a component table and a "join"
        struct RuleComponentTable{
                std::map<std::string, std::vector<std::size_t>> rule_component;
                std::map<std::size_t, std::vector<std::string>> component_rule;
                std::map<std::size_t, graph_type> components;
        } compTab;

       std::map<std::string, std::vector<std::vector<std::size_t>>> rule_instances;
       std::map<std::size_t, std::vector<typename ModelType::key_type>> ordering;

        //Grammar gamma;
        RuleSystem<typename ModelType::key_type> rule_system;
        std::map<gplex_key_type, std::pair<double, double>> geocell_progress;
        DGGML::VtkFileWriter<typename ModelType::graph_type> vtk_writer;
        std::string results_dir_name;
    };
}

#endif //DGGML_SIMULATIONRUNNER_H
