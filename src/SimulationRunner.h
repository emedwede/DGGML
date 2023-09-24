#ifndef DGGML_SIMULATIONRUNNER_H
#define DGGML_SIMULATIONRUNNER_H

#include <iostream>
#include <memory>
#include <map>
#include <set>
#include <random>
#include <chrono>
#include <string>
#include <filesystem>

#include "DggModel.hpp"
#include "YAGL_Graph.hpp"
#include "YAGL_Node.hpp"
#include "YAGL_Algorithms.hpp"
#include "Utlities/VtkWriter.hpp"
#include "ExpandedComplex2D.hpp"
#include "CartesianComplex2D.hpp"
#include "CellList.hpp"
#include "RuleSystem.hpp"
#include "Grammar.h"
#include "Utlities/MathUtils.hpp"
#include "CartesianHashFunctions.hpp"
#include "ApproximateSSA.hpp"
#include "AnalyzedGrammar.hpp"


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
        using rule_key_t = std::size_t;
        using cplex_key_t = typename DGGML::CartesianComplex2D<>::graph_type::key_type;

        explicit SimulationRunner(const ModelType& _model) :
        model(std::make_shared<ModelType>(_model)) {}

        //TODO: we probably want to make this part of the constructor instead or
        // even better yet, move most of these objects to be external to the class and injected upon
        // construction
        void initialize()
        {
            //Should the model initialize even happen here? or sooner
            std::cout << "Initializing " << model->name << "\n";
            model->initialize();

            std::cout << "Printing grammar rules...\n";
            model->gamma.print();

            std::cout << "Performing grammar analysis and building analyzed grammar data structure\n";
            grammar_analysis = AnalyzedGrammar<graph_type>(model->gamma);
            std::cout << "Grammar analysis complete\n";

            //building rule map sets
            for(auto& [key, value] : model->geoplex2D.graph.getNodeSetRef())
                rule_map.insert({key, {}});

            //order matters here, which indicates maybe I should have a
            //file writer class which initializes with the save directory?
            create_save_directory();
            write_cell_complex();
            write_system_graph(0);

            compute_single_component_matches();

            set_geocell_propensities();

            model->collect();
            model->print_metrics();

            // A cell list is used to accelerate the spatial geometric search, but
            // we could use other methods like a bounding volume hierarchy(ArborX), kd-trees etc
            // see books on collision detection like gpu-gems or Real-Time Collision Detection
            CellList test_cell_list(model->geoplex2D.reaction_grid, model->system_graph, rule_system);

            compute_all_rule_instances(test_cell_list);

            map_rule_instances_to_geocells();

            //TODO: deprecate and refactor
            //builds a list of interior geocells to iterate
            build_bucketsND(bucketsND, model->geoplex2D);
        }

        void run() {
            //return;

            for(auto i = 0; i <= 2; i++) {
                auto countNd = std::count_if(rule_map.begin(), rule_map.end(),
                                             [&](auto &iter) {
                    return model->geoplex2D.getGraph().findNode(iter.first)->second.getData().type == i;
                });
                std::cout << "Count" << (2 - i) << "D: " << countNd << "\n";
            }

            for(auto i = 0; i <= model->settings.NUM_STEPS; i++)
            {
                std::cout << "Running step " << i << "\n";

                double tot_time = 0.0;

                for(int d = 0; d < 3; d++)
                {
                    double dim_time = 0.0;

                    if(d == 1 || d == 2) continue;
                    std::cout << "Running the Hybrid ODES/SSA inner loop " << (2 - d) << "D phase\n";
                    for(auto& bucket : bucketsND[d])
                    {
                        auto k = bucket.first;
                        auto start = std::chrono::high_resolution_clock::now();
                        approximate_ssa(rule_system, grammar_analysis, rule_map, rule_instances, model, k, geocell_progress[k]);
                        return;
                        auto stop = std::chrono::high_resolution_clock::now();
                        auto duration =
                                std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
                        std::cout << "Cell " << k << " took " << duration.count()
                                  << " milliseconds and has a current tau "
                                  << geocell_progress[k].first << "\n";
                        dim_time += duration.count();
                    }

                    tot_time += dim_time;
                    std::cout << (2 - d) << "D took " << dim_time << " milliseconds\n";

                    std::cout << "Synchronizing work\n";
                }

                std::cout << "Running the checkpointer\n";
                write_system_graph(0);

                std::cout << "Total dimensional time is " << tot_time << " milliseconds\n";
                //time_count.push_back(tot_time);
                return;
            }
        }
    private:

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

        void compute_single_component_matches()
        {
            //single component matches are space invariant, but multi-component matches
            //are not
            std::vector<std::vector<key_type>> match_set;
            //TODO: I think we need to store the ordering or the rooted spanning tree
            for(auto& [k, pattern] : grammar_analysis.unique_components)
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
            for(auto& [k, pattern] : grammar_analysis.unique_components)
            {
                std::cout << "Component " << k << " has " << rule_system.count(k) << " instances\n";
            }
        }


        template<typename CellListType> //TODO: fix me the hacky template fix for a type we should know!
        void compute_all_rule_instances(CellListType& test_cell_list)
        {
            for(auto& [name, rule] : grammar_analysis.with_rules) {
                auto count = std::count_if(rule_instances.begin(), rule_instances.end(), [&](auto& iter) { return iter.second.name == name; });
                std::cout << "So far we have found " << count << " instances of with rule " << name << "\n";
            }
            for(auto& [name, rule] : grammar_analysis.solving_rules) {
                auto count = std::count_if(rule_instances.begin(), rule_instances.end(), [&](auto& iter) { return iter.second.name == name; });
                std::cout << "So far we have found " << count << " instances of solving rule " << name << "\n";
            }
            // This is a recursive backtracking function, and it is memory efficient because it does a DFS (inorder traversal)
            // so only one vector is need for finding the result
            // TODO: Check for bugs patterns > 2, currently some permutations may not be picked up correctly
            std::function<void(std::string, std::vector<std::size_t>&, int, std::vector<std::size_t>&, decltype(test_cell_list)&, std::size_t)> reaction_instance_backtracker =
                    [&](std::string name, std::vector<std::size_t>& result, int k, std::vector<std::size_t>& pattern, decltype(test_cell_list)& cell_list, std::size_t c)
                    {
                        if( k == pattern.size())
                        {
                            auto generated_key = instance_key_gen.get_key();
                            rule_instances[generated_key].name = name;
                            rule_instances[generated_key].category = "stochastic";
                            rule_instances[generated_key].components = result;
                            rule_instances[generated_key].anchor = rule_system[result[0]].anchor;
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

                                            auto a1 = rule_system[result[0]].anchor;
                                            auto a2 = m2.first.anchor;
                                            auto& p1 = model->system_graph.findNode(a1)->second.getData().position;
                                            auto& p2 = model->system_graph.findNode(a2)->second.getData().position;
                                            auto d = calculate_distance(p1, p2);

                                            // TODO: a constraint for nearness is needed, there are other ways to do it
                                            //  this is just a placeholder
                                            if(d < model->settings.MAXIMAL_REACTION_RADIUS)
                                            {
                                                result[k] = m2.second;
                                                reaction_instance_backtracker(name, result, k + 1, pattern, cell_list, c);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    };

            for(auto& m1 : rule_system)
            {
                for(auto& [name, pattern] : grammar_analysis.rule_component)
                {
                    //skip solving rules for now
                    if(grammar_analysis.with_rules.find(name) == grammar_analysis.with_rules.end())
                        continue;
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
            for(auto& name : grammar_analysis.rule_names) {
                auto count = std::count_if(rule_instances.begin(), rule_instances.end(), [&](auto& iter) { return iter.second.name == name; });
                std::cout << "Rule " << name << " has " << count << " instances\n";
                sum += count;
            }
            std::cout << "There are " << sum << " rule instances in total\n";
        }

        void map_rule_instances_to_geocells()
        {
            //depending on the complexity, map anchor nodes to reaction subcells or
            // alternative map reaction instances directly using a computed position
            // or anchor node. We could also use a participation list etc, but for now
            // we've reduced the problem to sorting a whole reaction by an anchor node
            for(auto& [key, inst] : rule_instances)
            {
                auto& node_data = model->system_graph.findNode(inst.anchor)->second.getData();

                double xp = node_data.position[0];
                double yp = node_data.position[1];
                int ic, jc;

                model->geoplex2D.reaction_grid.locatePoint(xp, yp, ic, jc);
                auto cardinal = model->geoplex2D.reaction_grid.cardinalCellIndex(ic, jc);
                // here would be the anchor list reduction, plus whatever other code needs to be added
                //anchor_list.insert({match.first.anchor, cardinal});
                auto max_cell = model->geoplex2D.cell_label[cardinal];
                rule_map[max_cell].push_back(key);
            }

            auto sum = 0;
            for(auto& [key, value] : rule_map)
            {
                sum += value.size();
            }
            std::cout << "Total rules mapped " << sum << "\n";
        }

        //TODO: maybe not shared, but unique?
        std::shared_ptr<ModelType> model;

        AnalyzedGrammar<graph_type> grammar_analysis;

        KeyGenerator<std::size_t> instance_key_gen;
        struct RuleInstType {
            std::string category;
            std::string name;
            std::vector<std::size_t> components;
            key_type anchor;
        };
        std::map<std::size_t, RuleInstType> rule_instances;
        std::map<cplex_key_t, std::vector<rule_key_t>> rule_map;
        std::map<std::size_t, std::vector<typename ModelType::key_type>> ordering;

        std::map<gplex_key_type, std::vector<key_type>> bucketsND[3];
        //Grammar gamma;
        RuleSystem<typename ModelType::key_type> rule_system;
        std::map<gplex_key_type, std::pair<double, double>> geocell_progress;
        DGGML::VtkFileWriter<typename ModelType::graph_type> vtk_writer;
        std::string results_dir_name;
    };
}

#endif //DGGML_SIMULATIONRUNNER_H
