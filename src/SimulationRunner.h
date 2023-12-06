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
#include "ComponentMatchMap.hpp"
#include "RuleMatchMap.hpp"
#include "Grammar.h"
#include "Utlities/MathUtils.hpp"
#include "CartesianHashFunctions.hpp"
#include "ApproximateSSA.hpp"
#include "AnalyzedGrammar.hpp"
#include "pattern_matching.hpp"
#include "phi_functions.hpp"

namespace DGGML {
    template<typename ModelType>
    class SimulationRunner {
    public:
        using model_type = ModelType;
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
            cell_list = CellList<graph_type>(model->geoplex2D.reaction_grid, &model->system_graph, &component_matches);
            compute_all_rule_matches();

            map_rule_matches_to_geocells();

            //TODO: deprecate and refactor
            //builds a list of interior geocells to iterate
            build_bucketsND(bucketsND, model->geoplex2D);
        }

        void run() {
            //return;

            for(auto i = 0; i <= 2; i++) {
                auto countNd = std::count_if(rule_map.begin(), rule_map.end(),
                                             [&](auto &iter) {
                    return model->geoplex2D.getGraph().findNode(iter.first)->second.getData().type == i && model->geoplex2D.getGraph().findNode(iter.first)->second.getData().interior;
                });
                std::cout << "Total number of expanded " << (2 - i) << "D cells: " << countNd << "\n";
            }

            for(auto i = 0; i <= model->settings.NUM_STEPS; i++)
            {
                std::cout << "Running step " << i << "\n";

                double tot_time = 0.0;

                for(int d = 0; d < 3; d++)
                {
                    double dim_time = 0.0;

                    //TODO: for 1D and 0D, we may not have to re-bin the components into cells due to diffusion
                    // if we make the depth of reaction cells we search for 1D and 0D deeper i.e. search 2 away vs. 1
                    if(d == 1 || d == 2) continue;
                    std::cout << "Running the Hybrid ODES/SSA inner loop " << (2 - d) << "D phase\n";
                    for(auto& bucket : bucketsND[d])
                    {
                        auto k = bucket.first;
                        auto start = std::chrono::high_resolution_clock::now();
                        approximate_ssa(component_matches, grammar_analysis,
                                        rule_map, rule_matches,
                                        model, k, geocell_progress[k], cell_list);
                        auto stop = std::chrono::high_resolution_clock::now();
                        auto duration =
                                std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
                        std::cout << "Cell " << k << " took " << duration.count()
                                  << " milliseconds and has a current tau "
                                  << geocell_progress[k].first << "\n";
                        dim_time += duration.count();
                        break;//return;
                    }

                    tot_time += dim_time;
                    std::cout << (2 - d) << "D took " << dim_time << " milliseconds\n";

                    std::cout << "Synchronizing work\n";
                }

                std::cout << "Running the checkpointer\n";
                write_system_graph(i+1);

                std::cout << "Total dimensional time is " << tot_time << " milliseconds\n";
                //time_count.push_back(tot_time);
                if (i == 10) return;
            }
            //return;
        }
    //private:

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
                    ComponentMatch<typename ModelType::key_type> inst;
                    inst.match = match;
                    inst.type = k;
                    inst.anchor = match[0];
                    component_matches.insert(inst);
                }
                std::cout << "Found " << matches.size() << " instances\n";
            }

            //for(auto& [k, p] : instances)
            for(auto& [k, pattern] : grammar_analysis.unique_components)
            {
                std::cout << "Component " << k << " has " << component_matches.count(k) << " instances\n";
            }
        }

        void compute_all_rule_matches()
        {
            for(auto& [name, rule] : grammar_analysis.with_rules) {
                auto count = std::count_if(rule_matches.begin(), rule_matches.end(), [&](auto& iter) { return iter.second.name == name; });
                std::cout << "So far we have found " << count << " instances of with rule " << name << "\n";
            }
            for(auto& [name, rule] : grammar_analysis.solving_rules) {
                auto count = std::count_if(rule_matches.begin(), rule_matches.end(), [&](auto& iter) { return iter.second.name == name; });
                std::cout << "So far we have found " << count << " instances of solving rule " << name << "\n";
            }
            // This is a recursive backtracking function, and it is memory efficient because it does a DFS (inorder traversal)
            // so only one vector is need for finding the result


            for(auto m1 = component_matches.begin(); m1 != component_matches.end(); m1++)
            {
                for(auto& [name, pattern] : grammar_analysis.rule_component)
                {
                    int k = 0;
                    std::vector<std::size_t> result;
                    result.resize(pattern.size());

                    if(pattern.size() && m1->second.type == pattern.front())
                    {
                        result[k] = m1->first;
                        k++;
                        auto c = cell_list.locate_cell(m1);
                        reaction_instance_backtracker(*this, name, result, k, pattern, c);
                    }
                }
            }

            auto sum  = 0;
            for(auto& name : grammar_analysis.rule_names) {
                auto count = std::count_if(rule_matches.begin(), rule_matches.end(), [&](auto& iter) { return iter.second.name == name; });
                std::cout << "Rule " << name << " has " << count << " instances\n";
                sum += count;
            }
            std::cout << "There are " << sum << " rule instances in total\n";
        }

        void map_rule_matches_to_geocells()
        {
            //depending on the complexity, map anchor nodes to reaction subcells or
            // alternative map reaction instances directly using a computed position
            // or anchor node. We could also use a participation list etc, but for now
            // we've reduced the problem to sorting a whole reaction by an anchor node
            for(auto& [key, inst] : rule_matches)
            {
                //auto max_cell = anchored_phi(inst, model->system_graph, model->geoplex2D);
                auto max_cell = min_dim_phi(inst, model->system_graph, model->geoplex2D, component_matches);
                rule_map[max_cell].push_back(key);
            }

            auto sum = 0;
            for(auto& [key, value] : rule_map)
            {
                if(model->geoplex2D.getGraph().findNode(key)->second.getData().interior) {
                    auto dim = model->geoplex2D.getGraph().findNode(key)->second.getData().type;
                    std::cout << "Rules mapped to " << (2 - dim) << "D cell " << key << ": " << value.size() << "\n";
                }
                sum += value.size();
            }
            std::cout << "Total rules mapped " << sum << "\n";
        }

        //TODO: maybe not shared, but unique?
        std::shared_ptr<ModelType> model;

        CellList<graph_type> cell_list;

        AnalyzedGrammar<graph_type> grammar_analysis;

        KeyGenerator<std::size_t> instance_key_gen;

        //TODO: I think rule_maps could steal instances, especially since in this formulation,
        // each instance exists in a shared memory state, the rule_matches map.
        RuleMatchMap<typename ModelType::key_type>  rule_matches;
        ComponentMatchMap<typename ModelType::key_type> component_matches;

        std::map<cplex_key_t, std::vector<rule_key_t>> rule_map; //rule keys for a cell could be a tree/map/set vs vector
        std::map<std::size_t, std::vector<typename ModelType::key_type>> ordering;

        std::map<gplex_key_type, std::vector<key_type>> bucketsND[3];

        std::map<gplex_key_type, std::pair<double, double>> geocell_progress;
        DGGML::VtkFileWriter<typename ModelType::graph_type> vtk_writer;
        std::string results_dir_name;
    };
}

#endif //DGGML_SIMULATIONRUNNER_H
