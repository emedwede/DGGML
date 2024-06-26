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
#include "HelperStructs.hpp"

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

//            std::cout << "Printing grammar rules...\n";
//            model->gamma.print();

            std::cout << "Performing grammar analysis and building analyzed grammar data structure\n";
            grammar_analysis = AnalyzedGrammar<graph_type>(model->gamma);
            std::cout << "Grammar analysis complete\n";

            //building rule map sets
            for(auto& [key, value] : model->geoplex2D.graph.getNodeSetRef())
                rule_map.insert({key, {}});

            model->checkpoint(0);
            compute_single_component_matches();
            set_geocell_propensities();

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

            for(auto i = 1; i <= model->settings.NUM_STEPS; i++)
            {
                std::cout << "Running step " << i << "\n";

                double tot_time = 0.0;

                for(int d = 0; d < 3; d++)
                {
                    double dim_time = 0.0;

                    //TODO: for 1D and 0D, we may not have to re-bin the components into cells due to diffusion
                    // if we make the depth of reaction cells we search for 1D and 0D deeper i.e. search 2 away vs. 1
                    //if(d == 1 || d == 2) continue;
                    std::cout << "Running the Hybrid ODES/SSA inner loop " << (2 - d) << "D phase\n";
                    for(auto& bucket : bucketsND[d])
                    {
                        auto k = bucket.first;
                        auto start = std::chrono::high_resolution_clock::now();
                        approximate_ssa(component_matches, grammar_analysis,
                                        rule_map, rule_matches,
                                        model, k, geocell_properties_list[k], cell_list);
                        auto stop = std::chrono::high_resolution_clock::now();
                        auto duration =
                                std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
//                        std::cout << "Cell " << k << " took " << duration.count()
//                                  << " milliseconds and has a current tau "
//                                  << geocell_properties_list[k].tau << "\n";
                        dim_time += duration.count();
                        //break;//return;
                    }

                    //TODO: needs more work. First shot at adding back in rejected matches to lower dims
                    // It only works as is because of the type of rules we have, but rejected matches could still be rejected
                    // for a final time if they were invalidated by a rewrite after they were added to the list
                    //TODO: one partial fix: check to see if components in a rule still exist
                    for(auto& bucket : bucketsND[d])
                    {
                        auto k = bucket.first;
                        auto& geocell_properties = geocell_properties_list[k];
//                        geocell_properties.print();
                        //reset the rules fired counter
                        geocell_properties.current_rules_fired = 0;
                        if(!geocell_properties.rejected_rule_matches.empty())
                        {
                            for(auto& [cell_id, rejected_match] : geocell_properties.rejected_rule_matches)
                            {
                                //determine if any of the components of the match are not found, if yes, it's invalid
                                //since they were deleted by some other rule firing after the fact
                                bool valid = true;
                                for(auto& comp : rejected_match.components)
                                {
                                    if(component_matches.find(comp) == component_matches.end())
                                    {
                                       valid = false;
                                    }
                                }
                                if(valid) {
                                    auto res = rule_matches.insert(rejected_match);
                                    //if successfully inserted, added to the rule_map
                                    if (res.second) {
                                        rule_map[k].push_back(res.first->first);
                                    }
                                }
                            }
                            geocell_properties.rejected_rule_matches.clear();
                        }
                        //TODO: check to see if lower dimension have any rules containing invalidated components
                        ExpandedComplex2D<>& geoplex2D = model->geoplex2D;
//                        std::cout << "geocell " << k << " has nbrs: { ";
                        //for(auto& nbr : geoplex2D.graph.out_neighbors(k))
                        for(auto& [nbr, value] : model->geoplex2D.graph.getNodeSetRef())
                        {
//                            std::cout << nbr << " ";
                            std::vector<std::size_t> invalid_rids; //rule ids
                            std::vector<std::size_t> invalid_idxs; //indices
                            for(auto idx = 0; idx < rule_map[nbr].size(); idx++)
                            {
                                auto rid = rule_map[nbr][idx];
                                if(auto search = rule_matches.find(rid); search != rule_matches.end())
                                {
                                    for(auto& comp : search->second.components)
                                    {
                                        if(geocell_properties.invalidated_components.find(comp) != geocell_properties.invalidated_components.end())
                                        {

                                            //the rule is invalid
//                                            std::cout << "rule id " << rid << " in cell " << nbr << " contains component " << comp << "\n";

                                            invalid_rids.push_back(rid);
                                            invalid_idxs.push_back(idx);
                                            break; // we're done on the first of
                                        }
                                    }
                                }
                            }
                            //first sort the invalid_idxs in descending order since we will remove from the
                            //back first
                            auto cpy = rule_map[nbr];
                            std::sort(invalid_idxs.begin(), invalid_idxs.end(), std::greater<std::size_t>());
                            for(auto iter = 0; iter < invalid_rids.size(); iter++)
                            {
                                rule_matches.erase(invalid_rids[iter]);
                                //erase back to front
//                                auto id = std::find(rule_map[nbr].begin(), rule_map[nbr].end(), invalid_rids[iter]);
//                                if(id != rule_map[nbr].end())
//                                    rule_map[nbr].erase(id);
//                                if(rule_map.find(nbr) == rule_map.end())
//                                    std::cout << "nbr " << nbr << " does not exists\n";
//                                if(invalid_idxs[iter] >= rule_map[nbr].size()) {
//                                    std::cout << "index " << invalid_idxs[iter] << " dne\n";
//                                    std::cout << "max size is " << rule_map[nbr].size() << "\n";
//                                    std::cout << "originally: ";
//                                    std::cout << "max size was " << cpy.size() << "and we had rules: \n";
//                                    for(auto ii = 0; ii < cpy.size(); ii++)
//                                        std::cout << "{ " << ii << ", " << cpy[ii] << "} ";
//                                    std::cout << "\n";
//                                    std::cout << "we wanted to remove: ";
//                                    for(auto& item : invalid_idxs)
//                                        std::cout << item << " ";
//                                    std::cout << "\n";
//                                }
                                rule_map[nbr].erase(rule_map[nbr].begin()+invalid_idxs[iter]);
                            }
                        }
//                        std::cout << " }\n";
                        geocell_properties.invalidated_components.clear();
                    }
                    tot_time += dim_time;
                    std::cout << (2 - d) << "D took " << dim_time << " milliseconds\n";

                    //std::cout << "Synchronizing work\n";
                }
                model->checkpoint(i);
                std::cout << "Total dimensional time is " << tot_time << " milliseconds\n";
                //time_count.push_back(tot_time);
                std::cout << model->system_graph << "\n";

                //rebuild the cell list for the next iteration
                cell_list.rebuild();
                //clear the rule map if it's not empty
                if(!rule_map.empty())
                {
                    for(auto& item : rule_map)
                        item.second.clear();
                }
                //remap rules to geocells
                map_rule_matches_to_geocells();
                //std::cin.get();
                //if(i == 1)
                    //std::cin.get();
                //if(i == 3) return;
            }
            //model->print_metrics();
            //return;
        }
    //private:

        void set_geocell_propensities()
        {
            std::cout << "Setting intial cell propensities to zero\n";
            for(auto& [key, value] : model->geoplex2D.graph.getNodeSetRef())
                geocell_properties_list[key];
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
                    std::vector<typename ModelType::key_type> order;
                    std::vector<typename ModelType::key_type> match;
                    for(auto& [key, value] : item)
                    {
                        order.push_back(key);
                        match.push_back(value);
                    }
                    ComponentMatch<typename ModelType::key_type> inst;
                    inst.match = match;
                    inst.order = order;
                    inst.type = k;
                    inst.anchor = match[0];
                    component_matches.insert(inst);
                }
                //std::cout << "Found " << matches.size() << " instances\n";
            }

            //for(auto& [k, p] : instances)
//            for(auto& [k, pattern] : grammar_analysis.unique_components)
//            {
//                //std::cout << "Component " << k << " has " << component_matches.count(k) << " instances\n";
//            }
        }

        void compute_all_rule_matches()
        {
//            for(auto& [name, rule] : grammar_analysis.with_rules) {
//                auto count = std::count_if(rule_matches.begin(), rule_matches.end(), [&](auto& iter) { return iter.second.name == name; });
//                std::cout << "So far we have found " << count << " instances of with rule " << name << "\n";
//            }
//            for(auto& [name, rule] : grammar_analysis.solving_rules) {
//                auto count = std::count_if(rule_matches.begin(), rule_matches.end(), [&](auto& iter) { return iter.second.name == name; });
//                std::cout << "So far we have found " << count << " instances of solving rule " << name << "\n";
//            }
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

//            auto sum  = 0;
//            for(auto& name : grammar_analysis.rule_names) {
//                auto count = std::count_if(rule_matches.begin(), rule_matches.end(), [&](auto& iter) { return iter.second.name == name; });
//                std::cout << "Rule " << name << " has " << count << " instances\n";
//                sum += count;
//            }
//            std::cout << "There are " << sum << " rule instances in total\n";
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
                //TODO: fix this error where we get an invalid max cell for some reason
//                for(auto i = 0; i < inst.components.size(); i++)
//                {
//                    if(auto search = component_matches.find(inst.components[i]); search == component_matches.end())
//                    {
//                        std::cout << "did not find component " << i << "\n";
//                    }
//                    else if (auto& m = component_matches[inst.components[i]].match; m.size() == 0)
//                    {
//                        std::cout << "max_cell: " << max_cell << ", for rule id " << key << " has size " << inst.components.size() << "\n";
//                        std::cout << "found component " << inst.components[i] << " in an invalid state\n";
//                    }
//                }
                rule_map[max_cell].push_back(key);
            }

//            auto sum = 0;
//            for(auto& [key, value] : rule_map)
//            {
//                if(model->geoplex2D.getGraph().findNode(key)->second.getData().interior) {
//                    auto dim = model->geoplex2D.getGraph().findNode(key)->second.getData().type;
//                    std::cout << "Rules mapped to " << (2 - dim) << "D cell " << key << ": " << value.size() << "\n";
//                }
//                sum += value.size();
//            }
//            std::cout << "Total rules mapped " << sum << "\n";
        }

        //TODO: maybe not shared, but unique?
        std::shared_ptr<ModelType> model;

        CellList<graph_type> cell_list;

        AnalyzedGrammar<graph_type> grammar_analysis;

        KeyGenerator<std::size_t> instance_key_gen;

        //TODO: I think rule_maps could steal instances, especially since in this formulation,
        // each instance exists in a shared memory state, the rule_matches map.
        using node_key_type = typename ModelType::key_type;
        RuleMatchMap<node_key_type>  rule_matches;
        ComponentMatchMap<node_key_type> component_matches;

        std::map<cplex_key_t, std::vector<rule_key_t>> rule_map; //rule keys for a cell could be a tree/map/set vs vector
        std::map<std::size_t, std::vector<typename ModelType::key_type>> ordering;

        std::map<gplex_key_type, std::vector<key_type>> bucketsND[3];

        GeocellPropertiesList<gplex_key_type, node_key_type> geocell_properties_list;
        //DGGML::VtkFileWriter<typename ModelType::graph_type> vtk_writer;
    };
}

#endif //DGGML_SIMULATIONRUNNER_H
