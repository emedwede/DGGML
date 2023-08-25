#ifndef PLANT_MODEL_HPP
#define PLANT_MODEL_HPP

#include <iostream>

#include "PlantTypes.hpp"
#include "PlantUtils.hpp"
#include "PlantGrammar.hpp"
#include "PlantSSA.hpp"

#include "DggModel.hpp"
#include "YAGL_Graph.hpp"
#include "YAGL_Node.hpp"
#include "YAGL_Algorithms.hpp"
#include "Utlities/VtkWriter.hpp"

#include "ExpandedComplex2D.hpp"

#include "CartesianComplex2D.hpp"

#include "CartesianHashFunctions.hpp"

#include "RuleSystem.hpp"

#include "Utlities/MathUtils.hpp"

#include "CellList.hpp"

#include <map>
#include <random>
#include <chrono>
#include <string>
#include <filesystem>

namespace DGGML
{
    // Slowly deleting this class until it's seperated in the code base
    /*
}
    template <typename InterfaceType>
    class PlantModel : public DggModel<InterfaceType> {
    public:

        void run() override {
            //TODO: move the simulation algorithm to its own class
            //is this where we run the simulation?
            for(auto i = 1; i <= settings.NUM_STEPS; i++)
            {
                std::cout << "Running step " << i << std::endl;
                 
                // I think I only need to sort nodes by the highest dimensional cell to then
                // map rules to the approrpiate cell? 
                using cplex_key_t = typename DGGML::CartesianComplex2D<>::graph_type::key_type;
                //bucket of dimensional keys
                std::vector<cplex_key_t> bucket2d;
                std::vector<cplex_key_t> bucket1d;
                std::vector<cplex_key_t> bucket0d;

                for(auto& [key, value] : geoplex2D.graph.getNodeSetRef())
                {
                    if(value.getData().type == 0)
                        bucket2d.push_back(key);
                    else if(value.getData().type == 1)
                        bucket1d.push_back(key);
                    else if(value.getData().type == 2)
                        bucket0d.push_back(key);
                }
                std::cout << "Bucket2D: " << bucket2d.size() << "\n";
                std::cout << "Bucket1D: " << bucket1d.size() << "\n";
                std::cout << "Bucket0D: " << bucket0d.size() << "\n";

                using rule_key_t = std::size_t;             
                //first attempt at a dimensional mapping function
                std::map<cplex_key_t, std::vector<rule_key_t>> rule_map;
                for(const auto& item : bucket2d)
                    rule_map.insert({item, {}});
                for(const auto& item : bucket1d)
                    rule_map.insert({item, {}});
                for(const auto& item : bucket0d)
                    rule_map.insert({item, {}});
                
                //create a cell list
                CellList test_cell_list(geoplex2D.reaction_grid, system_graph, rule_system);
                auto growing_matches = 0;
                for(auto c = 0; c < test_cell_list.totalNumCells(); c++)
                {
                    int imin, imax, jmin, jmax;
                    test_cell_list.getCells(c, imin, imax, jmin, jmax); 
                    for(auto i = imin; i < imax; i++)
                    {
                        for(auto j = jmin; j < jmax; j++)
                        {
                            auto nbr_idx = test_cell_list.cardinalCellIndex(i, j);
                            for(const auto& match1 : test_cell_list.data[c])
                            {
                                for(const auto& match2 : test_cell_list.data[nbr_idx])
                                {
                                    if(match1.first.type == Rule::G && match2.first.type == Rule::G 
                                            && match1.second != match2.second)
                                    {
                                        auto a1 = match1.first.anchor;
                                        auto a2 = match2.first.anchor;
                                        auto& p1 = system_graph.findNode(a1)->second.getData().position;
                                        auto& p2 = system_graph.findNode(a2)->second.getData().position;
                                        auto d = calculate_distance(p1, p2); 
                                        //std::cout << "Distance: " << d << "\n";
                                        if(d < test_cell_list.grid._dx)
                                            growing_matches++;
                                    }
                                }
                            }
                        }
                    }
                }
                std::cout << "Potential Growing End Collision Matches Found: " << growing_matches << "\n"; 
                break;

                //reduces the potential multiset of anchor nodes into a set and
                //maps anchors to reaction subcells
                std::unordered_map<key_type, cplex_key_t> anchor_list;
                for(const auto& match : rule_system)
                {
                    auto& node_data = system_graph.findNode(match.first.anchor)->second.getData();

                    double xp = node_data.position[0];
                    double yp = node_data.position[1];
                    int ic, jc;
                    
                    geoplex2D.reaction_grid.locatePoint(xp, yp, ic, jc);
                    auto cardinal = geoplex2D.reaction_grid.cardinalCellIndex(ic, jc);

                    anchor_list.insert({match.first.anchor, cardinal});
                }
                std::cout << "Size of anchor_list: " << anchor_list.size() << "\n";
                
                std::vector<std::size_t> anchored_phi_count = {0, 0, 0, 0};
                for(const auto& match : rule_system)
                {
                    //current max not min, because dimensions are labeled as 2d -> 0 (TODO: fix)
                    auto& instance = match.first.match;
                    auto max_cardinal = anchor_list[match.first.anchor];
                    auto max_dim = geoplex2D.dim_label[max_cardinal];
                    auto max_cell = geoplex2D.cell_label[max_cardinal];
                    for(auto& l : instance)
                    {
                        auto& participation = rule_system.inverse_index.find(l)->second;
                        for(auto& p : participation)
                        {
                            auto cur_cardinal = anchor_list[rule_system[p].anchor];
                            auto cur_dim = geoplex2D.dim_label[cur_cardinal];
                            if(cur_dim > max_dim)
                            {
                                max_dim = cur_dim;
                                max_cell = geoplex2D.cell_label[cur_cardinal];
                            }
                        }
                    }
                    rule_map[max_cell].push_back(match.second);
                    anchored_phi_count[max_dim]++;
                }

                auto d_tot = 0;
                for(auto i = 0; i < anchored_phi_count.size(); i++) 
                {
                    std::cout << "Dim " << i << ": " << anchored_phi_count[i] << "\n";
                    d_tot += anchored_phi_count[i];
                }
                std::cout << "Total mapped: " << d_tot << "\n";

                //rule map total 
                auto r_tot = 0;
                for(const auto& [k, v] : rule_map)
                {
                    r_tot += v.size();
                    auto dim = geoplex2D.getGraph().findNode(k)->second.getData().type;
                }
                std::cout << "Total rules: " << r_tot << "\n";

                //scoped printing section for testing
                //TODO: make these scoped prints unit tests for phi mapping function
                {std::unordered_map<std::size_t, std::size_t> counts;
                for(auto& c : bucket2d) counts.insert({c, 0});
                for(auto& [key, value] : cell_list) counts[value]++;
                auto fold_op = 
                    [](auto val, auto& item){ return val+item.second; };
                auto bucket_sum = 
                    std::accumulate(counts.begin(), counts.end(), 0, fold_op);
                std::cout << bucket_sum << " nodes sorted\n";}

                //scoped block print, every rule should be assigned
                {auto sum = 0;
                for(auto& item : rule_map) sum += item.second.size();
                std::cout << "Rules assigned: " << sum << "\n";
                std::cout << "Rules existing: " << rule_system.size() << "\n";}
                std::map<gplex_key_type, std::vector<key_type>> bucketsND[3];
                std::size_t complementND[3] = {0, 0, 0};
                
                DGGML::build_bucketsND(bucketsND, geoplex2D);
                
                double tot_time = 0.0;

                for(int i = 0; i < 3; i++)
                {
                    double dim_time = 0.0;
                    
                    //if(i == 0 || i == 1) continue;
                    std::cout << "Running the Hybrid ODES/SSA inner loop " << (2 - i) << "D phase\n";
                    for(auto& bucket : bucketsND[i])
                    {
                        auto k = bucket.first; 
                        auto start = std::chrono::high_resolution_clock::now();
                        geocell_progress[k] = 
                            plant_model_ssa_new(rule_system, k, rule_map, anchor_list, 
                                    geoplex2D, system_graph, settings, geocell_progress[k]);
                        auto stop = std::chrono::high_resolution_clock::now();
                        auto duration = 
                        std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
                        std::cout << "Cell " << k << " took " << duration.count() 
                            << " milliseconds and has a current tau " 
                            << geocell_progress[k].first << "\n";
                        dim_time += duration.count();
                    }
                
                    tot_time += dim_time;
                    std::cout << (2 - i) << "D took " << dim_time << " milliseconds\n";
                                
                    std::cout << "Synchronizing work\n";
                }

                //TODO: this is where a barrier would be for a parallel code
                std::cout << "Running the checkpointer\n";
                //TODO: The checkpointer to save time steps
                vtk_writer.save(system_graph, title+std::to_string(i));
                
                std::cout << "Total dimensional time is " << tot_time << " milliseconds\n";
                time_count.push_back(tot_time);
            }
        }


    private:
        RuleSystem<Plant::mt_key_type> rule_system;
        Parameters settings;
        CartesianComplex2D<> cplex2D;
};
*/
} //end namespace DGGML

#endif 
