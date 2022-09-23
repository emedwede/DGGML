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
#include "VtkWriter.hpp"

#include "ExpandedComplex2D.hpp"

#include "CartesianComplex2D.hpp"

#include "CartesianHashFunctions.hpp"

#include "RuleSystem.hpp"

#include <map>
#include <random>
#include <chrono>
#include <string>
#include <filesystem>

namespace Cajete
{
    struct Parameters
    {
        std::string EXPERIMENT_NAME;
        double DELTA;
        double DELTA_DELTA_T;
        int NUM_INTERNAL_STEPS;
        std::size_t CELL_NX;
        std::size_t CELL_NY;
        double CELL_DX;
        double CELL_DY;
        bool GHOSTED;
        std::size_t NUM_MT;
        double MT_MIN_SEGMENT_INIT;
        double MT_MAX_SEGMENT_INIT;
        std::size_t NUM_STEPS;
        double LENGTH_DIV_FACTOR;
        double DIV_LENGTH;
        double DIV_LENGTH_RETRACT;
        double V_PLUS;
        double V_MINUS;
        double SIGMOID_K;
        double TOTAL_TIME;
        double MAXIMAL_REACTION_RADIUS;
        double DELTA_T_MIN;
        double RHO_TEST_RATE; //a tunable test paramater for MT dynamics
    };
  
    template <typename DataType>
    void print_numpy_array_stats(DataType& data, std::string var_name)
    {
        std::cout << var_name << " = np.asarray([";
            for(auto i = 0; i < data.size(); i++)
            {
                if(i != 0 && i % 20 == 0) std::cout << "\n";
                if(i != data.size() - 1)
                    std::cout << data[i] << ", ";
                else
                    std::cout << data[i];
            }
            std::cout << "]);\n";
    }

    template <typename ParamType, typename InterfaceType>
    void set_parameters(ParamType& settings, InterfaceType& interface)
    {
        std::string_view temp = interface["META"]["EXPERIMENT"];
        settings.EXPERIMENT_NAME = static_cast<std::string>(temp);

        std::cout << settings.EXPERIMENT_NAME << "+++\n";
                settings.CELL_NX = std::size_t(interface["SETTINGS"]["CELL_NX"]); 
        settings.CELL_NY = std::size_t(interface["SETTINGS"]["CELL_NY"]);
        
        settings.CELL_DX = double(interface["SETTINGS"]["CELL_DX"]);
        settings.CELL_DY = double(interface["SETTINGS"]["CELL_DY"]);
        settings.GHOSTED = bool(interface["SETTINGS"]["GHOSTED"]);

        settings.NUM_MT = std::size_t(interface["SETTINGS"]["NUM_MT"]);
        settings.MT_MIN_SEGMENT_INIT = double(interface["SETTINGS"]["MT_MIN_SEGMENT_INIT"]);
        settings.MT_MAX_SEGMENT_INIT = double(interface["SETTINGS"]["MT_MAX_SEGMENT_INIT"]);

        settings.LENGTH_DIV_FACTOR = double(interface["SETTINGS"]["LENGTH_DIV_FACTOR"]);
        settings.DIV_LENGTH = double(interface["SETTINGS"]["DIV_LENGTH"]);
        settings.DIV_LENGTH_RETRACT = double(interface["SETTINGS"]["DIV_LENGTH_RETRACT"]);
        settings.V_PLUS = double(interface["SETTINGS"]["V_PLUS"]);
        settings.V_MINUS = double(interface["SETTINGS"]["V_MINUS"]);

        settings.SIGMOID_K = double(interface["SETTINGS"]["SIGMOID_K"]);
        
        settings.MAXIMAL_REACTION_RADIUS = settings.DIV_LENGTH*2.0;
        
        //Simulate until the specified unit time
        settings.TOTAL_TIME = double(interface["SETTINGS"]["TOTAL_TIME"]);
        settings.NUM_INTERNAL_STEPS = 5;
        //Delta should be big, but not to big. In this case, the maximum amount of time it would
        //take one MT to grow a single unit of MT
        settings.DELTA = 
            0.25*settings.MAXIMAL_REACTION_RADIUS / std::max(settings.V_PLUS, settings.V_MINUS);
        //The internal step of the solver should be at least this small
        settings.DELTA_DELTA_T = settings.DELTA / settings.NUM_INTERNAL_STEPS; 
        settings.DELTA_T_MIN = settings.DELTA_DELTA_T;
        settings.NUM_STEPS = settings.TOTAL_TIME / settings.DELTA;

        settings.RHO_TEST_RATE = double(interface["EXPERIMENTAL"]["RHO_TEST_RATE"]);
    }

    //Models are inteded to be designed based on the 
    //DggModel specification. Right now it's very loose
    //and capable of handling almost anything
    template <typename InterfaceType>
    class PlantModel : public DggModel<InterfaceType> {
    public:
        using key_type = Plant::mt_key_type;
        using gplex_key_type = typename Cajete::ExpandedComplex2D<>::graph_type::key_type;
        using data_type = Plant::MT_NodeData;
        using graph_type = YAGL::Graph<key_type, data_type>;
        using node_type = typename graph_type::node_type;

        void init(InterfaceType& interface) override {

            std::cout << "\n\n-----------------------------------------------------------------------\n";
            //TODO: implement timers to monitor start up phase
            std::cout << "Initializing the plant model simulation\n";
            
            std::cout << "Parsing the input interface and setting configuration settings\n";
            //TODO: handle the interface input
            set_parameters(settings, interface); 
            
            std::cout << "Cleaning up old results folder if it exists and creating a new one\n";
            results_dir_name = settings.EXPERIMENT_NAME + "_results";
            std::filesystem::remove_all(results_dir_name);
            std::filesystem::create_directory(results_dir_name);
            
            std::cout << "Generating the expanded cell complex\n";
            geoplex2D.init(settings.CELL_NX, 
                    settings.CELL_NY, 
                    settings.CELL_DX, 
                    settings.CELL_DY, 
                    settings.GHOSTED); //ghosted
            //std::cout << geoplex2D;
           
            std::cout << "Setting intial cell propensities to zero\n";
            for(auto& [key, value] : geoplex2D.graph.getNodeSetRef())
                geocell_progress[key] = {0.0, 0.0};

            //Save expanded cell complex graph
            Cajete::VtkFileWriter<typename Cajete::ExpandedComplex2D<>::types::graph_type> writer;
            writer.save(geoplex2D.getGraph(), results_dir_name+"/factory_geoplex");
            
            std::cout << "Initializing the system graph\n";
            Plant::microtubule_uniform_scatter(system_graph, geoplex2D, settings); 
            
            //std::cout << "Generating the grammar\n";
            //TODO: implement a grammar setup phase
            
            //TODO: we should fold the below code into some interface 
            //function to apply the rule set
            
            //precomputes all of the single component matches//
            
            //bucket to temporarily interfact the matcher rule
            std::vector<Cajete::Plant::mt_key_type> node_set;
            for(const auto& [key, value] : system_graph.getNodeSetRef())
                node_set.push_back(key);
            
            auto matches = 
                Cajete::Plant::microtubule_growing_end_matcher(system_graph, node_set);
            
            for(auto& item : matches)
                rule_system.push_back({std::move(item), Cajete::Rule::G});
            
            matches = 
                Cajete::Plant::microtubule_retraction_end_matcher(system_graph, node_set);
            
            for(auto& item : matches)
                rule_system.push_back({std::move(item), Cajete::Rule::R});
            
            std::cout << "matches: " << rule_system.size() << "\n";
        }

        void run() override {
            std::cout << "Running the plant model simulation\n";
            
            Cajete::VtkFileWriter<graph_type> vtk_writer;
            std::vector<std::size_t> con_com;
            con_com.push_back(YAGL::connected_components(system_graph));
            std::vector<std::size_t> total_nodes;
            std::vector<std::size_t> type_counts[5];
            std::vector<double> time_count; 
            std::size_t junction_count = 0;
            std::size_t positive_count = 0;
            std::size_t negative_count = 0;
            std::size_t zipper_count = 0;
            std::size_t intermediate_count = 0;
            for(auto iter = system_graph.node_list_begin(); 
                    iter != system_graph.node_list_end(); iter++) {
                    auto itype = iter->second.getData().type;
                    if(itype == Plant::negative)
                        negative_count++;
                    if(itype == Plant::positive)
                        positive_count++;
                    if(itype == Plant::intermediate)
                        intermediate_count++;
                    if(itype == Plant::junction)
                        junction_count++;
                    if(itype == Plant::zipper)
                        zipper_count++;
            }
            type_counts[Plant::negative].push_back(negative_count);
            type_counts[Plant::positive].push_back(positive_count);
            type_counts[Plant::intermediate].push_back(intermediate_count);
            type_counts[Plant::junction].push_back(junction_count);
            type_counts[Plant::zipper].push_back(zipper_count);
            
            total_nodes.push_back(negative_count +
                    positive_count + intermediate_count + junction_count + zipper_count);
                
            std::string title = results_dir_name+"/simulation_step_";
            std::cout << "Saving the initial state of the system graph\n";
            vtk_writer.save(system_graph, title+std::to_string(0));
            
            //TODO: move the simulation algorithm to its own class
            //is this where we run the simulation?
            for(auto i = 1; i <= settings.NUM_STEPS; i++)
            {
                std::cout << "Running step " << i << std::endl;
                 
                // I think I only need to sort nodes by the highest dimensional cell to then
                // map rules to the approrpiate cell? 
                using cplex_key_t = typename Cajete::CartesianComplex2D<>::graph_type::key_type;
                //bucket of dimensional keys
                std::vector<cplex_key_t> bucket2d;
                std::vector<cplex_key_t> bucket1d;

                for(auto& [key, value] : geoplex2D.graph.getNodeSetRef())
                {
                    if(value.getData().type == 0)
                        bucket2d.push_back(key);
                    else if(value.getData().type == 1)
                        bucket1d.push_back(key);
                }
                std::cout << "Bucket2D: " << bucket2d.size() << "\n";
                std::cout << "Bucket1D: " << bucket1d.size() << "\n";
                std::unordered_map<key_type, gplex_key_type> cell_list;
                
                for(auto& [key, value] : system_graph.getNodeSetRef())
                {
                    auto& node_data = value.getData();
                    double xp = node_data.position[0];
                    double yp = node_data.position[1];
                    int ic, jc;
                    
                    geoplex2D.coarse_grid.locatePoint(xp, yp, ic, jc);
                    geoplex2D.coarse_cell_to_fine_lattice(ic, jc); 
                    auto cardinal = geoplex2D.fine_grid.cardinalLatticeIndex(ic, jc);
                    cell_list.insert({key, cardinal}); 
                }
                
                //scoped printing section for testing
                {std::unordered_map<std::size_t, std::size_t> counts;
                for(auto& c : bucket2d) counts.insert({c, 0});
                for(auto& [key, value] : cell_list) counts[value]++;
                auto fold_op = 
                    [](auto val, auto& item){ return val+item.second; };
                auto bucket_sum = 
                    std::accumulate(counts.begin(), counts.end(), 0, fold_op);
                std::cout << bucket_sum << " nodes sorted\n";}

                using rule_key_t = std::size_t;
                std::map<cplex_key_t, std::vector<rule_key_t>> rule_map;
                for(const auto& item : bucket2d)
                    rule_map.insert({item, {}});
                for(const auto& item : bucket1d)
                    rule_map.insert({item, {}});
                
                for(const auto& match : rule_system)
                {
                    auto& instance = match.first.match;
                    //prints ids
                    //for(auto& l : instance)
                    //    std::cout << cell_list.find(l)->second << " ";
                    //std::cout << "\n";
                    std::set<Plant::mt_key_type> intersection;
                    for(auto& l : instance)
                        intersection.insert(cell_list[l]);
                    if(intersection.size() == 1)
                    {
                        rule_map[*intersection.begin()].push_back(match.second);
                        std::cout << "maps to highest dim cell\n";
                    }
                    else
                    {
                        std::cout << "intersection size: " 
                            << intersection.size() << "\n";
                        std::set<std::size_t> edge_set;
                        std::set<std::size_t> vert_set;
                        std::set<std::size_t> found_set;
                        auto& geoplex_graph = geoplex2D.getGraph();
                        for(auto& item : intersection)
                        {
                            for(auto iter = geoplex_graph.out_neighbors_begin(item);
                                    iter != geoplex_graph.out_neighbors_end(item);
                                    iter++)
                            {
                                auto id = *iter;
                                auto inode = geoplex_graph.findNode(id)->second;
                                auto inode_data = inode.getData();
                                if(inode_data.type == 1 && inode_data.interior)
                                {
                                    auto search = found_set.find(id);
                                    if(search != found_set.end())
                                        edge_set.insert(id);
                                    else 
                                        found_set.insert(id);

                                    for(auto jter = 
                                            geoplex_graph.out_neighbors_begin(id);
                                            jter != 
                                            geoplex_graph.out_neighbors_end(id);
                                            jter++)
                                    {
                                        auto jd = *jter;
                                        auto jnode = 
                                            geoplex_graph.findNode(jd)->second;
                                        auto jnode_data = jnode.getData();
                                        if(jnode_data.type == 2 && jnode_data.interior)
                                        {
                                            auto search = found_set.find(jd);
                                            if(search != found_set.end())
                                                vert_set.insert(jd);
                                            else
                                                found_set.insert(jd);
                                        }
                                    }
                                }
                            }
                        }
                        if(edge_set.size() == 0 && vert_set.size() == 1)
                        {
                            std::cout << "maps to 0D cell\n";
                            rule_map[*vert_set.begin()].push_back(match.second); 
                        }
                        else if(edge_set.size() == 1)
                        {
                            std::cout << "maps to 1D cell\n";
                            rule_map[*edge_set.begin()].push_back(match.second);
                        }
                        else
                            std::cout << "Currently undefined scenario\n";
                        //else if(edge_set.size() > 1)
                        //    std::cout << "maps to 0D cell\n";
                        
                    }
                }
                
                //scoped block print, every rule should be assigned
                {auto sum = 0;
                for(auto& item : rule_map) sum += item.second.size();
                std::cout << "Rules assigned: " << sum << "\n";
                std::cout << "Rules existing: " << rule_system.size() << "\n";}
                std::map<gplex_key_type, std::vector<key_type>> bucketsND[3];
                std::size_t complementND[3] = {0, 0, 0};
                
                double dim_time = 0.0;
                double tot_time = 0.0;
                
                Cajete::build_bucketsND(bucketsND, geoplex2D);
                
                std::cout << "Running the Hybrid ODES/SSA inner loop 2D phase\n";
                for(auto& bucket : bucketsND[0])
                {
                    
                    auto k = bucket.first; 
                    auto start = std::chrono::high_resolution_clock::now();
                    geocell_progress[k] = 
                        plant_model_ssa_new(rule_system, k, rule_map, cell_list, 
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
                std::cout << "2D took " << dim_time << " milliseconds\n";
                
                //std::cout << "----------------\n";
                //std::cout << "CC: " << YAGL::connected_components(system_graph); 
                //std::cout << "\n---------------\n";
                                
                std::cout << "Synchronizing work\n";
                
                dim_time = 0.0;
                std::cout << "Running the Hybrid ODES/SSA inner loop 1D phase\n";
                for(auto& bucket : bucketsND[1])
                {
                    auto k = bucket.first;
                    auto start = std::chrono::high_resolution_clock::now();
                    geocell_progress[k] = 
                        plant_model_ssa_new(rule_system, k, rule_map, cell_list, 
                            geoplex2D, system_graph, settings, geocell_progress[k]);
                    auto stop = std::chrono::high_resolution_clock::now();
                    auto duration = 
                    std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
                    std::cout << "Cell " << k << " took " << 
                        duration.count() << " milliseconds\n";
                    dim_time += duration.count();
                }
                
                tot_time += dim_time;
                std::cout << "1D took " << dim_time << " milliseconds\n";

                std::cout << "Synchronizing work\n";
                //TODO: this is where a barrier would be for a parallel code
                for(auto& item : bucketsND) item.clear();
                
                dim_time = 0.0;
                std::cout << "Running the Hybrid ODES/SSA inner loop 0D phase\n";
                for(auto& bucket : bucketsND[2])
                {
                    auto k = bucket.first; 
                    auto start = std::chrono::high_resolution_clock::now();
                    geocell_progress[k] = 
                        plant_model_ssa_new(rule_system, k, rule_map, cell_list, 
                            geoplex2D, system_graph, settings, geocell_progress[k]);
                    //plant_model_ssa(bucket, geoplex2D, system_graph, settings);
                    auto stop = std::chrono::high_resolution_clock::now();
                    auto duration = 
                    std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
                    std::cout << "Cell " << k << " took " 
                        << duration.count() << " milliseconds\n";
                    dim_time += duration.count();
                }
                tot_time += dim_time;
                std::cout << "0D took " << dim_time << " milliseconds\n";
        
                std::cout << "Synchronizing work\n";
                

                //TODO: this is where a barrier would be for a parallel code
                std::cout << "Running the checkpointer\n";
                //TODO: The checkpointer to save time steps
                vtk_writer.save(system_graph, title+std::to_string(i));
                
                con_com.push_back(YAGL::connected_components(system_graph)); 
                
                std::size_t junction_count = 0;
                std::size_t positive_count = 0;
                std::size_t negative_count = 0;
                std::size_t zipper_count = 0;
                std::size_t intermediate_count = 0;
                for(auto iter = system_graph.node_list_begin(); 
                        iter != system_graph.node_list_end(); iter++) {
                        auto itype = iter->second.getData().type;
                        if(itype == Plant::negative)
                            negative_count++;
                        if(itype == Plant::positive)
                            positive_count++;
                        if(itype == Plant::intermediate)
                            intermediate_count++;
                        if(itype == Plant::junction)
                            junction_count++;
                        if(itype == Plant::zipper)
                            zipper_count++;
                }
                type_counts[Plant::negative].push_back(negative_count);
                type_counts[Plant::positive].push_back(positive_count);
                type_counts[Plant::intermediate].push_back(intermediate_count);
                type_counts[Plant::junction].push_back(junction_count);
                type_counts[Plant::zipper].push_back(zipper_count);
                
                total_nodes.push_back(negative_count +
                        positive_count + intermediate_count + junction_count + zipper_count);
                
                std::cout << "Total dimensional time is " << tot_time << " milliseconds\n";
                time_count.push_back(tot_time);
            }
            std::cout << "-----------------------------------------------------------------------\n\n";
            
            /*
            print_numpy_array_stats(con_com, "con_com");
            print_numpy_array_stats(type_counts[Plant::negative], "negative");
            print_numpy_array_stats(type_counts[Plant::positive], "positive");
            print_numpy_array_stats(type_counts[Plant::intermediate], "intermediate");
            print_numpy_array_stats(type_counts[Plant::junction], "junction");
            print_numpy_array_stats(type_counts[Plant::zipper], "zipper");
            print_numpy_array_stats(total_nodes, "total_nodes");
            print_numpy_array_stats(time_count, "time_count");
            */
        }


    private:
        RuleSystem<Plant::mt_key_type> rule_system;
        std::map<gplex_key_type, std::pair<double, double>> geocell_progress;
        Parameters settings;
        CartesianComplex2D<> cplex2D;
        ExpandedComplex2D<> geoplex2D;
        YAGL::Graph<Plant::mt_key_type, Plant::MT_NodeData> system_graph; 
        std::string results_dir_name;
};

} //end namespace Cajete

#endif 
