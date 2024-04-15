#ifndef DGGML_CMAMODEL_H
#define DGGML_CMAMODEL_H
#include <fstream>

#include "DGGML.h"
#include "PlantUtils.hpp"
#include "cmaRules.hpp"
#include "MathUtils.hpp"
#include "simdjson.h"
#include "ExpandedComplex2D.hpp"
#include "YAGL_Algorithms.hpp"
#include "Parameters.hpp"

namespace CMA {
    using graph_grammar_t = DGGML::Grammar<Plant::graph_type>;
    class cmaModel : public DGGML::Model<graph_grammar_t> {
    public:
        void initialize() override
        {
            std::cout << "Initializing the plant model simulation\n";
            //settings.set_default();
            name = settings.EXPERIMENT_NAME;
            std::cout << "Creating the grammar\n";

            if(settings.ENABLE_GROWTH)
                create_with_mt_growing_rule(gamma, system_graph,settings);

            if(settings.ENABLE_GROWTH)
                create_ode_mt_growing_rule(gamma, system_graph,settings);

            if(settings.ENABLE_RETRACTION)
                create_with_mt_retraction_rule(gamma, system_graph,settings);

            if(settings.ENABLE_RETRACTION)
                create_ode_mt_retraction_rule(gamma, system_graph,settings);

            if(settings.ENABLE_STANDARD_BOUNDARY_CATASTROPHE)
                create_with_standard_boundary_catastrophe(gamma, system_graph,settings);

            if(settings.ENABLE_CLASP_BOUNDARY_CATASTROPHE)
                create_with_clasp_boundary_catastrophe(gamma, system_graph,settings);

            if(settings.ENABLE_INTERMEDIATE_CIC)
                create_with_intermediate_cic(gamma, system_graph,settings);

            if(settings.ENABLE_NEGATIVE_CIC)
                create_with_negative_cic(gamma, system_graph,settings);

            if(settings.ENABLE_POSITIVE_CIC)
                create_with_positive_cic(gamma, system_graph,settings);

            if(settings.ENABLE_CROSSOVER)
                create_with_crossover_rule(gamma, system_graph, settings);

            if(settings.ENABLE_UNCROSSOVER)
                create_with_uncrossover_rule(gamma, system_graph, settings);

            if(settings.ENABLE_ZIPPERING)
                create_with_zippering_rules(gamma, system_graph, settings);

            if(settings.CLASP_ENABLE_ENTRY)
                create_with_clasp_entry_rule(gamma, system_graph, settings);

            if(settings.CLASP_ENABLE_EXIT)
                create_with_clasp_exit_rule(gamma, system_graph, settings);

            if(settings.CLASP_ENABLE_DETACHMENT)
                create_with_clasp_detachment(gamma, system_graph, settings);

            if(settings.CLASP_ENABLE_CAT)
                create_with_clasp_catastrophe(gamma, system_graph, settings);

            if(settings.ENABLE_MT_DESTRUCTION)
                create_with_destruction_rule(gamma, system_graph, settings);

            if(settings.ENABLE_CREATION)
                create_with_creation_rule(gamma, system_graph, settings);

            if(settings.ENABLE_RECOVERY)
                create_with_recovery_rule(gamma, system_graph, settings);

            std::cout << "Generating the expanded cell complex\n";
            geoplex2D.init(settings.CELL_NX,
                           settings.CELL_NY,
                           settings.CELL_DX,
                           settings.CELL_DY,
                           settings.GHOSTED,
                           settings.MAXIMAL_REACTION_RADIUS); //ghosted
            //std::cout << geoplex2D;
            std::cout << "Initializing the system graph\n";
            Plant::microtubule_uniform_scatter(system_graph, geoplex2D, settings, gen);
        }

        struct MyMetrics {
            std::vector<std::size_t> con_com;
            std::vector<std::size_t> total_nodes;
            std::vector<std::size_t> type_counts[5];
            std::vector<double> time_count;
            std::vector<double> correlation_avg_global;
            std::vector<double> correlation_avg_local;
            std::vector<double> correlation_dist;
            std::vector<std::size_t> angular_histogram;

            std::size_t junction_count;
            std::size_t positive_count;
            std::size_t negative_count;
            std::size_t zipper_count;
            std::size_t intermediate_count;

            void reset_count()
            {
                junction_count = 0;
                positive_count = 0;
                negative_count = 0;
                zipper_count = 0;
                intermediate_count = 0;
                correlation_dist.clear();
                angular_histogram.clear();

            }
        } metrics;

        void collect(std::size_t step) override
        {
            metrics.reset_count();
            metrics.con_com.push_back(YAGL::connected_components(system_graph));

            for(auto iter = system_graph.node_list_begin();
                iter != system_graph.node_list_end(); iter++) {
                auto& itype = iter->second.getData();
                if(std::holds_alternative<Plant::Negative>(itype.data))
                    metrics.negative_count++;
                if(std::holds_alternative<Plant::Positive>(itype.data))
                    metrics.positive_count++;
                if(std::holds_alternative<Plant::Intermediate>(itype.data))
                    metrics.intermediate_count++;
                if(std::holds_alternative<Plant::Junction>(itype.data))
                    metrics.junction_count++;
                if(std::holds_alternative<Plant::Zipper>(itype.data))
                    metrics.zipper_count++;
            }
            metrics.type_counts[0].push_back(metrics.negative_count);
            metrics.type_counts[1].push_back(metrics.positive_count);
            metrics.type_counts[2].push_back(metrics.intermediate_count);
            metrics.type_counts[3].push_back(metrics.junction_count);
            metrics.type_counts[4].push_back(metrics.zipper_count);

            metrics.total_nodes.push_back(metrics.negative_count + metrics.positive_count +
            metrics.intermediate_count + metrics.junction_count + metrics.zipper_count);
            metrics.correlation_dist = compute_two_point_correlation_alpha(system_graph, settings);
            auto size = (double)metrics.correlation_dist.size();
            auto avg = std::accumulate(metrics.correlation_dist.begin(), metrics.correlation_dist.end(), 0.0)/size;
            metrics.correlation_avg_global.push_back(avg);
            auto local_size = metrics.correlation_dist.size()/3;
            auto local_avg =
                    std::accumulate(metrics.correlation_dist.begin(), metrics.correlation_dist.begin()+local_size, 0.0)
                    /((double)local_size);
            metrics.correlation_avg_local.push_back(local_avg);
            metrics.angular_histogram = compute_orientation_histogram(system_graph, settings);
        }

        void print_metrics(std::size_t step) override
        {
            // Open the file for appending, create if it doesn't exist
            std::string filename = settings.RESULTS_DIR+"/"+"metrics" + std::to_string(step) +".json";
            std::ofstream outputFile(filename, std::ios::app | std::ios::out);
            outputFile << "{";
            print_json_array_stats(outputFile, metrics.con_com, "con_com");
            outputFile << ",\n";
            print_json_array_stats(outputFile, metrics.type_counts[0], "negative");
            outputFile << ",\n";
            print_json_array_stats(outputFile, metrics.type_counts[1], "positive");
            outputFile << ",\n";
            print_json_array_stats(outputFile, metrics.type_counts[2], "intermediate");
            outputFile << ",\n";
            print_json_array_stats(outputFile, metrics.type_counts[3], "junction");
            outputFile << ",\n";
            print_json_array_stats(outputFile, metrics.type_counts[4], "zipper");
            outputFile << ",\n";
            print_json_array_stats(outputFile, metrics.total_nodes, "total_nodes");
            outputFile << ",\n";
            print_json_array_stats(outputFile, metrics.correlation_avg_global, "correlation_avg_global");
            outputFile << ",\n";
            print_json_array_stats(outputFile, metrics.correlation_avg_local, "correlation_avg_local");
            outputFile << ",\n";
            print_json_array_stats(outputFile, metrics.correlation_dist, "correlation_dist");
            outputFile << ",\n";
            print_json_array_stats(outputFile, metrics.angular_histogram, "angular_hist");
            outputFile << "\n}";
            //print_numpy_array_stats_csv(outputFile, metrics.time_count, "time_count");
            outputFile.close();
        }

        // Checkpoint uses basic file writing functionality provided by the VTK file writer
        void checkpoint(std::size_t step) override
        {
            auto results_dir_name = settings.RESULTS_DIR;
            // on the first step of the simulation
            if(step == 0)
            {
                std::cout << "Running the checkpoint for step " << step << "\n";
                // Create the local save directory
                std::cout << "Cleaning up old results folder if it exists and creating a new one\n";
                std::filesystem::remove_all(results_dir_name);
                std::filesystem::create_directory(results_dir_name);

                // Save expanded cell complex graph
                //DGGML::VtkFileWriter<typename DGGML::ExpandedComplex2D<>::types::graph_type> writer;
                //writer.save(model->geoplex2D.getGraph(), results_dir_name+"/factory_geoplex");
                DGGML::GridFileWriter grid_writer;
                grid_writer.save({geoplex2D.reaction_grid,geoplex2D.dim_label},
                                 results_dir_name+"/expanded_cell_complex");

                std::string title = results_dir_name+"/simulation_step_";
                std::cout << "Saving the initial state of the system graph\n";
                DGGML::VtkFileWriter<graph_type> vtk_writer;
                vtk_writer.save(system_graph, title+std::to_string(step));
                collect(step);
            }
            //save every 10 steps
            if(step != 0 && step % settings.CHECKPOINT_FREQUENCY == 0 || step == settings.NUM_STEPS)
            {
                std::cout << "Running the checkpoint for step " << step << "\n";
                std::cout << "Saving the system graph\n";
                std::string title = results_dir_name+"/simulation_step_";
                DGGML::VtkFileWriter<graph_type> vtk_writer;
                vtk_writer.save(system_graph, title+std::to_string(step));
                std::cout << "Collecting metrics\n";
                collect(step);
                std::cout << "Writing the collected metrics to a file\n";
                print_metrics(step);
            }
            //std::cin.get();
        }

        void set_parameters(simdjson::ondemand::document& interface)
        {
            settings.set_parameters(interface);
        }

        //TODO: separate core settings out into the base class
        Parameters settings;
    };
}


#endif //DGGML_CMAMODEL_H
