#ifndef DGGML_CMAMODEL_H
#define DGGML_CMAMODEL_H
#include <fstream>

#include "DGGML.h"
#include "PlantUtils.hpp"
#include "cmaRules.hpp"
#include "MathUtils.hpp"
//#include "simdjson.h"
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
            settings.set_default();
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
            std::vector<std::vector<double>> angular_correlation;

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
            }
        } metrics;

        void collect() override
        {
            metrics.reset_count();
            metrics.con_com.push_back(YAGL::connected_components(system_graph));

            for(auto iter = system_graph.node_list_begin();
                iter != system_graph.node_list_end(); iter++) {
                auto itype = iter->second.getData();
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
            //metrics.angular_correlation.push_back(compute_two_point_correlation_alpha(system_graph, settings));
        }

        void print_metrics() override
        {
            // Open the file for appending, create if it doesn't exist
            std::string filename = settings.RESULTS_DIR+"/"+"metrics.csv";
            std::ofstream outputFile(filename, std::ios::app | std::ios::out);
            print_numpy_array_stats_csv(outputFile, metrics.con_com, "con_com");
            print_numpy_array_stats_csv(outputFile, metrics.type_counts[0], "negative");
            print_numpy_array_stats_csv(outputFile, metrics.type_counts[1], "positive");
            print_numpy_array_stats_csv(outputFile, metrics.type_counts[2], "intermediate");
            print_numpy_array_stats_csv(outputFile, metrics.type_counts[3], "junction");
            print_numpy_array_stats_csv(outputFile, metrics.type_counts[4], "zipper");
            print_numpy_array_stats_csv(outputFile, metrics.total_nodes, "total_nodes");
            auto start = 0;
            auto end = metrics.angular_correlation.size()-1;
            //print_numpy_array_stats_csv(outputFile, metrics.angular_correlation[start], "angle_corr_start");
            //print_numpy_array_stats_csv(outputFile, metrics.angular_correlation[end], "angle_corr_end");
            auto result = compute_two_point_correlation_alpha(system_graph, settings);
            print_numpy_array_stats_csv(outputFile, result, "angle_corr_end");
            auto hist = compute_orientation_histogram(system_graph, settings);
            print_numpy_array_stats_csv(outputFile, hist, "orientation_hist");
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
                collect();
            }
            //save every 10 steps
            if(step != 0 && step % settings.CHECKPOINT_FREQUENCY == 0)
            {
                std::cout << "Running the checkpoint for step " << step << "\n";
                std::cout << "Saving the system graph\n";
                std::string title = results_dir_name+"/simulation_step_";
                DGGML::VtkFileWriter<graph_type> vtk_writer;
                vtk_writer.save(system_graph, title+std::to_string(step));
                collect();
            }
            if(step == settings.NUM_STEPS)
            {
                std::cout << "Writing the final metrics to a file\n";
                print_metrics();
            }
            //std::cin.get();
        }

//        void set_parameters(simdjson::ondemand::document& interface)
//        {
//            settings.set_parameters(interface);
//        }

        //TODO: separate core settings out into the base class
        Parameters settings;
    };
}


#endif //DGGML_CMAMODEL_H
