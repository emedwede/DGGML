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
        using key_type = Plant::mt_key_type; //TODO: Plant namespace needs to be removed from here
        using gplex_key_type = typename DGGML::ExpandedComplex2D<>::graph_type::key_type;
        using data_type = Plant::MT_NodeData;
        using graph_type = YAGL::Graph<key_type, data_type>;
        using node_type = typename graph_type::node_type;
        explicit SimulationRunner(const ModelType& model) : model(std::make_shared<ModelType>(model)) {}

        //May want to make this part of the constructor instead
        void initialize()
        {
            //Should the model initialize even happen here? or sooner
            std::cout << "Initializing " << model->name << "\n";
            model->initialize();

            //order matters here, which indicates maybe I should have a
            //file writer class which initializes with the save directory?
            create_save_directory();
            write_cell_complex();

            set_geocell_propensities();

            model->gamma.print();
            //compute_matches();
        }

        void run() {

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

//        void compute_matches()
//        {
//            std::vector<std::vector<Plant::mt_key_type>> match_set;
//            //TODO: I think we need to store the ordering or the rooted spanning tree
//            for(auto& pattern : model->gamma.minimal_set)
//            {
//                auto matches = YAGL::subgraph_isomorphism2(pattern, model->system_graph);
//                for(auto& item : matches)
//                {
//                    std::vector<Plant::mt_key_type> match;
//                    for(auto& [key, value] : item)
//                    {
//                        match.push_back(value);
//                        std::cout << "{" << key << " -> " << value << "} ";
//                    } std::cout << "\n";
//                    rule_system.push_back({std::move(match), DGGML::Rule::G});
//                }
//                std::cout << "Found " << matches.size() << " instances\n";
//            }
//        }
        std::shared_ptr<ModelType> model;
        //Grammar gamma;
        RuleSystem<Plant::mt_key_type> rule_system;
        std::map<gplex_key_type, std::pair<double, double>> geocell_progress;
        std::string results_dir_name;
    };
}

#endif //DGGML_SIMULATIONRUNNER_H
