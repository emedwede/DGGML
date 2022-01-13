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

#include <map>
#include <random>

namespace Cajete
{
    struct Parameters
    {
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
    };
    
    template <typename ParamType>
    void set_parameters(ParamType& settings)
    {
        settings.DELTA = 0.1;
        settings.NUM_INTERNAL_STEPS = 10;
        settings.DELTA_DELTA_T = settings.DELTA / settings.NUM_INTERNAL_STEPS;
        
        settings.CELL_NX = 3;
        settings.CELL_NY = 3;
        
        settings.CELL_DX = 15.0;
        settings.CELL_DY = 15.0;
        settings.GHOSTED = true;

        settings.NUM_MT = 128;
        settings.MT_MIN_SEGMENT_INIT = 0.5;
        settings.MT_MAX_SEGMENT_INIT = 1.0;
        settings.NUM_STEPS = 50;

        settings.LENGTH_DIV_FACTOR = 1.2;
        settings.DIV_LENGTH = 2.0;
        settings.DIV_LENGTH_RETRACT = -0.2*settings.DIV_LENGTH;
        settings.V_PLUS = 1.0;
        settings.V_MINUS = settings.V_PLUS / 2.0;

        settings.SIGMOID_K = 10.0; 
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

        void init(InterfaceType interface) override {

            std::cout << "\n\n-----------------------------------------------------------------------\n";
            //TODO: implement timers to monitor start up phase
            std::cout << "Initializing the plant model simulation\n";
            
            std::cout << "Parsing the input interface and setting configuration settings\n";
            //TODO: handle the interface input
            set_parameters(settings); 

            std::cout << "Generating the expanded cell complex\n";
            geoplex2D.init(settings.CELL_NX, 
                    settings.CELL_NY, 
                    settings.CELL_DX, 
                    settings.CELL_DY, 
                    settings.GHOSTED); //ghosted
            std::cout << geoplex2D;
            
            //Save expanded cell complex graph
            Cajete::VtkFileWriter<typename Cajete::ExpandedComplex2D<>::types::graph_type> writer;
            writer.save(geoplex2D.getGraph(), "factory_geoplex");
            
            std::cout << "Initializing the system graph\n";
            Plant::microtubule_unit_scatter(system_graph, geoplex2D, settings); 
            
            //std::cout << "Generating the grammar\n";
            //TODO: implement a grammar setup phase
    
        }

        void run() override {
            std::cout << "Running the plant model simulation\n";
            
            Cajete::VtkFileWriter<graph_type> vtk_writer;
                
            std::string title = "results/factory_test_step_";
            std::cout << "Saving the initial state of the system graph\n";
            vtk_writer.save(system_graph, title+std::to_string(0));

            //TODO: move the simulation algorithm to its own class
            //is this where we run the simulation?
            for(auto i = 1; i <= settings.NUM_STEPS; i++)
            {
                std::cout << "Running step " << i << std::endl;
                
                std::map<gplex_key_type, std::vector<key_type>> bucketsND[3];
                std::size_t complementND[3] = {0, 0, 0};

                
                std::cout << "Binning the graph into 2D partitions\n";
                //TODO: optimize this to work only for 2D
                Cajete::expanded_cartesian_complex_sort_stl(bucketsND, complementND, geoplex2D, system_graph);

                std::cout << "Running the Hybrid ODES/SSA inner loop 2D phase\n";
                for(auto& bucket : bucketsND[0])
                {
                   auto k = bucket.first; //check to see if this is a domain to simulate
                   if(geoplex2D.getGraph().findNode(k)->second.getData().interior)
                   {
                        plant_model_ssa(bucket, geoplex2D, system_graph, settings);
                   }
                }
                
                //TODO: remove, right now connected_components should remain constant with only 
                //growth rules
                std::cout << "----------------\n";
                std::cout << "CC: " << YAGL::connected_components(system_graph); 
                std::cout << "\n---------------\n";
                
                std::cout << "Synchronizing work\n";
                //TODO: this is where a barrier would be for a parallel code
                for(auto& item : bucketsND) item.clear();
                

                std::cout << "Binning the graph into 1D partitions\n";
                //TODO: optimize this to work only for 1D
                Cajete::expanded_cartesian_complex_sort_stl(bucketsND, complementND, geoplex2D, system_graph);

                std::cout << "Running the Hybrid ODES/SSA inner loop 1D phase\n";
                for(auto& bucket : bucketsND[1])
                {
                   auto k = bucket.first; //check to see if this is a domain to simulate
                   if(geoplex2D.getGraph().findNode(k)->second.getData().interior)
                   {
                        std::cout << bucket.second.size() << " objects in bucket " << k << "\n";
                        plant_model_ssa(bucket, geoplex2D, system_graph, settings);
                   }
                }
               
                std::cout << "Synchronizing work\n";
                //TODO: this is where a barrier would be for a parallel code
                for(auto& item : bucketsND) item.clear();
                
                std::cout << "Binning the graph into 0D partitions\n";
                //TODO: optimize this so it only works for 0D
                Cajete::expanded_cartesian_complex_sort_stl(bucketsND, complementND, geoplex2D, system_graph);

                std::cout << "Running the Hybrid ODES/SSA inner loop 0D phase\n";
                for(auto& bucket : bucketsND[2])
                {
                   auto k = bucket.first; //check to see if this is a domain to simulate
                   if(geoplex2D.getGraph().findNode(k)->second.getData().interior)
                   {
                        plant_model_ssa(bucket, geoplex2D, system_graph, settings);
                   }
                }

                std::cout << "Synchronizing work\n";
                //TODO: this is where a barrier would be for a parallel code
                for(auto item : bucketsND) item.clear();
        
                std::cout << "Running the checkpointer\n";
                //TODO: The checkpointer to save time steps
                vtk_writer.save(system_graph, title+std::to_string(i));
            }
            std::cout << "-----------------------------------------------------------------------\n\n";

        }


    private:
        Parameters settings;
        CartesianComplex2D<> cplex2D;
        ExpandedComplex2D<> geoplex2D;
        YAGL::Graph<Plant::mt_key_type, Plant::MT_NodeData> system_graph; 
};

} //end namespace Cajete

#endif 
