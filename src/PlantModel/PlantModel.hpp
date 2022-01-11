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
            
            std::cout << "Generating the expanded cell complex\n";
            geoplex2D.init(1, 1, 25.0, 25.0, true); //ghosted
            std::cout << geoplex2D;
            
            //Save expanded cell complex graph
            Cajete::VtkFileWriter<typename Cajete::ExpandedComplex2D<>::types::graph_type> writer;
            writer.save(geoplex2D.getGraph(), "factory_geoplex");
            
            std::size_t num_mt = 8;
            std::cout << "Initializing the system graph\n";
            Plant::microtubule_unit_scatter(system_graph, geoplex2D, num_mt); 
            
            //std::cout << "Generating the grammar\n";
            //TODO: implement a grammar setup phase
    
        }

        void run() override {
            std::cout << "Running the plant model simulation\n";
            
            std::size_t num_steps = 25;
            Cajete::VtkFileWriter<graph_type> vtk_writer;
                
            std::string title = "results/factory_test_step_";
            std::cout << "Saving the initial state of the system graph\n";
            vtk_writer.save(system_graph, title+std::to_string(0));

            //TODO: move the simulation algorithm to its own class
            //is this where we run the simulation?
            for(auto i = 1; i <= num_steps; i++)
            {
                std::cout << "Running step " << i << std::endl;
                
                std::map<gplex_key_type, std::vector<key_type>> bucketsND[3];
                std::size_t complementND[3] = {0, 0, 0};

                
                std::cout << "Binning the graph into 2D partitions\n";
                //TODO: optimize this to work only for 2D
                Cajete::expanded_cartesian_complex_sort_stl(bucketsND, complementND, geoplex2D, system_graph);

                std::cout << "Running the Hybrid ODES/SSA inner loop 2D phase\n";
                //TODO: implement the inner loop of the SSA for 2D
                for(auto& bucket : bucketsND[0])
                {
                   auto k = bucket.first; //check to see if this is a domain to simulate
                   if(geoplex2D.getGraph().findNode(k)->second.getData().interior)
                   {
                        plant_model_ssa(bucket, geoplex2D, system_graph);
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
                //TODO: implement the inner loop of the SSA for 1D
                
                std::cout << "Synchronizing work\n";
                //TODO: this is where a barrier would be for a parallel code
                for(auto item : bucketsND) item.clear();

                std::cout << "Binning the graph into 0D partitions\n";
                //TODO: optimize this so it only works for 0D
                Cajete::expanded_cartesian_complex_sort_stl(bucketsND, complementND, geoplex2D, system_graph);

                std::cout << "Running the Hybrid ODES/SSA inner loop 0D phase\n";
                //TODO: implement the inner loop of the SSA for 0D
                
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
        CartesianComplex2D<> cplex2D;
        ExpandedComplex2D<> geoplex2D;
        YAGL::Graph<Plant::mt_key_type, Plant::MT_NodeData> system_graph; 
};

} //end namespace Cajete

#endif 
