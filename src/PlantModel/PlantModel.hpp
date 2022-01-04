#ifndef PLANT_MODEL_HPP
#define PLANT_MODEL_HPP

#include <iostream>

#include "PlantTypes.hpp"
#include "PlantUtils.hpp"

#include "DggModel.hpp"
#include "YAGL_Graph.hpp"
#include "YAGL_Node.hpp"
#include "VtkWriter.hpp"

#include "ExpandedComplex2D.hpp"

#include "CartesianComplex2D.hpp"

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
            geoplex2D.init(1, 1, 15.0, 15.0, true); //ghosted
            std::cout << geoplex2D;
            
            //Save expanded cell complex graph
            Cajete::VtkFileWriter<typename Cajete::ExpandedComplex2D<>::types::graph_type> writer;
            writer.save(geoplex2D.getGraph(), "factory_geoplex");

            std::cout << "Initializing the system graph\n";
            Plant::microtubule_unit_scatter(system_graph, geoplex2D, 5); 
            
            std::cout << "Generating the grammar\n";
            //TODO: implement a grammar setup phase
    
        }

        void run() override {
            std::cout << "Running the plant model simulation\n";
            
            std::size_t num_steps = 1;
            Cajete::VtkFileWriter<graph_type> vtk_writer;
                

            std::cout << "Saving the initial state of the system graph\n";
            vtk_writer.save(system_graph, "factory_test_step_0");

            //TODO: move the simulation algorithm to its own class
            //is this where we run the simulation?
            for(auto i = 1; i <= num_steps; i++)
            {
                std::cout << "Running step " << i << std::endl;

                std::cout << "Binning the graph into 2D partitions\n";
                //TODO: implement a 2D partioning scheme using the cell complex
                                
                std::cout << "Running the Hybrid ODES/SSA inner loop 2D phase\n";
                //TODO: implement the inner loop of the SSA for 2D
                
                std::cout << "Synchronizing work\n";
                //TODO: this is where a barrier would be for a parallel code
                
                std::cout << "Binning the graph into 1D partitions\n";
                //TODO: implement a 1D partioning scheme using the cell complex
                
                std::cout << "Running the Hybrid ODES/SSA inner loop 1D phase\n";
                //TODO: implement the inner loop of the SSA for 1D
                
                std::cout << "Synchronizing work\n";
                //TODO: this is where a barrier would be for a parallel code

                std::cout << "Binning the graph into 0D partitions\n";
                //TODO: implement a 2D partioning scheme using the cell complex
                
                std::cout << "Running the Hybrid ODES/SSA inner loop 0D phase\n";
                //TODO: implement the inner loop of the SSA for 0D
                
                std::cout << "Synchronizing work\n";
                //TODO: this is where a barrier would be for a parallel code
                
                std::cout << "Running the checkpointer\n";
                //TODO: The checkpointer to save time steps
                //vtk_writer.save(system_graph, "factory_test");
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
