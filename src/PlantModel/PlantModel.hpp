#ifndef PLANT_MODEL_HPP
#define PLANT_MODEL_HPP

#include <iostream>

#include "PlantTypes.hpp"

#include "DggModel.hpp"
#include "YAGL_Graph.hpp"
#include "YAGL_Node.hpp"
#include "VtkWriter.hpp"

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
            std::cout << "Initializing the plant model simulation\n";
            std::random_device random_device;
            std::mt19937 random_engine(random_device());
            std::uniform_real_distribution<double> distribution_global(1.0, 9.0);
            std::uniform_real_distribution<double> distribution_local(0.5, 1.0);
   
            std::size_t num_mt = 8;
            std::size_t segments = 3;
            for(auto i = 0; i < num_mt; i++) 
            {
                auto x_c = distribution_global(random_engine);
                auto y_c = distribution_global(random_engine);
                auto z_c = 0.0; //distribution_global(random_engine);

                auto x_r = x_c + distribution_local(random_engine);
                auto y_r = y_c + distribution_local(random_engine); 
                auto z_r = 0.0; //z_c + distribution_local(random_engine);
                
                auto x_l = x_c - distribution_local(random_engine);
                auto y_l = y_c - distribution_local(random_engine); 
                auto z_l = 0.0; //z_c - distribution_local(random_engine);
                
                node_type node_l(i*segments, {{x_l, y_l, z_l}, {0.0, 0.0, 0.0}, Plant::negative});
                node_type node_c(i*segments+1, {{x_c, y_c, z_c}, {0.0, 0.0, 0.0}, Plant::intermediate});
                node_type node_r(i*segments+2, {{x_r, y_r, z_r}, {0.0, 0.0, 0.0}, Plant::positive});

                system_graph.addNode(node_l);
                system_graph.addNode(node_c);
                system_graph.addNode(node_r);

                system_graph.addEdge(node_l, node_c);
                system_graph.addEdge(node_r, node_c);
            }
        }

        void run() override {
            std::cout << "Running the plant model simulation\n";
            
            std::size_t num_steps = 100;

            //is this where we run the simulation?
            for(auto i = 0; i < num_steps; i++)
            {
                std::cout << "Running step " << i << std::endl;
                Cajete::VtkFileWriter<graph_type> vtk_writer;
                vtk_writer.save(system_graph, "factory_test");
            }
        }


    private:
        YAGL::Graph<Plant::mt_key_type, Plant::MT_NodeData> system_graph; 
};

} //end namespace Cajete

#endif 
