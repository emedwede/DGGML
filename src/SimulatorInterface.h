//
// Created by Eric Medwedeff on 7/20/23.
//

#ifndef DGGML_SIMULATORINTERFACE_H
#define DGGML_SIMULATORINTERFACE_H

#include "SimulationRunner.h"

namespace DGGML {

    template<typename ModelType>
    class SimulatorInterface {
    public:
        void setModel(const ModelType& model)
        {
            simRunner = std::make_unique<SimulationRunner<ModelType>>(model);
        }

        void simulate()
        {
            if(simRunner)
            {
                simRunner->initialize();
                simRunner->run();
            }
            else
            {
                std::cout << "Error: no model set!";
            }
        }
    private:
        std::unique_ptr<SimulationRunner<ModelType>> simRunner;
    };

} // DGGML

#endif //DGGML_SIMULATORINTERFACE_H
