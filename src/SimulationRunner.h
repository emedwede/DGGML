#ifndef DGGML_SIMULATIONRUNNER_H
#define DGGML_SIMULATIONRUNNER_H

#include <iostream>
#include <memory>

namespace DGGML {
    template<typename ModelType>
    class SimulationRunner {
    public:
        explicit SimulationRunner(const ModelType& model) : model(std::make_shared<ModelType>(model)) {}

        void run() {
            if(model)
                std::cout << "Running " << model->name << "\n";
            else
                std::cout << "Model not set!\n";
        }
    private:
        std::shared_ptr<ModelType> model;
    };
}

#endif //DGGML_SIMULATIONRUNNER_H
