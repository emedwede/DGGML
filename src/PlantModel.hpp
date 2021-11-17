#ifndef PLANT_MODEL_HPP
#define PLANT_MODEL_HPP

#include <iostream>

#include "DggModel.hpp"

namespace Cajete
{
    //Models are inteded to be designed based on the 
    //DggModel specification. Right now it's very loose
    //and capable of handling almost anything
    template <typename InterfaceType>
    class PlantModel : public DggModel<InterfaceType> {
    public:
        void init(InterfaceType interface) const override {
            std::cout << "Initializing the plant model simulation\n";
        }

        void run() const override {
            std::cout << "Running the plant model simulation\n";
        }
};

} //end namespace Cajete

#endif 
