#ifndef DGG_FACTORY_HPP
#define DGG_FACTORY_HPP

#include <iostream>
#include "Utlities/MemoryManager.hpp"

namespace DGGML
{

	//Very simple templated factory pattern interface

	template <template<typename> typename ModelType, typename InterfaceType>
	class DggFactory {
    public:
        using model_type = ModelType<InterfaceType>;
        using interface_type = InterfaceType;

        ~DggFactory() {}
        
        model_type* create(interface_type& interface) const {
            //allocate the model using my memory manager wrapper
            return DGGML::MemoryManager::allocate_std<model_type>(1);
            //return new model_type(); 
        }

        void execute(interface_type& interface) const {
            model_type* model = this->create(interface);

            model->init(interface);

            model->run();
            
            //deallocate using my memory manager wrapper
            DGGML::MemoryManager::deallocate_std(model);
            //delete model;
        }
	};

} //end namespace DGGML

#endif
