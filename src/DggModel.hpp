#ifndef DGG_MODEL_HPP
#define DGG_MODEL_HPP

namespace DGGML
{

//Simple specification for a DGG Model base type
//
//It's template on a generic interface type, this 
//is inteded to ensure the model can be set by any
//command line parser
template <typename InterfaceType>
class DggModel {
    public:
        using interface_type = InterfaceType;

        virtual ~DggModel() = default;

        virtual void init(InterfaceType& interface) = 0;

        virtual void run() = 0;
};

} //end namespace DGGML

#endif 
