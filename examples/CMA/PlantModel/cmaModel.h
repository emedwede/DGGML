#ifndef DGGML_CMAMODEL_H
#define DGGML_CMAMODEL_H

#include "DGGML.h"

namespace CMA {
    using graph_grammar_t = DGGML::Grammar;
    class cmaModel : public DGGML::Model<graph_grammar_t> {
    public:
        void initialize() const override
        {
            std::cout << "Initializing model\n";
        }
    };
}


#endif //DGGML_CMAMODEL_H
