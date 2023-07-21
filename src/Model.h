#ifndef DGGML_MODEL_H
#define DGGML_MODEL_H

namespace DGGML {
    template<typename GraphGrammarType>
    class Model {
    public:
        virtual void initialize() = 0;

        virtual ~Model() = default;

        std::string name;
        GraphGrammarType gamma;
    };
}

#endif //DGGML_MODEL_H
