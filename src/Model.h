#ifndef DGGML_MODEL_H
#define DGGML_MODEL_H

namespace DGGML {
    template<typename GraphGrammarType>
    class Model {
    public:
        std::string name;
        GraphGrammarType gamma;

        using graph_grammar_type = GraphGrammarType;
        using graph_type = typename GraphGrammarType::graph_type;
        using key_type = typename GraphGrammarType::key_type;

        KeyGenerator<std::size_t> gen;
        graph_type system_graph;
        ExpandedComplex2D<> geoplex2D;

        virtual void initialize() = 0;

        virtual void collect(std::size_t step) { std::cout << "Base collect does nothing by default\n"; };
        virtual void print_metrics(std::size_t step) { std::cout << "Base print_metrics does nothing by default\n"; }
        virtual void checkpoint(std::size_t step) { std::cout << "Base checkpointing does nothing by default\n"; }

        virtual ~Model() = default;

    };
}

#endif //DGGML_MODEL_H
