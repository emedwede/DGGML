#ifndef CAJETE_CELL_COMPLEX_HPP
#define CAJETE_CELL_COMPLEX_HPP

namespace Cajete 
{
    template <typename GraphType>
    class CellComplex
    {
        public:
            using graph_type = GraphType;
            using node_type = typename graph_type::node_type;
            using edge_type = typename graph_type::edge_type;
            using data_type = typename graph_type::data_type;

        private:
            GraphType graph;
    };
} // end namespace Cajete

#endif
