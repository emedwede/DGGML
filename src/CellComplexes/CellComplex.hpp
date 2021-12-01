#ifndef CAJETE_CELL_COMPLEX_2D_HPP
#define CAJETE_CELL_COMPLEX_2D_HPP

#include "CellComplexTypes.hpp"
#include "BrickGrid2D.hpp"

namespace Cajete 
{
    template <typename GraphType = CplexGraph2D_t>
    class CellComplex2D
    {
        public:
            using graph_type = GraphType;
            using node_type = typename graph_type::node_type;
            using edge_type = typename graph_type::edge_type;
            using data_type = typename graph_type::data_type;
       
            CellComplex2D(BrickGrid2D& brick_grid, double fat_rad)
            {
                grid = &brick_grid;
                f_r = fat_rad;

                initailize();
            }

            void initailize()
            {
                fill_cell_complex2DBrick();
            }

            void fill_cell_complex2DBrick() 
            {
             
            }
        private:

            double f_r; //border fattening radius
            double eps; //criteria for well separatedness

            CplexGraph2D_t graph;
            BrickGrid2D* grid;
    };
} // end namespace Cajete

#endif
