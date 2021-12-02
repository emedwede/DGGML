#ifndef CAJETE_CARTESIAN_COMPLEX_2D_HPP
#define CAJETE_CARTESIAN_COMPLEX_2D_HPP

#include "CellComplexTypes.hpp"
#include "CartesianGrid2D.hpp"

namespace Cajete 
{
    template <typename GraphType = CplexGraph2D_t>
    class CartesianComplex2D
    {
        public:
            using graph_type = GraphType;
            using node_type = typename graph_type::node_type;
            using edge_type = typename graph_type::edge_type;
            using data_type = typename graph_type::data_type;
            
            CartesianComplex2D(std::size_t n, std::size_t m, double d_x, double d_y)
            {
                ppc = 2; //Fixed for 2D lattice

                nx = n*ppc+1; 
                ny = m*ppc+1; 
                dx = d_x; 
                dy = d_y;
                
                grid.init(0.0, 0.0, n*dx, m*dy, nx, ny);

                build();
            }

                        const CplexGraph2D_t& getGraph() 
            {
                return graph;
            }

            friend std::ostream& operator<<(std::ostream& os, CartesianComplex2D& cplex)
            {
                //os << cplex.grid;
                os << cplex.graph;
                
                return os;
            }
            
        private:
            
            void build()
            {
                std::cout << "Building the cell complex for a 2D cartesian grid\n";
                
                //We do extra work because redundant edges may be attempted to be added
                for(int j = 1; j < ny; j += 2)
                {
                    for(int i = 1; i < nx; i += 2)
                    {
                        // Generate all the edges from center of 2D cell to center of 1D cell
                        create_edge(i, j, 0, i, j+1, 1);
                        create_edge(i, j, 0, i, j-1, 1);
                        create_edge(i, j, 0, i+1, j, 1);
                        create_edge(i, j, 0, i-1, j, 1);

                        // Generate all the edge from center of 1D cell to center of 0D cells 
                        create_edge(i+1, j, 1, i+1, j+1, 2);
                        create_edge(i+1, j, 1, i+1, j-1, 2);

                        create_edge(i-1, j, 1, i-1, j+1, 2);
                        create_edge(i-1, j, 1, i-1, j-1, 2);
                        
                        create_edge(i, j-1, 1, i+1, j-1, 2);
                        create_edge(i, j-1, 1, i-1, j-1, 2);
                        
                        create_edge(i, j+1, 1, i+1, j+1, 2);
                        create_edge(i, j+1, 1, i-1, j+1, 2);
                    }
                }
            }
            

            void create_edge(int i_a, int j_a, std::size_t a_t, int i_b, int j_b, std::size_t b_t)
            {
                        std::size_t a = grid.cardinalCellIndex(i_a, j_a);
                        std::size_t b = grid.cardinalCellIndex(i_b, j_b);

                        double px_a, py_a, pz_a, px_b, py_b, pz_b;
                        pz_a = pz_b = 0.0; //set the z components to be zero in 2D
                        
                        //convert the cardinal indices to grid coordinates
                        grid.cardinalToPoint(px_a, py_a, a);
                        grid.cardinalToPoint(px_b, py_b, b);
                        
                        //Add the 2D cell if it does not exist, otherwise assign
                        graph.addNode({a, {a_t, {px_a, py_a, pz_a}, {{0, 0}, {0, 0}, {0, 0}, {0, 0}}}});
                        //Add the 1D cell if it does not exist, otherwise assign 
                        graph.addNode({b, {b_t, {px_b, py_b, pz_b}, {{0, 0}, {0, 0}, {0, 0}, {0, 0}}}});
                        
                        graph.addEdge(graph.findNode(a)->second, graph.findNode(b)->second);
            }

            double f_r; //border fattening radius
            double eps; //criteria for well separatedness
            
            std::size_t num_2D_cells;
            std::size_t num_1D_cells;
            std::size_t num_0D_cells;
            std::size_t total_cells;
            std::size_t nx; //number of cells in x direction
            std::size_t ny; //number of cells in y direction 
            
            std::size_t ppc; //number of discrete points to form the subcells

            double dx; //width of cell in the x direction 
            double dy; //width of cell in the y direction 
            
            CplexGraph2D_t graph;
            CartesianGrid2D grid;
    };

    } // end namespace Cajete

#endif
