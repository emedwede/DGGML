#ifndef DGGML_CARTESIAN_COMPLEX_2D_HPP
#define DGGML_CARTESIAN_COMPLEX_2D_HPP

#include "CellComplexTypes.hpp"
#include "CartesianGrid2D.hpp"

namespace DGGML
{
    template <typename GraphType = CplexGraph2D_t>
    class CartesianComplex2D
    {
        public:
            using graph_type = GraphType;
            using node_type = typename graph_type::node_type;
            using edge_type = typename graph_type::edge_type;
            using data_type = typename graph_type::data_type;
            
            CartesianComplex2D() = default;

            CartesianComplex2D(std::size_t n, std::size_t m, double d_x, double d_y, bool g_c = false)
            {
                init(n, m, d_x, d_y, g_c);
            }
            
            //TODO: fix incompatibility between grid parameters and cell complex cell size
            void init(std::size_t n, std::size_t m, double d_x, double d_y, bool g_c = false)
            {
                ppc = 2; //Fixed for 2D lattice
                 
                ghost_cells = g_c; //indicates that we need to account for ghost cells
                if(g_c)
                {
                    n = n+2; //pad around boundary in dim  
                    m = m+2; //pad around boundary in dim
                } 

                nx = n*ppc;//+1; 
                ny = m*ppc;//+1; 

                dx = d_x; 
                dy = d_y;

                fine_grid.init(0.0, 0.0, n*dx, m*dy, nx, ny);
                coarse_grid.init(0.0, 0.0, n*dx, m*dy, n, m);

                nx = fine_grid._px;
                ny = fine_grid._py;

                if(ghost_cells) 
                {
                    min_x = coarse_grid._min_x + coarse_grid._dx;
                    min_y = coarse_grid._min_y + coarse_grid._dy;
                    max_x = coarse_grid._max_x - coarse_grid._dx;
                    max_y = coarse_grid._max_y - coarse_grid._dy;
                } else
                {
                    min_x = coarse_grid._min_x;
                    min_y = coarse_grid._min_y;
                    max_x = coarse_grid._max_x;
                    max_y = coarse_grid._max_y;
                }

                num_interior_0D_cells = 0;
                num_exterior_0D_cells = 0;
                num_interior_1D_cells = 0;
                num_exterior_1D_cells = 0;
                num_interior_2D_cells = 0;
                num_exterior_2D_cells = 0;

                build();
            }

            const CplexGraph2D_t& getGraph() 
            {
                return graph;
            }
            
            const CartesianGrid2D& getFineGrid()
            {
                return fine_grid;
            }
            
            const CartesianGrid2D& getCoarseGrid()
            {
                return coarse_grid;
            }
            
            //converts the coarse cell ids to the the fine lattice point ids
            void coarse_cell_to_fine_lattice(int& ic, int& jc)
            {
                ic = 2*ic+1;
                jc = 2*jc+1;
            }

            std::size_t get0dInteriorCellCount()
            {
                return num_interior_0D_cells;
            }

            std::size_t get0dExteriorCellCount()
            {
                return num_exterior_0D_cells;
            }
            
            std::size_t get0dTotalCellCount()
            {
                return num_interior_0D_cells + num_exterior_0D_cells;
            }

            std::size_t get1dInteriorCellCount()
            {
                return num_interior_1D_cells;
            }
            
            std::size_t get1dExteriorCellCount()
            {
                return num_exterior_1D_cells;
            }
            
            std::size_t get1dTotalCellCount()
            {
                return num_interior_1D_cells + num_exterior_1D_cells;
            }

            std::size_t get2dInteriorCellCount()
            {
                return num_interior_2D_cells;
            }

            std::size_t get2dExteriorCellCount()
            {
                return num_exterior_2D_cells;
            }
            
            std::size_t get2dTotalCellCount()
            {
                return num_interior_2D_cells + num_exterior_2D_cells;
            }
             
            std::size_t getTotalInteriorCellCount()
            {
                return get0dInteriorCellCount() + get1dInteriorCellCount() + get2dInteriorCellCount();
            }

            std::size_t getTotalExteriorCellCount()
            {
                return get0dExteriorCellCount() + get1dExteriorCellCount() + get2dExteriorCellCount();
            }

            std::size_t getTotalCellCount() 
            {
                return get0dTotalCellCount() + get1dTotalCellCount() + get2dTotalCellCount();
            }
            
            //TODO: generate lists of interior cells 
            
            friend std::ostream& operator<<(std::ostream& os, CartesianComplex2D& cplex)
            {
                os << "\nCoarse Grid: \n" << cplex.coarse_grid;
                os << "\nFine Grid: \n" << cplex.fine_grid;
                os << cplex.graph;
                
                return os;
            }
             
        public:
            
            void build()
            {
                //std::cout << "Building the cell complex for a 2D cartesian grid\n";
                
                //We do extra work because redundant edges may be attempted to be added
                for(int j = 1; j < fine_grid._py; j += 2)
                {
                    for(int i = 1; i < fine_grid._px; i += 2)
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

                //next compute the number of interior edges
                //there are two versions one is for a ghosted complex, and that is much easier
                //to compute
                for(auto iter = graph.node_list_begin(); iter != graph.node_list_end(); iter++)
                {
                    auto key = iter->first;
                    auto& data = iter->second.getData(); //needs to be auto& since auto doesn't deduce a reference type
                
                    auto num_2D_nbrs = 0;
                    //std::cout << "Corners: for point {" << data.position[0] << ", " << data.position[1] << "}:\n";
                    //for(auto i = 0; i < 4; i++) 
                    //{
                    //    std::cout << "{" << data.corners[i][0] << ", " << data.corners[i][1] << "}\n";
                    //}

                    if(data.type == 0) // found a 2D cell
                    {
                        if(!ghost_cells) // if not ghost cells method 1, count all 2D 
                        {
                            num_interior_2D_cells++;
                            data.interior = true;
                        } else { //else method 2
                            if(data.ghosted == true) //if this node is ghosted then it in the exterior 
                            {
                                num_exterior_2D_cells++;
                                data.interior = false;
                            } else {
                                num_interior_2D_cells++;
                                data.interior = true;
                            }        
                        }
                    }


                    if(data.type == 1) //found an 1D cell type 
                    {
                        if(!ghost_cells)
                        {
                            //we now need to iterate over its connections and
                            //count how many 2D cells it is connected to,
                            //alternatievely, we could just check if it's nbr list
                            //size is 3. If it is then it's an exterior. If it's 4
                            //then it is interior. The first way is more generic
                            for(auto j = graph.out_neighbors_begin(iter->second);
                                j != graph.out_neighbors_end(iter->second); j++)
                            {
                                if(graph.findNode(*j)->second.getData().type == 0)
                                {
                                    num_2D_nbrs++;
                                }
                            } 

                            if(num_2D_nbrs > 1)
                            {
                                num_interior_1D_cells++;
                                data.interior = true;
                            } else //we found an exterior 1D cell
                            {
                                num_exterior_1D_cells++;
                                data.interior = false;
                            }
                        } else {
                            if(data.ghosted == true)
                            {
                                num_exterior_1D_cells++;
                                data.interior = false;
                            } else {
                                num_interior_1D_cells++;
                                data.interior = true;
                            }    
                        }
                    }
                    
                    auto num_1D_nbrs = 0;

                    if(data.type == 2)
                    {
                        if(!ghost_cells)
                        {
                            for(auto j = graph.out_neighbors_begin(iter->second);
                                j != graph.out_neighbors_end(iter->second); j++)
                            {
                                if(graph.findNode(*j)->second.getData().type == 1)
                                {
                                    num_1D_nbrs++;
                                }
                            }
                            if(num_1D_nbrs == 4) {
                                num_interior_0D_cells++;
                                data.interior = true;
                            } else {
                                num_exterior_0D_cells++;
                                data.interior = false;
                            }
                        } else {
                            if(data.ghosted == true)
                            {
                                num_exterior_0D_cells++;
                                data.interior = false;
                            } else {
                                num_interior_0D_cells++;
                                data.interior = true;
                            }
                        }
                    }
                }
            }
            

            void create_edge(int i_a, int j_a, std::size_t a_t, int i_b, int j_b, std::size_t b_t)
            {
                        std::size_t a = fine_grid.cardinalLatticeIndex(i_a, j_a);
                        std::size_t b = fine_grid.cardinalLatticeIndex(i_b, j_b);

                        double px_a, py_a, pz_a, px_b, py_b, pz_b;
                        pz_a = pz_b = 0.0; //set the z components to be zero in 2D
                        
                        //convert the cardinal indices to fine_grid coordinates
                        fine_grid.cardinalLatticeToPoint(px_a, py_a, a);
                        fine_grid.cardinalLatticeToPoint(px_b, py_b, b);
                        
                        //ghosted
                        bool ghosted_a = false;
                        bool ghosted_b = false;
                        bool interior_a = true;
                        bool interior_b = true;
                        if(ghost_cells)
                        {
                            auto w = ppc + 1;
                            //in the outer envelope a.k.a. the ghosted region
                            if((i_a < w || i_a >= nx - w) || (j_a < w || j_a >= ny - w))
                            {
                                ghosted_a = true;
                                interior_a = false;
                            }
                            if((i_b < w || i_b >= nx - w) || (j_b < w || j_b >= ny - w))
                            {
                                ghosted_b = true;
                                interior_b = false;
                            }

                        }
                        
                        if(a_t == 0)
                        {
                            double c0x, c0y, c1x, c1y, c2x, c2y, c3x, c3y;

                            //lower left
                            fine_grid.cardinalLatticeToPoint(c0x, c0y, fine_grid.cardinalLatticeIndex(i_a-1, j_a-1));
                            //lower right
                            fine_grid.cardinalLatticeToPoint(c1x, c1y, fine_grid.cardinalLatticeIndex(i_a+1, j_a-1));
                            //upper left
                            fine_grid.cardinalLatticeToPoint(c2x, c2y, fine_grid.cardinalLatticeIndex(i_a-1, j_a+1));
                            //upper right
                            fine_grid.cardinalLatticeToPoint(c3x, c3y, fine_grid.cardinalLatticeIndex(i_a+1, j_a+1));

                            //Add the cell if it does not exist, otherwise assign
                            graph.addNode({a, {a_t, {px_a, py_a, pz_a}, 
                                {{c0x, c0y}, {c1x, c1y}, {c2x, c2y}, {c3x, c3y}}, ghosted_a, interior_a}});
                        } 
                        else //corners are not expanded
                        {
                            //Add the cell if it does not exist, otherwise assign
                            graph.addNode({a, {a_t, {px_a, py_a, pz_a}, 
                                {{px_a, py_a}, {px_a, py_a}, {px_a, py_a}, {px_a, py_a}}, ghosted_a, interior_a}});
                        }

                        if(b_t == 0) 
                        {
                            double c0x, c0y, c1x, c1y, c2x, c2y, c3x, c3y;
                            fine_grid.cardinalLatticeToPoint(c0x, c0y, fine_grid.cardinalLatticeIndex(i_b-1, j_b-1));
                            fine_grid.cardinalLatticeToPoint(c1x, c1y, fine_grid.cardinalLatticeIndex(i_b+1, j_b-1));
                            fine_grid.cardinalLatticeToPoint(c2x, c2y, fine_grid.cardinalLatticeIndex(i_b-1, j_b+1));
                            fine_grid.cardinalLatticeToPoint(c3x, c3y, fine_grid.cardinalLatticeIndex(i_b+1, j_b+1));

                            //Add the cell if it does not exist, otherwise assign
                            graph.addNode({b, {b_t, {px_b, py_b, pz_b}, 
                                {{c0x, c0y}, {c1x, c1y}, {c2x, c2y}, {c3x, c3y}}, ghosted_b, interior_b}});
                        } 
                        else //corners are not expanded
                        {
                            //Add the cell if it does not exist, otherwise assign 
                            graph.addNode({b, {b_t, {px_b, py_b, pz_b}, 
                                {{px_b, py_b}, {px_b, py_b}, {px_b, py_b}, {px_b, py_b}}, ghosted_b, interior_b}});
                        
                        }
                        //Add the cell if it does not exist, otherwise assign
                        //graph.addNode({a, {a_t, {px_a, py_a, pz_a}, 
                        //        {{0, 0}, {0, 0}, {0, 0}, {0, 0}}, ghosted_a}});
                        //Add the cell if it does not exist, otherwise assign 
                        //graph.addNode({b, {b_t, {px_b, py_b, pz_b}, 
                        //        {{0, 0}, {0, 0}, {0, 0}, {0, 0}}, ghosted_b}});
                        
                        graph.addEdge(a, b);
            }

            bool ghost_cells;
            double min_x, min_y, max_x, max_y;
            std::size_t num_interior_0D_cells;
            std::size_t num_exterior_0D_cells;
            std::size_t num_interior_1D_cells;
            std::size_t num_exterior_1D_cells;
            std::size_t num_interior_2D_cells;
            std::size_t num_exterior_2D_cells;

            std::size_t total_cells;
            std::size_t nx; //number of cells in x direction
            std::size_t ny; //number of cells in y direction 
            
            std::size_t ppc; //number of discrete points to form the subcells

            double dx; //width of cell in the x direction 
            double dy; //width of cell in the y direction 
            
            CplexGraph2D_t graph;
            CartesianGrid2D fine_grid;
            CartesianGrid2D coarse_grid;
    };

    } // end namespace DGGML

#endif
