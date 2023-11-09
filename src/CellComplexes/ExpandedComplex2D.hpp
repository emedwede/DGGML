#ifndef EXPANDED_CELL_COMPLEX_2D_HPP
#define EXPANDED_CELL_COMPLEX_2D_HPP

#include "CartesianComplex2D.hpp"

namespace DGGML
{
    template <typename GraphType = CplexGraph2D_t>
    class ExpandedComplex2D : public CartesianComplex2D<CplexGraph2D_t>
    {
        public:
            
            using types = CartesianComplex2D<GraphType>;
            
            CartesianGrid2D reaction_grid;
            std::vector<int> dim_label;
            std::vector<std::size_t> cell_label;
           ExpandedComplex2D() = default; 

           ExpandedComplex2D(std::size_t n, std::size_t m, double d_x, double d_y, bool g_c = false, double _epsilon = 0.0) : CartesianComplex2D(n, m, d_x, d_y, g_c)
           {
                init(n, m, d_x, d_y, g_c, _epsilon);                 
           }     
            
           void init(std::size_t n, std::size_t m, double d_x, double d_y, bool g_c = false, double _epsilon = 0.0)
           {
                CartesianComplex2D<GraphType>::init(n, m, d_x, d_y, g_c);
                if(_epsilon == 0.0)
                    epsilon = 0.03*dx; //currently default epsilon, make it a parameter related to reaction distance later
                else
                    epsilon = _epsilon;
                
                if(g_c)
                {
                    n = n+2; //pad around boundary in dim  
                    m = m+2; //pad around boundary in dim
                } 
                std::size_t n_r_c = std::floor( ( dx ) / epsilon );
                std::size_t m_r_c = std::floor( ( dy ) / epsilon );

                auto approx_epsilon_x = ( dx ) / n_r_c;
                auto approx_epsilon_y = ( dy ) / m_r_c;
                auto n_r = n_r_c * n;
                auto m_r = m_r_c * m;
                //std::cout << "approx_epsilon_x: " << approx_epsilon_x << ", approx_epsilon_y: " << approx_epsilon_y << ", epsilon: " << epsilon << "\n";
                reaction_grid.init(0.0, 0.0, n*d_x, m*d_y, n_r, m_r);
                dim_label.resize(reaction_grid.totalNumCells());
                cell_label.resize(reaction_grid.totalNumCells());
                std::fill(dim_label.begin(), dim_label.end(), 0);
                build();
           }
           
           //TODO: is this overload the right fix?
           CplexGraph2D_t& getGraph()
           {
               return graph;
           }
            
           //TODO: build function is not well designed and needs rework 
           //
           //essentially, the 1D cell assignment would reassign 0D cells if it happened after
           void build() 
           {
               //we need to loop over every cell in the complex and expand
               //it based on it's classification and if the complex has been
               //ghosted
               //do all 0D after
                for(auto iter = graph.node_list_begin(); iter != graph.node_list_end(); iter++) 
                {
                    auto key = iter->first;
                    auto& data = iter->second.getData();
                    
                    if( data.type == 0 && data.interior == true) //found an interior 2D cell 
                    {
                        //we only need to deal with expanding the interior 
                        double cx = data.position[0];
                        double cy = data.position[1];
                        label(data.corners[lower_left][0], data.corners[lower_left][1],
                            data.corners[upper_right][0], data.corners[upper_right][1], 0, key);
                    }
                }
               for(auto iter = graph.node_list_begin(); iter != graph.node_list_end(); iter++) 
               {
                    auto key = iter->first;
                    auto& data = iter->second.getData();
                    
                    //TODO: expand interior edges near the exterior correctly
                    //we don't need to expand the highest dimensional cell 
                    if(data.type == 1) //found a 1D cell 
                    {
                        //we need to determine if the cell is a horizontal or vertical edge
                        //or if it is near the exterior. we can ignore the ghosted region 
                        if(data.interior == true) 
                        {
                            //We already know this expanded complex can only have two 0D nbrs connected to 1D
                            double x1, x2, x3, y1, y2, y3; std::size_t c_count = 0;
                            for(auto jter = graph.out_neighbors_begin(key); jter != graph.out_neighbors_end(key); jter++)
                            {    
                                auto& node_data = graph.findNode(*jter)->second.getData();
                                if(node_data.type == 2)
                                {
                                    if(c_count == 0)
                                    {
                                        x1 = node_data.position[0];
                                        y1 = node_data.position[1];
                                        c_count++;    
                                    } else {
                                        x3 = node_data.position[0];
                                        y3 = node_data.position[1];
                                    }
                                }
                            }
                            x2 = data.position[0]; y2 = data.position[1];
                            
                            auto gamma_epsilon = 2.0*epsilon;//TODO: fix scaling 2.0*epsilon;
                            if(x1 == x2 && x2 == x3) //x is the axis of collinearity
                            {
                                data.corners[lower_left][0] = x3 - gamma_epsilon;
                                data.corners[lower_left][1] = y3;
                                data.corners[lower_right][0] = x3 + gamma_epsilon;
                                data.corners[lower_right][1] = y3;
                                data.corners[upper_left][0] = x1 - gamma_epsilon;
                                data.corners[upper_left][1] = y1;
                                data.corners[upper_right][0] = x1 + gamma_epsilon;
                                data.corners[upper_right][1] = y1;
                             
                            }
                            if(y1 == y2 && y2 == y3) //y is the axis of collinearity
                            {
                                data.corners[lower_left][0] = x3;
                                data.corners[lower_left][1] = y3 - gamma_epsilon;
                                data.corners[lower_right][0] = x3;
                                data.corners[lower_right][1] = y3 + gamma_epsilon;
                                data.corners[upper_left][0] = x1;
                                data.corners[upper_left][1] = y1 - gamma_epsilon;
                                data.corners[upper_right][0] = x1;
                                data.corners[upper_right][1] = y1 + gamma_epsilon;
                            }
                            
//                            for(auto i = 0; i < 4; i++)
//                            {
//                                std::cout << i << ": { ";
//                                for(auto j = 0; j < 2; j++)
//                                {
//                                    std::cout << data.corners[i][j] << " ";
//                                }
//                                std::cout << "}\n";
//                            }
                            double min_x, min_y, max_x, max_y;
                            min_x = data.corners[lower_left][0];
                            min_y = data.corners[lower_left][1];
                            max_x = min_x;
                            max_y = min_y;
                        
                            for(auto i = 0; i < 4; i++)
                            {
                                min_x = ( min_x <= data.corners[i][0] ) ? min_x : data.corners[i][0];
                                max_x = ( max_x >= data.corners[i][0] ) ? max_x : data.corners[i][0];
                                min_y = ( min_y <= data.corners[i][1] ) ? min_y : data.corners[i][1];
                                max_y = ( max_y >= data.corners[i][1] ) ? max_y : data.corners[i][1];
                            }
                            //std::cout << "{ min_x: " << min_x << ", min_y: " << min_y << ", max_x: " << max_x << ", max_y: " << max_y << " }\n";
                            //TODO: need to fix, works temporarily to map the cells to the right cell
                            label(min_x+0.1*epsilon, min_y+0.1*epsilon, max_x-0.1*epsilon, max_y-0.1*epsilon, 1, key);
                        }else {
                            
                        }
                    }
                }
                
                //do all 0D after
                for(auto iter = graph.node_list_begin(); iter != graph.node_list_end(); iter++) 
                {
                    auto key = iter->first;
                    auto& data = iter->second.getData();
                    
                    if( data.type == 0 && data.interior == false) //found an exterior 2D ghost cell 
                    {
                        //we only need to deal with expanding the interior 
                        double cx = data.position[0];
                        double cy = data.position[1];
                        label(data.corners[lower_left][0], data.corners[lower_left][1],
                            data.corners[upper_right][0], data.corners[upper_right][1], 3, key);
                    } 
                    if( data.type == 2 && data.interior == true) //found a 0D cell 
                    {
                        //we only need to deal with expanding the interior 
                        double cx = data.position[0];
                        double cy = data.position[1];
                        auto gamma_epsilon = 3.5*epsilon; //TODO: fix scaling here
                        //TODO: when we do 2.0 epsilon, it's off by 1 since I think we fall onto
                        // the boundary of the other cell and it get's labeled thanks to how cardinal cell
                        // is included
                        //TODO: breaks for rectangular grids, because the expansion width should vary per dimension
                        data.corners[lower_left][0] = cx - gamma_epsilon;
                        data.corners[lower_left][1] = cy - gamma_epsilon;
                        data.corners[lower_right][0] = cx + gamma_epsilon;
                        data.corners[lower_right][1] = cy - gamma_epsilon;
                        data.corners[upper_left][0] = cx - gamma_epsilon;
                        data.corners[upper_left][1] = cy + gamma_epsilon;
                        data.corners[upper_right][0] = cx + gamma_epsilon;
                        data.corners[upper_right][1] = cy + gamma_epsilon;
                        label(data.corners[lower_left][0], data.corners[lower_left][1],
                            data.corners[upper_right][0], data.corners[upper_right][1], 2, key);
                    } 
                }
           }
            

           void label(double ll_x, double ll_y, double ur_x, double ur_y, int d, std::size_t key) 
           {
                int min_i, min_j, max_i, max_j;
                reaction_grid.locatePoint(ll_x, ll_y, min_i, min_j);
                reaction_grid.locatePoint(ur_x, ur_y, max_i, max_j);
                //if(min_i > max_i) std::swap(min_i, max_i);
                //if(min_j > min_i) std::swap(max_i, max_j);
                //std::cout << min_i << " " << min_j << " " << max_i << " " << max_j << "\n";
                for(auto i = min_i; i <= max_i; i++)
                {
                    for(auto j = min_j; j <= max_j; j++)
                    {
                        auto cardinal = reaction_grid.cardinalCellIndex(i, j);
                        dim_label[cardinal] = d;
                        cell_label[cardinal] = key;
                    }
                }
                //std::cout << "pass\n";

           }
            friend std::ostream& operator<<(std::ostream& os, ExpandedComplex2D& geoplex)
            {
                os << "\n---Geoplex---\n";
                os << "\nCoarse Grid: \n" << geoplex.coarse_grid;
                os << "\nFine Grid: \n" << geoplex.fine_grid;
                os << geoplex.graph;
                
                return os;
            }
             
        private:
           double epsilon;
        };
} //ending namespace DGGML

#endif
