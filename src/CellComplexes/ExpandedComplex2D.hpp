#ifndef EXPANDED_CELL_COMPLEX_2D_HPP
#define EXPANDED_CELL_COMPLEX_2D_HPP

#include "CartesianComplex2D.hpp"

namespace Cajete
{
    template <typename GraphType = CplexGraph2D_t>
    class ExpandedComplex2D : public CartesianComplex2D<CplexGraph2D_t>
    {
        public:
           ExpandedComplex2D(std::size_t n, std::size_t m, double d_x, double d_y, bool g_c = false) : CartesianComplex2D(n, m, d_x, d_y, g_c)
           {
                epsilon = 0.1*dx; //currently default epsilon, make it a parameter related to reaction distance later
                build();
           }     

           //TODO: is this overload the right fix?
           CplexGraph2D_t& getGraph()
           {
               return graph;
           }

           void build() 
           {
               //we need to loop over every cell in the complex and expand
               //it based on it's classification and if the complex has been
               //ghosted
               for(auto iter = graph.node_list_begin(); iter != graph.node_list_end(); iter++) 
               {
                    auto key = iter->first;
                    auto& data = iter->second.getData();
                    
                    //we don't need to expand the highest dimensional cell 
                    if(data.type == 1) //found a 1D cell 
                    {
                        //we need to determine if the cell is a horizontal or vertical edge
                        //or if it is near the exterior. we can ignore the ghosted region 
                        if(!ghost_cells && data.interior == true) 
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
                            
                            if(x1 == x2 && x2 == x3) //x is the axis of collinearity 
                            {
                                std::cout << "Collinear in x\n";
                                data.corners[lower_left][0] = x3 - epsilon;
                                data.corners[lower_left][1] = y3;
                                data.corners[lower_right][0] = x3 + epsilon;
                                data.corners[lower_right][1] = y3;
                                data.corners[upper_left][0] = x1 - epsilon;
                                data.corners[upper_left][1] = y1;
                                data.corners[upper_right][0] = x1 + epsilon;
                                data.corners[upper_right][1] = y1;
                             
                            }
                            if(y1 == y2 && y2 == y3) //y is the axis of collinearity
                            {
                                std::cout << "Collinear in y\n";
                                data.corners[lower_left][0] = x3;
                                data.corners[lower_left][1] = y3 - epsilon;
                                data.corners[lower_right][0] = x3;
                                data.corners[lower_right][1] = y3 + epsilon;
                                data.corners[upper_left][0] = x1;
                                data.corners[upper_left][1] = y1 - epsilon;
                                data.corners[upper_right][0] = x1;
                                data.corners[upper_right][1] = y1 + epsilon;
                            }
                        }else {
                            
                        }
                    }
                    if( data.type == 2 && data.interior == true) //found a 0D cell 
                    {
                        //we only need to deal with expanding the interior 
                        double cx = data.position[0];
                        double cy = data.position[1];

                        data.corners[lower_left][0] = cx - epsilon;
                        data.corners[lower_left][1] = cy - epsilon;
                        data.corners[lower_right][0] = cx + epsilon;
                        data.corners[lower_right][1] = cy - epsilon;
                        data.corners[upper_left][0] = cx - epsilon;
                        data.corners[upper_left][1] = cy + epsilon;
                        data.corners[upper_right][0] = cx + epsilon;
                        data.corners[upper_right][1] = cy + epsilon;
                    } 
               }
           }
        private:
           double epsilon;
        };
} //ending namespace Cajete

#endif
