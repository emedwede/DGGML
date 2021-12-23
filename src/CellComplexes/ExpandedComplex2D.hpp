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
                build();
           }     

           void build() 
           {
               //we need to loop over every cell in the complex and expand
               //it based on it's classification and if the complex has been
               //ghosted
               for(auto iter = graph.node_list_begin(); iter != graph.node_list_end(); iter++) 
               {
                    auto key = iter->first;
                    auto data = iter->second.getData();

                    if(data.type == 2) //found a 2D cell 
                    {
                        //regardless of ghost cells, we full expand 2D based on the grid

                    } else if (data.type == 1) //found a 1D cell 
                    {
                        //we need to determine if the cell is a horizontal or vertical edge
                        //or if it is near the exterior. we can ignore the ghosted region 
                        if(!ghost_cells) 
                        {
                        
                        } else {
                            
                        }
                    } else if( data.type == 0 ) //found a 0D cell 
                    {
                        //we only need to deal with expanding the interior 
                        
                    } 
               }
           }
        };
} //ending namespace Cajete

#endif
