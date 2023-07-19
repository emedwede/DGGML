#ifndef __DGGML_BRICKGRID2D_HPP
#define __DGGML_BRICKGRID2D_HPP

#include <cmath>
#include <iostream>

namespace DGGML
{

class BrickGrid2D 
{
    public: 
        double _min_x;
        double _min_y;
        double _max_x;
        double _max_y;
        double _dx;
        double _dy;
        double _rdx;
        double _rdy;
        std::size_t _num_r; //number of rows
        //Even starts at "zero"
        std::size_t _num_ee; //number of elements in an even row
        std::size_t _num_eo; //number of elements in an odd row
       
        BrickGrid2D() {}
        
        //We don't have to set the number of elements in
        //even and odd rows as the number of rows, but for
        //convenience we do
        BrickGrid2D(const double min_x, const double min_y, 
                  const double max_x, const double max_y,
                  const double delta_x, const double delta_y)
            : _min_x(min_x)
            , _min_y(min_y)
            , _max_x(max_x)
            , _max_y(max_y)
        {
            //The grid construction is smart in the sense that it will
            //make the grid cells a bit larger so that we can have our
            //desired grid size
            _num_ee = cartesianCellsBetween(max_x, min_x, 1.0/delta_x);
            _num_eo = _num_ee+1;
            _num_r  = cartesianCellsBetween(max_y, min_y, 1.0/delta_y);
            _dx = (max_x - min_x)/_num_ee;
            _dy = (max_y - min_y)/_num_r;
            _rdx = 1.0/_dx;
            _rdy = 1.0/_dy;
        }

        ~BrickGrid2D() {}

        std::size_t totalNumCells() const {
            //something like this, in cartesian it's nx*ny*nz
            //for brick grid it's number of n rows + number n+1 rows
            if(_num_r%2) { //ODD
                return  ((_num_r-1)/2)*_num_ee+((_num_r-1)/2)*_num_eo+_num_ee; 
            } else{ //EVEN
                return (_num_r/2)*_num_ee+(_num_r/2)*_num_eo; 
            }
        }
       
        //Computes the number of interior 2D cells
        std::size_t compute_num_2D_zones() {
            return totalNumCells();
        }

        //Computes the number of interior edges between cells
        std::size_t compute_num_1D_zones() {
            if(_num_r%2) {
                return ((_num_r-1)/2)*(_num_ee-1)+((_num_r-1)/2)*_num_ee+(_num_ee-1)+(_num_r-1)*(2*_num_ee); 
            } else {
                return (_num_r/2)*(_num_ee-1)+(_num_r/2)*_num_ee+(_num_r-1)*(2*_num_ee);
            }
        }
        
        //Compute the number of interior vertices, i.e. where cell edges meet
        std::size_t compute_num_0D_zones() {
            return (2*_num_ee-1)*(_num_r-1);
        }


        int cardinalCellIndex(const int i, const int j) const {
            if(j%2) { 
                return ((j-1)/2)*_num_ee+((j-1)/2)*_num_eo+_num_ee+i;
            } else { //EVEN
                return (j/2)*_num_ee+(j/2)*_num_eo+i;
            }
        }
        
        void ijCellIndex(const int cardinal,const int& i, const int& j) {
           //TODO implement conversion 
        }

        void locatePoint(const double xp, const double yp, int& ic, int& jc) const {
            //we first need to determine the row then we can locate the point
            jc = floor((yp - _min_y) / _dy);
            if(jc%2) { //ODD
                if(xp < _dx/2) {
                    ic = 0;
                } else if (xp > _max_x - _dx/2) {
                    ic = _num_eo - 1;
                } else {
                    ic = floor((xp - (_min_x + _dx/2))/_dx)+1;
                }
            } else { //EVEN
                ic = floor((xp - _min_x) / _dx);
            }
        }
        
        void minMaxCellCorners(int ic, int jc, 
                double& x_min, double& y_min,
                double& x_max, double& y_max) const  
        {
            if(jc % 2 == 0) {
                x_min = ic*_dx;
                y_min = jc*_dy;
                x_max = (ic+1)*_dx;
                y_max = (jc+1)*_dy;
            } else {
                if(ic == 0) {
                    x_min = ic*_dx;
                    y_min = jc*_dy;
                    x_max = (ic+0.5)*_dx;
                    y_max = (jc+1)*_dy;
                } else if(ic == _num_eo-1) {
                    x_min = (ic-0.5)*_dx;
                    y_min = jc*_dy;
                    x_max = ic*_dx;
                    y_max = (jc+1)*_dy;
                } else {
                    x_min = (ic-0.5)*_dx;
                    y_min = jc*_dy;
                    x_max = (ic+0.5)*_dx;
                    y_max = (jc+1)*_dy;
                }
            }    
        }

        int cartesianCellsBetween(const double max, const double min, const double rdelta) const {
            return floor((max-min)*rdelta);
        }
};
    
std::ostream& operator<<(std::ostream& os, const BrickGrid2D& grid)
{
    os << "Staggered 2D Grid: {\n"
       << "\t {xmin: " << grid._min_x << ", ymin: " << grid._min_y
       << ", xmax: " << grid._max_x << ", ymax: " << grid._max_y << "}\n"
       << "\t {Number of even row cells: " << grid._num_ee << "}\n"
       << "\t {Number of  odd row cells: " << grid._num_eo << "}\n"
       << "\t {Number of rows: " << grid._num_r << "}\n"
       << "}\n";
    return os;
}

} // end namespace DGGML

#endif
