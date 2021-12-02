#ifndef __CAJETE_CARTESIANGRID2D_HPP
#define __CAJETE_CARTESIANGRID2D_HPP

#include <cmath>
#include <iostream>

namespace Cajete 
{

class CartesianGrid2D {
    public:
        double _min_x;
        double _min_y;
        double _max_x;
        double _max_y;
        double _dx;
        double _dy;
        double _rdx;
        double _rdy;
        int _nx;
        int _ny;

        CartesianGrid2D () {}
       
        void init(const double min_x, const double min_y,
                  const double max_x, const double max_y,
                  const std::size_t nx, const std::size_t ny) {
            _min_x = min_x;
            _min_y = min_y;
            _max_x = max_x;
            _max_y = max_y;
            
            _nx = nx; 
            _ny = ny;

            _dx = (max_x - min_x) / _nx;
            _dy = (max_y - min_y) / _ny;

            _rdx = 1.0 / _dx;
            _rdy = 1.0 / _dy;
        }
        
        int totalNumCells() const {return _nx * _ny;}

        int cellsBetween(const double max, const double min, const double rdelta) const {
            return floor((max-min)*rdelta);
        }

        void locatePoint( const double xp, const double yp, int& ic, int& jc ) const
        {
            // Since we use a floor function a point on the outer boundary
            // will be found in the next cell, causing an out of bounds error
            ic = cellsBetween( xp, _min_x, _rdx );
            ic = ( ic == _nx ) ? ic - 1 : ic;
            jc = cellsBetween( yp, _min_y, _rdy );
            jc = ( jc == _ny ) ? jc - 1 : jc;
        }
        
        
        void ijCellIndex(const int cardinal, int& ic, int& jc) const 
        {
           ic = cardinal % _nx;
           jc = ( cardinal - (cardinal % _nx) ) / _nx;
        }

        int cardinalCellIndex(const int i, const int j) const
        {
            return (j*_nx) + i;
        }

        void cardinalToPoint(double& xp, double& yp, const int cardinal) const 
        {
            int i, j;
            ijCellIndex(cardinal, i, j);

            xp = i*_dx;
            yp = j*_dy;
        }

        friend std::ostream& operator<<(std::ostream& os, CartesianGrid2D& grid)
        {
            os << "Grid: {\n"
                << "\t{ minx: " << grid._min_x << " , " 
                << "min_y: " << grid._min_y << " , "
                << "max_x: " << grid._max_x << " , "
                << "max_y: " << grid._max_y << " }\n"
                << "\t{ dx: " << grid._dx  << ", dy: " << grid._dy << " }\n"
                << "\t{ nx: " << grid._nx  << ", ny: " << grid._ny << " }\n"
                << "\t{ total number of cells: " << grid.totalNumCells() << " }\n"
                << "}\n";
                return os;
        }
};

} // end namespace Cajete

#endif
