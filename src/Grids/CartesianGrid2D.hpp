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
        int _nx; //number of cells in the x direction
        int _ny; //number of cells in the y direction
        int _px; //number of grid points in the x direction 
        int _py; //number of grid points in the y direction

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
            _px = nx+1;
            _py = ny+1;

            _dx = (max_x - min_x) / _nx;
            _dy = (max_y - min_y) / _ny;

            _rdx = 1.0 / _dx;
            _rdy = 1.0 / _dy;
        }
        
        //total number of cells formed from the lattice
        int totalNumCells() const {return _nx * _ny;}
       
        //total number of lattice grid points
        int totalNumPoints() const {return _px * _py;}

        int cellsBetween(const double max, const double min, const double rdelta) const {
            return floor((max-min)*rdelta);
        }
        
        //locate which cell the poin belongs to
        void locatePoint( const double xp, const double yp, int& ic, int& jc ) const
        {
            // Since we use a floor function a point on the outer boundary
            // will be found in the next cell, causing an out of bounds error
            ic = cellsBetween( xp, _min_x, _rdx );
            ic = ( ic == _nx ) ? ic - 1 : ic;
            jc = cellsBetween( yp, _min_y, _rdy );
            jc = ( jc == _ny ) ? jc - 1 : jc;
        }
        
        //returns the ij index of a cell from a cardinal index 
        void ijCellIndex(const int cardinal, int& ic, int& jc) const 
        {
           ic = cardinal % _nx;
           jc = ( cardinal - (cardinal % _nx) ) / _nx;
        }
        
        //returns the ij index of lattice point from a caridinal index
        void ijLatticeIndex(const int cardinal, int& ic, int& jc) const 
        {
            ic = cardinal % _px;
            jc = (cardinal - (cardinal % _px)) / _px;
        }

        //returns the cardinal cell index from ij 
        int cardinalCellIndex(const int i, const int j) const
        {
            return (j*_nx) + i;
        }
        
        //returns the cardianl lattice point index from ij 
        int cardinalLatticeIndex(const int i, const int j) const 
        {
            return (j*_px) + i;
        }

        //converts the cardinal index of a cell to it's central point
        void cardinalCellToPoint(double& xp, double& yp, const int cardinal) const 
        {
            int i, j;
            ijCellIndex(cardinal, i, j);

            xp = i*_dx + _dx/2.0;
            yp = j*_dy + _dy/2.0;
        }
        
        //converts the cardinal index of a lattice point to it's exact point 
        void cardinalLatticeToPoint(double& xp, double& yp, const int cardinal) const 
        {
            int i, j;
            ijLatticeIndex(cardinal, i, j);
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
                << "\t{ total number of lattice points: " << grid.totalNumPoints() << " }\n"
                << "}\n";
                return os;
        }
};

} // end namespace Cajete

#endif
