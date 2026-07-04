//
//  DistGrid.hpp
//  DRAGON/Refinement
//  User-Facing Header file
//
//  Created by Bobbie Markwick on 12/06/2026.
//

#ifndef AMRGrid_hpp
#define AMRGrid_hpp

#include "Grid.hpp"
#include "Config.h"
#include <vector>
#include <memory>


struct DistGrid1D: public Grid1D{
    
    DistGrid1D(int nx, double dx, int ghosts=3, bool root=true);
    DistGrid1D(const DistGrid1D&) = delete; //No copying
    DistGrid1D& operator=(const DistGrid1D&) = delete;

    //Advance Forward in time
    void advance(double dt, bool check_cfl = true) override;
    
private:
    std::vector<std::unique_ptr<DistGrid1D>> children;//The subgrids that actually do the processing
    int ncx; //The number of children we have
    double size_x; //Total width
    //Synchronization with children
    void loadFromChildren();
    void pushToChildren();
};


struct DistGrid2D: public Grid2D{
    
    DistGrid2D(int nx, int ny, double dx, double dy, int ghosts=3, bool root=true);
    DistGrid2D(const DistGrid2D&) = delete; //No copying
    DistGrid2D& operator=(const DistGrid2D&) = delete;

    //Advance Forward in time
    void advance(double dt, bool check_cfl = true) override;
    
private:
    std::vector<std::unique_ptr<DistGrid2D>> children;//The subgrids that actually do the processing
    int ncx, ncy;//The number of children we have in each dimension
    double size_x, size_y;//Total physical size of the Grid
    //Synchronization with Children
    void loadFromChildren();
    void pushToChildren();
};


struct DistGrid3D: public Grid3D{
    
    DistGrid3D(int nx, int ny, int nz, double dx, double dy, double dz, int ghosts=3, bool root=true);
    DistGrid3D(const DistGrid3D&) = delete; //No copying
    DistGrid3D& operator=(const DistGrid3D&) = delete;

    
    //Advance Forward in time
    void advance(double dt, bool check_cfl = true) override;
    
private:
    std::vector<std::unique_ptr<DistGrid3D>> children;//The subgrids that actually do the processing
    int ncx, ncy, ncz;//The number of children we have in each dimension
    double size_x, size_y, size_z; //The total width of the grid

    //Synchronising with children
    void loadFromChildren();
    void pushToChildren();
};


#endif /* AMRGrid_hpp */
