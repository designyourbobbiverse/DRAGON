//
//  AMRGrid.hpp
//  DRAGON
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


struct AMRGrid1D{
    
    AMRGrid1D(int nx, double dx, int ghosts=3);
    AMRGrid1D(const AMRGrid1D&) = delete; //No copying
    AMRGrid1D& operator=(const AMRGrid1D&) = delete;

    
    //Grid access
    //Can take inputs <0 or >= n to access ghost cells
    PrimitiveState& operator[](int i);
    const PrimitiveState& operator[](int i) const;
    int getSize() const, getGhosts() const;

    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
    
    //Advance Forward in time
    void advance(double dt, bool check_cfl = true);
    
private:
    Grid1D data;//The Grid containing all fo the reference data
    std::vector<std::unique_ptr<AMRGrid1D>> children;//The subgrids that actually do the processing
    int ncx; //The number of children we have
    double size_x, dx; //Total width, width per cell
    int ghosts;
    //Synchronization with children
    void loadFromChildren();
    void pushToChildren();
};


struct AMRGrid2D{
    
    AMRGrid2D(int nx, int ny, double dx, double dy, int ghosts=3);
    AMRGrid2D(const AMRGrid2D&) = delete; //No copying
    AMRGrid2D& operator=(const AMRGrid2D&) = delete;

    
    //Grid access
    //Can take inputs <0 or >= n to access ghost cells
    PrimitiveState& operator[](int i,int j);
    const PrimitiveState& operator[](int i,int j) const;
    int getSizeX() const, getSizeY() const, getGhosts() const;

    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
    
    //Advance Forward in time
    void advance(double dt, bool check_cfl = true);
    
private:
    Grid2D data;//The Grid containing all fo the reference data
    std::vector<std::unique_ptr<AMRGrid2D>> children;//The subgrids that actually do the processing
    int ncx, ncy;//The number of children we have in each dimension
    double size_x, size_y;//Total physical size of the Grid
    double dx, dy;//Physical size per cell
    int ghosts;
    //Synchronization with Children
    void loadFromChildren();
    void pushToChildren();
};


struct AMRGrid3D{
    
    AMRGrid3D(int nx, int ny, int nz, double dx, double dy, double dz, int ghosts=3);
    AMRGrid3D(const AMRGrid3D&) = delete; //No copying
    AMRGrid3D& operator=(const AMRGrid3D&) = delete;

    
    //Grid access
    //Can take inputs <0 or >= n to access ghost cells
    PrimitiveState& operator[](int i,int j,int k);
    const PrimitiveState& operator[](int i,int j, int k) const;
    int getSizeX() const, getSizeY() const, getSizeZ() const, getGhosts() const;

    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
    
    //Advance Forward in time
    void advance(double dt, bool check_cfl = true);
    
private:
    Grid3D data;//The Grid containing all fo the reference data
    std::vector<std::unique_ptr<AMRGrid3D>> children;//The subgrids that actually do the processing
    int ncx, ncy, ncz;//The number of children we have in each dimension
    double size_x, size_y, size_z; //The total width of the grid
    double dx, dy, dz; //Spacing of an individual cell
    int ghosts;

    //Synchronising with children
    void loadFromChildren();
    void pushToChildren();
};


#endif /* AMRGrid_hpp */
