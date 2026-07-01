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


struct DistGrid1D: public Grid{
    
    DistGrid1D(int nx, double dx, int ghosts=3, bool root=true);
    DistGrid1D(const DistGrid1D&) = delete; //No copying
    DistGrid1D& operator=(const DistGrid1D&) = delete;

    
    //Grid access
    //Can take inputs <0 or >= n to access ghost cells
    PrimitiveState& operator[](int i);
    const PrimitiveState& operator[](int i) const;
    int getSize() const, getGhosts() const;

    //Advance Forward in time
    void advance(double dt, bool check_cfl = true);
    
private:
    Grid1D data;//The Grid containing all fo the reference data
    std::vector<std::unique_ptr<DistGrid1D>> children;//The subgrids that actually do the processing
    int ncx; //The number of children we have
    double size_x, dx; //Total width, width per cell
    int ghosts;
    //Synchronization with children
    void loadFromChildren();
    void pushToChildren();
};


struct DistGrid2D: public Grid{
    
    DistGrid2D(int nx, int ny, double dx, double dy, int ghosts=3, bool root=true);
    DistGrid2D(const DistGrid2D&) = delete; //No copying
    DistGrid2D& operator=(const DistGrid2D&) = delete;

    
    //Grid access
    //Can take inputs <0 or >= n to access ghost cells
    PrimitiveState& operator[](int i,int j);
    const PrimitiveState& operator[](int i,int j) const;
    int getSizeX() const, getSizeY() const, getGhosts() const;
    #ifdef MHD
    //Access Edge Magnetic potentials. Only Az is used in 2D
    //A[i,j] is on the corner between w[i-1,j-1] and w[i,j]
    ExtendedArray2D<vec3>& _A(){ return data._A();}
    void initialize_B_fields();
    #endif

    //Advance Forward in time
    void advance(double dt, bool check_cfl = true);
    
private:
    Grid2D data;//The Grid containing all fo the reference data
    std::vector<std::unique_ptr<DistGrid2D>> children;//The subgrids that actually do the processing
    int ncx, ncy;//The number of children we have in each dimension
    double size_x, size_y;//Total physical size of the Grid
    double dx, dy;//Physical size per cell
    int ghosts;
    //Synchronization with Children
    void loadFromChildren();
    void pushToChildren();
};


struct DistGrid3D: public Grid{
    
    DistGrid3D(int nx, int ny, int nz, double dx, double dy, double dz, int ghosts=3, bool root=true);
    DistGrid3D(const DistGrid3D&) = delete; //No copying
    DistGrid3D& operator=(const DistGrid3D&) = delete;

    
    //Grid access
    //Can take inputs <0 or >= n to access ghost cells
    PrimitiveState& operator[](int i,int j,int k);
    const PrimitiveState& operator[](int i,int j, int k) const;
    int getSizeX() const, getSizeY() const, getSizeZ() const, getGhosts() const;
    #ifdef MHD
    //Access Edge Magnetic potentials.
    //A[i,j,k] is the corner w[i-1/2,j-1/2,k-1/2] to each of the 3 adjacent corners of w[i,j,k]
    ExtendedArray3D<vec3>& _A(){ return data._A();}
    void initialize_B_fields();
    #endif
    
    //Advance Forward in time
    void advance(double dt, bool check_cfl = true);
    
private:
    Grid3D data;//The Grid containing all fo the reference data
    std::vector<std::unique_ptr<DistGrid3D>> children;//The subgrids that actually do the processing
    int ncx, ncy, ncz;//The number of children we have in each dimension
    double size_x, size_y, size_z; //The total width of the grid
    double dx, dy, dz; //Spacing of an individual cell
    int ghosts;

    //Synchronising with children
    void loadFromChildren();
    void pushToChildren();
};


#endif /* AMRGrid_hpp */
