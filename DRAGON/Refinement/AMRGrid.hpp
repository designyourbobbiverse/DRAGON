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
    ~AMRGrid1D();
    AMRGrid1D(const AMRGrid1D&) = delete; //No cop.ying
    AMRGrid1D& operator=(const AMRGrid1D&) = delete;

    
    //Grid access
    PrimitiveState& operator[](int i);
    const PrimitiveState& operator[](int i) const;
    int getSize() const, getGhosts() const;

    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
    
    //Advance Forward in time
    void advance(double dt, bool check_cfl = true);
    
private:
    Grid1D data;
    std::vector<std::unique_ptr<AMRGrid1D>> children;
    int ghosts, ncx;
    double size_x, dx;
    
    void loadFromChildren();
    void pushToChildren();
};


struct AMRGrid2D{
    
    AMRGrid2D(int nx, int ny, double dx, double dy, int ghosts=3);
    ~AMRGrid2D();
    AMRGrid2D(const AMRGrid2D&) = delete; //No cop.ying
    AMRGrid2D& operator=(const AMRGrid2D&) = delete;

    
    //Grid access
    PrimitiveState& operator[](int i,int j);
    const PrimitiveState& operator[](int i,int j) const;
    int getSizeX() const, getSizeY() const, getGhosts() const;

    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
    
    //Advance Forward in time
    void advance(double dt, bool check_cfl = true);
    
private:
    Grid2D data;
    std::vector<std::unique_ptr<AMRGrid2D>> children;
    int ghosts, ncx, ncy;
    double size_x, size_y, dx, dy;
    
    void loadFromChildren();
    void pushToChildren();
};


struct AMRGrid3D{
    
    AMRGrid3D(int nx, int ny, int nz, double dx, double dy, double dz, int ghosts=3);
    ~AMRGrid3D();
    AMRGrid3D(const AMRGrid3D&) = delete; //No cop.ying
    AMRGrid3D& operator=(const AMRGrid3D&) = delete;

    
    //Grid access
    PrimitiveState& operator[](int i,int j,int k);
    const PrimitiveState& operator[](int i,int j, int k) const;
    int getSizeX() const, getSizeY() const, getSizeZ() const, getGhosts() const;

    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
    
    //Advance Forward in time
    void advance(double dt, bool check_cfl = true);
    
private:
    Grid3D data;
    std::vector<std::unique_ptr<AMRGrid3D>> children;
    int ghosts, ncx, ncy, ncz;
    double size_x, size_y, size_z, dx, dy, dz;
    
    void loadFromChildren();
    void pushToChildren();
};


#endif /* AMRGrid_hpp */
