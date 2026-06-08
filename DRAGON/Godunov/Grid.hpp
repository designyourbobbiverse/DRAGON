//
//  Godunov.hpp
//  DRAGON
//  User-Facing Header file
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#ifndef Godunov_hpp
#define Godunov_hpp

#include "FluidElement.hpp"
#include "Boundary.hpp"


struct Grid1D{
    Grid1D(int size, double dx, int ghosts=0);
    ~Grid1D();
    Grid1D(const Grid1D&) = delete; //No copying
    Grid1D& operator=(const Grid1D&) = delete;
    
    //Advance forward in time
    void advance(double dt);
    void god_sweep(double dt, Grid1D& _L, Grid1D& _R);

    //Grid access
    PrimitiveState& operator[](int k);
    const PrimitiveState& operator[](int k) const;
    int getSize(), getGhosts();
    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
private:
    PrimitiveState* w;
    int ghosts, size;
    double dx;
};

struct Grid2D{

    Grid2D(int nx, int ny, double dx, double dy, int ghosts=0);
    ~Grid2D();
    Grid2D(const Grid2D&) = delete; //No copying
    Grid2D& operator=(const Grid2D&) = delete;
    
    //Advance Forward in time
    void advance(double dt);
    
    //Grid access
    PrimitiveState& operator[](int i,int j);
    const PrimitiveState& operator[](int i,int j) const;
    int getSizeX(), getSizeY(), getGhosts();
    
    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();

private:
    PrimitiveState* w;
    int ghosts, nx, ny;
    double dx, dy;

    void advanceX(double dt);
    void advanceY(double dt);
    
};

struct Grid3D{
    Grid3D(int nx, int ny, int nz, double dx, double dy, double dz, int ghosts=0);
    ~Grid3D();
    Grid3D(const Grid3D&) = delete; //No copying
    Grid3D& operator=(const Grid3D&) = delete;

    //Advance Forward in time
    void advance(double dt);

    PrimitiveState& operator[](int i,int j,int k);
    const PrimitiveState& operator[](int i,int j,int k) const;
    int getSizeX(), getSizeY(), getSizeZ(), getGhosts();

    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
private:
    PrimitiveState* w;
    int ghosts, nx, ny, nz;
    double dx, dy, dz;
    
    void advanceX(double dt);
    void advanceY(double dt);
    void advanceZ(double dt);
};


#endif
