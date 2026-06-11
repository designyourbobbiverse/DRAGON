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
#include "Config.h"

struct Grid1D{
    Grid1D(int size, double dx, int ghosts=1);
    ~Grid1D();
    Grid1D(const Grid1D&) = delete; //No copying
    Grid1D& operator=(const Grid1D&) = delete;
    
    //Advance forward in time
    void advance(double dt);
    void god_sweep(double dt, Grid1D& _L, Grid1D& _R);

    //Grid access
    PrimitiveState& operator[](int k);
    const PrimitiveState& operator[](int k) const;
    int getSize() const, getGhosts() const;
    double dx;
    
    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
private:
    PrimitiveState* w;
    int ghosts, size;
};

struct Grid2D{

    Grid2D(int nx, int ny, double dx, double dy, int ghosts=1);
    ~Grid2D();
    Grid2D(const Grid2D&) = delete; //No copying
    Grid2D& operator=(const Grid2D&) = delete;
    
    //Grid access
    PrimitiveState& operator[](int i,int j);
    const PrimitiveState& operator[](int i,int j) const;
    int getSizeX() const, getSizeY() const, getGhosts() const;
    double dx, dy;
    
    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
    
    //Advance Forward in time
    void advance(double dt);
    void advance_split(double dt);
    void advance_unsplit(double dt);

private:
    PrimitiveState* w;
    int ghosts, nx, ny;
    int sweep_step = 0;

    void advanceX(double dt);
    void advanceY(double dt);
    void advanceXY(double dt);
    void computeHalfStates_X(Grid2D& _L, Grid2D& _R, double dt);
    void computeHalfStates_Y(Grid2D& _L, Grid2D& _R, double dt);
};

struct Grid3D{
    Grid3D(int nx, int ny, int nz, double dx, double dy, double dz, int ghosts=1);
    ~Grid3D();
    Grid3D(const Grid3D&) = delete; //No copying
    Grid3D& operator=(const Grid3D&) = delete;

    //Grid Access
    PrimitiveState& operator[](int i,int j,int k);
    const PrimitiveState& operator[](int i,int j,int k) const;
    int getSizeX() const, getSizeY() const, getSizeZ() const, getGhosts() const;
    double dx, dy, dz;

    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
    
    //Advance Forward in time
    void advance(double dt);
    void advance_split(double dt);
    void advance_unsplit(double dt);
    
    
private:
    PrimitiveState* w;
    int ghosts, nx, ny, nz;
    int sweep_step = 0;

    void advanceX(double dt);
    void advanceY(double dt);
    void advanceZ(double dt);
    void advanceXYZ(double dt);
    void computeHalfStates_X(Grid3D& _L, Grid3D& _R, double dt);
    void computeHalfStates_Y(Grid3D& _L, Grid3D& _R, double dt);
    void computeHalfStates_Z(Grid3D& _L, Grid3D& _R, double dt);


};


#endif
