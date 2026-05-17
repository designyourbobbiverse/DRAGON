//
//  Godunov.hpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#ifndef Godunov_hpp
#define Godunov_hpp

#include "FluidElement.hpp"

struct Grid1D{
    int size;

    Grid1D(int size, int ghosts=0);
    ~Grid1D();
    Grid1D(const Grid1D&) = delete; //No copying
    Grid1D& operator=(const Grid1D&) = delete;
    
    
    void advance(double dt, double dx);
    void advance(double dt, double dx, Grid1D& _L, Grid1D& _R);
    
    PrimitiveState& operator[](int k);
    const PrimitiveState& operator[](int k) const;

    
private:
    PrimitiveState* w;
    int ghosts;
    
    
};

struct Grid2D{
    int nx,ny;
    double dx, dy;

    Grid2D(int nx, int ny, int ghosts=0);
    ~Grid2D();
    Grid2D(const Grid2D&) = delete; //No copying
    Grid2D& operator=(const Grid2D&) = delete;
    
    PrimitiveState& operator[](int i,int j);
    const PrimitiveState& operator[](int i,int j) const;

    void advance(double dt);

private:
    PrimitiveState* w;
    int ghosts;

    void advanceX(double dt);
    void advanceY(double dt);
    
};

struct Grid3D{
    int nx,ny,nz;
    double dx, dy, dz;
    
    Grid3D(int nx, int ny, int nz, int ghosts=0);
    ~Grid3D();
    Grid3D(const Grid3D&) = delete; //No copying
    Grid3D& operator=(const Grid3D&) = delete;

    PrimitiveState& operator[](int i,int j,int k);
    const PrimitiveState& operator[](int i,int j,int k) const;
    
    void advance(double dt);
    
private:
    PrimitiveState* w;
    int ghosts;
    
    void advanceX(double dt);
    void advanceY(double dt);
    void advanceZ(double dt);
};


#endif /* Godunov_hpp */
