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
#include "ExtendedArray.hpp"
#include "Boundary.hpp"
#include "Config.h"


struct Grid1D{
private:
    ExtendedArray1D<PrimitiveState> w;
public:
    double dx;

    Grid1D(int size, double dx, int ghosts=1);
    ~Grid1D();
    Grid1D(const Grid1D&) = delete; //No cop.ying
    Grid1D& operator=(const Grid1D&) = delete;

    //Grid access
    PrimitiveState& operator[](int k);
    const PrimitiveState& operator[](int k) const;
    int getSize() const, getGhosts() const;
    
    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
    
    //Advance forward in time
    void advance(double dt, bool check_cfl=true);
    void god_sweep(double dt, ExtendedArray1D<PrimitiveState>& _L, ExtendedArray1D<PrimitiveState>& _R);
};

struct Grid2D{
private:
    ExtendedArray2D<PrimitiveState> w;
public:
    Grid2D(int nx, int ny, double dx, double dy, int ghosts=1);
    ~Grid2D();
    Grid2D(const Grid2D&) = delete; //No cop.ying
    Grid2D& operator=(const Grid2D&) = delete;
    
    //Grid access
    PrimitiveState& operator[](int i,int j);
    const PrimitiveState& operator[](int i,int j) const;
    int getSizeX() const, getSizeY() const, getGhosts() const;
    double dx, dy;
    
    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
    
    //Advance Forward in time
    void advance(double dt, bool check_cfl = true);
    void advance_split(double dt, bool check_cfl = true);
    void advance_unsplit(double dt, bool check_cfl = true);

private:
    int sweep_step = 0;

    void advanceX(double dt);
    void advanceY(double dt);
    void advanceXY(double dt);
};

struct Grid3D{
private:
    ExtendedArray3D<PrimitiveState>  w;
public:
    double dx, dy, dz;

    Grid3D(int nx, int ny, int nz, double dx, double dy, double dz, int ghosts=1);
    ~Grid3D();
    Grid3D(const Grid3D&) = delete; //No cop.ying
    Grid3D& operator=(const Grid3D&) = delete;

    //Grid Access
    PrimitiveState& operator[](int i,int j,int k);
    const PrimitiveState& operator[](int i,int j,int k) const;
    int getSizeX() const, getSizeY() const, getSizeZ() const, getGhosts() const;

    //Boundary
    Boundary::BoundaryList boundary = Boundary::BoundaryList();
    
    //Advance Forward in time
    void advance(double dt, bool check_cfl = true);
    void advance_split(double dt, bool check_cfl = true);
    void advance_unsplit(double dt, bool check_cfl = true);
    
    
private:
    int sweep_step = 0;

    void advanceX(double dt);
    void advanceY(double dt);
    void advanceZ(double dt);
    void advanceXYZ(double dt);
};


#endif
