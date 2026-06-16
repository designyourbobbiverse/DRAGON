//
//  Godunov/Splitting.hpp
//  DRAGON
//
//  Created by Bobbie Markwick on 16/05/2026.
//  Implementation based mostly on Toro (2009). https://doi.org/10.1007/b79761
//

#include "Grid.hpp"
#include "Boundary.hpp"
#include "CFL.hpp"
#include <math.h>
#include <cassert>


static int validGhosts(int g){
#ifdef MUSCL_Hancock
    return std::max(g, 2);
#else
    return std::max(g, 1);
#endif
}
//MARK: Array Wrappers
Grid2D::Grid2D(int nx_, int ny_, double dx_, double dy_, int g_): dx(dx_), dy(dy_), w(nx_, ny_,validGhosts(g_)) {}
PrimitiveState& Grid2D::operator[](int i, int j) { return w[i,j]; }
const PrimitiveState& Grid2D::operator[](int i, int j) const { return w[i,j]; }
int Grid2D::getSizeX() const { return w.getSizeX(); }
int Grid2D::getSizeY() const { return w.getSizeY(); }
int Grid2D::getGhosts() const { return w.getGhosts(); }
Grid2D::~Grid2D(){ w.~ExtendedArray2D(); }


Grid3D::Grid3D(int nx_, int ny_, int nz_, double dx_, double dy_, double dz_, int g_): dx(dx_), dy(dy_), dz(dz_),
    w(nx_, ny_, nz_, validGhosts(g_)) {}
PrimitiveState& Grid3D::operator[](int i, int j, int k) { return w[i,j,k]; }
const PrimitiveState& Grid3D::operator[](int i, int j, int k) const { return w[i,j,k]; }
int Grid3D::getSizeX() const { return w.getSizeX(); }
int Grid3D::getSizeY() const { return w.getSizeY(); }
int Grid3D::getSizeZ() const { return w.getSizeZ(); }
int Grid3D::getGhosts() const { return w.getGhosts(); }
Grid3D::~Grid3D(){ w.~ExtendedArray3D(); }


//MARK: 2D Split

void Grid2D::advance_split(double dt, bool check_cfl){
    while(dt > Timestep_Tolerance){
        //CFL Time Constraint
        double t1 = check_cfl ? std::min(dt,CFL::cfl_time(*this)) : dt;
        dt -= t1;
        

        //Advance, alternating which step comes first
        if (sweep_step++ % 2 == 0) {
            advanceX(t1/2);
            advanceY(t1);
            advanceX(t1/2);
        } else {
            advanceY(t1/2);
            advanceX(t1);
            advanceY(t1/2);
        }
    }
}
//MARK: 3D Split
void Grid3D::advance_split(double dt, bool check_cfl){
    while(dt > Timestep_Tolerance){
        //CFL Time Constraint
        double t1 = check_cfl ? std::min(dt,CFL::cfl_time(*this)) : dt;
        dt -= t1;
        //Advance, rotating step orders
        switch (sweep_step++ % 6) {
        case 0:
            advanceX(t1/2);
            advanceY(t1/2);
            advanceZ(t1);
            advanceY(t1/2);
            advanceX(t1/2);
            break;
        case 1:
            advanceZ(t1/2);
            advanceX(t1/2);
            advanceY(t1);
            advanceX(t1/2);
            advanceZ(t1/2);
            break;
        case 2:
            advanceY(t1/2);
            advanceZ(t1/2);
            advanceX(t1);
            advanceZ(t1/2);
            advanceY(t1/2);
            break;
        case 3:
            advanceZ(t1/2);
            advanceY(t1/2);
            advanceX(t1);
            advanceY(t1/2);
            advanceZ(t1/2);
            break;
        case 4:
            advanceX(t1/2);
            advanceZ(t1/2);
            advanceY(t1);
            advanceZ(t1/2);
            advanceX(t1/2);
            break;
        case 5:
            advanceY(t1/2);
            advanceX(t1/2);
            advanceZ(t1);
            advanceX(t1/2);
            advanceY(t1/2);
            break;
        }
    }
}


//MARK: 2D Component Sweeps

void Grid2D::advanceX(double dt){
   boundary.apply(*this);
    
    int nx = w.getSizeX(), ny = w.getSizeY(), ghosts = w.getGhosts();
    Grid1D _w(nx, dx, ghosts);
    ExtendedArray1D<PrimitiveState> _B1(nx, ghosts), _B2(nx, ghosts);
    
    for(int j=-ghosts; j<ny+ghosts; j++){
        for(int i=-ghosts; i<nx+ghosts; i++) _w[i] = w[i,j];

        _w.god_sweep(dt, _B1, _B2);
        
        for(int i=-ghosts; i<nx+ghosts; i++) w[i,j] = _w[i];
    }
}
void Grid2D::advanceY(double dt){
    boundary.apply(*this);

    int nx = w.getSizeX(), ny = w.getSizeY(), ghosts = w.getGhosts();
    Grid1D _w(ny, dy, ghosts);
    ExtendedArray1D<PrimitiveState> _B1(ny, ghosts), _B2(ny, ghosts);
    
    for(int i=-ghosts; i<nx+ghosts; i++){
        for(int j=-ghosts; j<ny+ghosts; j++)  _w[j] = w[i,j].swapXY();

        _w.god_sweep(dt, _B1, _B2);
        
        for(int j=-ghosts; j<ny+ghosts; j++)  w[i,j] = _w[j].swapXY();
    }
}



//MARK: 3D Component Sweeps

void Grid3D::advanceX(double dt){
    boundary.apply(*this);

    int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
    Grid1D _w(nx,dx,ghosts);
    ExtendedArray1D<PrimitiveState> _B1(nx, ghosts), _B2(nx, ghosts);
    
    for(int k=-ghosts; k<nz+ghosts; k++){
        for(int j=-ghosts; j<ny+ghosts; j++){
            for(int i=-ghosts; i<nx+ghosts; i++) _w[i] = w[i,j,k];
            
            _w.god_sweep(dt, _B1, _B2);
            
            for(int i=-ghosts; i<nx+ghosts; i++) w[i,j,k] = _w[i];
        }
    }
}
void Grid3D::advanceY(double dt){
    boundary.apply(*this);

    int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
    Grid1D _w(ny,dy,ghosts);
    ExtendedArray1D<PrimitiveState> _B1(ny, ghosts), _B2(ny, ghosts);
    
    for(int k=-ghosts; k<nz+ghosts; k++){
        for(int i=-ghosts; i<nx+ghosts; i++) {
            for(int j=-ghosts; j<ny+ghosts; j++) _w[j] = w[i,j,k].swapXY();

            _w.god_sweep(dt, _B1, _B2);
            
            for(int j=-ghosts; j<ny+ghosts; j++) w[i,j,k] = _w[j].swapXY();
       }
    }
}
void Grid3D::advanceZ(double dt){
    boundary.apply(*this);

    int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
    Grid1D _w(nz,dz,ghosts);
    ExtendedArray1D<PrimitiveState> _B1(nz, ghosts), _B2(nz, ghosts);
    
    for(int i=-ghosts; i<nx+ghosts; i++) {
        for(int j=-ghosts; j<ny+ghosts; j++) {
            for(int k=-ghosts; k<nz+ghosts; k++) _w[j] = w[i,j,k].swapXZ();
            
            _w.god_sweep(dt, _B1, _B2);
            
            for(int k=-ghosts; k<nz+ghosts; k++) w[i,j,k] = _w[j].swapXZ();
       }
    }
}

