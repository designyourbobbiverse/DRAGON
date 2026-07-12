//
//  Splitting.cpp
//  DRAGON/Hydro/Godunov
//
//  Created by Bobbie Markwick on 16/05/2026.
//  Implementation based mostly on Toro (2009). https://doi.org/10.1007/b79761
//

#include "Grid.hpp"
#include "Boundary.hpp"
#include "CFL.hpp"
#include "DragonWing.hpp"
#include "Riemann.hpp"
#include "TVD.hpp"
#include <math.h>
#include <cassert>


static int validGhosts(int g){
#if defined(MHD) && defined(CTU)
    return std::max(g, 3);
#elif defined(MUSCL_Hancock)
    return std::max(g, 2);
#else
    return std::max(g, 1);
#endif
}
//MARK: Array Wrappers
Grid1D::Grid1D(int s_, double dx_, int g_): w(s_, validGhosts(g_)), dx(dx_) { }
PrimitiveState& Grid1D::operator[](int k) { return w[k]; }
const PrimitiveState& Grid1D::operator[](int k) const { return w[k]; }
int Grid1D::getSize() const { return w.getSize(); }
int Grid1D::getGhosts() const { return w.getGhosts(); }

Grid2D::Grid2D(int nx_, int ny_, double dx_, double dy_, int g_):  w(nx_, ny_,validGhosts(g_)),
#ifdef MHD
    A(nx_+1, ny_+1,w.getGhosts()),
#endif
    dx(dx_), dy(dy_) { }
PrimitiveState& Grid2D::operator[](int i, int j) { return w[i,j]; }
const PrimitiveState& Grid2D::operator[](int i, int j) const { return w[i,j]; }
int Grid2D::getSizeX() const { return w.getSizeX(); }
int Grid2D::getSizeY() const { return w.getSizeY(); }
int Grid2D::getGhosts() const { return w.getGhosts(); }


Grid3D::Grid3D(int nx_, int ny_, int nz_, double dx_, double dy_, double dz_, int g_): w(nx_, ny_, nz_, validGhosts(g_)) ,
#ifdef MHD
    A(nx_+1, ny_+1, nz_+1, w.getGhosts()),
#endif
    dx(dx_), dy(dy_), dz(dz_) {    }
PrimitiveState& Grid3D::operator[](int i, int j, int k) { return w[i,j,k]; }
const PrimitiveState& Grid3D::operator[](int i, int j, int k) const { return w[i,j,k]; }
int Grid3D::getSizeX() const { return w.getSizeX(); }
int Grid3D::getSizeY() const { return w.getSizeY(); }
int Grid3D::getSizeZ() const { return w.getSizeZ(); }
int Grid3D::getGhosts() const { return w.getGhosts(); }



//MARK: Godunov Sweep
namespace Godunov{
void sweep(ExtendedArray1D<PrimitiveState>& w, double dt_dx){
    const int size = w.getSize(), ghosts = w.getGhosts();
    //Compute Fluxes
    ConservativeState fL, fR;
    PrimitiveState _LR = w[-ghosts], _RL, _RR;
    
    //Compute the leftmost flux
    TVD::MUSCL(w[-ghosts], _RL, w[-ghosts+1], _RR, w[-ghosts+2], dt_dx);
    fL = Riemann(_LR, _RL).flux_X(dt_dx);
    
    for(int i=-ghosts+1; i<size+ghosts-1; i++) {
        _LR = _RR;//Move Right state from previous cycle
        //Reconstruct Half-States (if applicable)
        if(i+2 < size+ghosts){
            TVD::MUSCL(w[i], _RL, w[i+1], _RR, w[i+2], dt_dx);
        } else {
            _RL = w[i+1];
        }
        //Compute & Apply Flux
        fR = Riemann(_LR, _RL).flux_X(dt_dx);
        w[i] += (fL - fR) * (dt_dx); //Apply flux to cell
        fL = fR; //Right flux on this cell must equal Left flux on next cell
    }
    //Compute & Apply the rightmost flux
    _LR = _RR;
    _RL = w[size+ghosts-1];
    fR = Riemann(_LR, _RL).flux_X(dt_dx);
    w[size+ghosts-2] += (fL - fR) * (dt_dx); //Apply flux to cell
}
}

//MARK: 1D Advance
void Grid1D::split_step(double dt){ Grid1D::unsplit_step(dt); }
void Grid1D::unsplit_step(double dt){
    boundary.apply(*this);
    Godunov::sweep(w, dt/dx);
}

//MARK: 2D Split
void Grid2D::split_step(double t1){
    //Advance (Strang Split), alternating which step comes first
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
//MARK: 3D Split
void Grid3D::split_step(double t1){
    //Advance  (Strang Split), rotating step orders
    switch (sweep_step++ % 6) {
    case 0: //Cyclic XYZ
        advanceX(t1/2);
        advanceY(t1/2);
        advanceZ(t1);
        advanceY(t1/2);
        advanceX(t1/2);
        break;
    case 1: //Cyclic ZXY
        advanceZ(t1/2);
        advanceX(t1/2);
        advanceY(t1);
        advanceX(t1/2);
        advanceZ(t1/2);
        break;
    case 2: //Cyclic YZX
        advanceY(t1/2);
        advanceZ(t1/2);
        advanceX(t1);
        advanceZ(t1/2);
        advanceY(t1/2);
        break;
    case 3: //Anticyclic ZYX
        advanceZ(t1/2);
        advanceY(t1/2);
        advanceX(t1);
        advanceY(t1/2);
        advanceZ(t1/2);
        break;
    case 4: //Anticyclic XZY
        advanceX(t1/2);
        advanceZ(t1/2);
        advanceY(t1);
        advanceZ(t1/2);
        advanceX(t1/2);
        break;
    case 5: //Anticyclic YXZ
        advanceY(t1/2);
        advanceX(t1/2);
        advanceZ(t1);
        advanceX(t1/2);
        advanceY(t1/2);
        break;
    }
}


//MARK: 2D Component Sweeps
void Grid2D::advanceX(double dt){
    boundary.apply(*this); //Apply Boundary Conditions before every sweep
    const int nx = w.getSizeX(), ny = w.getSizeY(), ghosts = w.getGhosts();
        auto __B = DRAGONWING::requestPrimitiveArrays(1, nx, ghosts);
    ExtendedArray1D<PrimitiveState>& _w = *__B[0];

    for(int j=-ghosts; j<ny+ghosts; j++){
        for(int i=-ghosts; i<nx+ghosts; i++) _w[i] = w[i,j]; //Copy to a 1D array

        Godunov::sweep(_w, dt/dx);
        
        for(int i=-ghosts; i<nx+ghosts; i++) w[i,j] = _w[i]; //Copy 1D array back to grid
    }
}
void Grid2D::advanceY(double dt){
    boundary.apply(*this); //Apply Boundary Conditions before every sweep
    const int nx = w.getSizeX(), ny = w.getSizeY(), ghosts = w.getGhosts();
        auto __B = DRAGONWING::requestPrimitiveArrays(1, ny, ghosts);
    ExtendedArray1D<PrimitiveState>& _w = *__B[0];

    for(int i=-ghosts; i<nx+ghosts; i++){
        for(int j=-ghosts; j<ny+ghosts; j++)  _w[j] = w[i,j].swappedXY(); //Dimension swap + copy to a 1D array
        
        Godunov::sweep(_w, dt/dy);
        
        for(int j=-ghosts; j<ny+ghosts; j++)  w[i,j] = _w[j].swappedXY(); //Dimension swap back + copy back to grid
    }
}


//MARK: 3D Component Sweeps
void Grid3D::advanceX(double dt){
    boundary.apply(*this); //Apply Boundary Conditions before every sweep
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
        auto __B = DRAGONWING::requestPrimitiveArrays(1, nx, ghosts);
    ExtendedArray1D<PrimitiveState>& _w = *__B[0];
    
    for(int k=-ghosts; k<nz+ghosts; k++){
        for(int j=-ghosts; j<ny+ghosts; j++){
            for(int i=-ghosts; i<nx+ghosts; i++) _w[i] = w[i,j,k]; //Copy to a 1D array
            
            Godunov::sweep(_w, dt/dx); //Sweep through the 1D array
            
            for(int i=-ghosts; i<nx+ghosts; i++) w[i,j,k] = _w[i]; //Copy back to grid
        }
    }
}

void Grid3D::advanceY(double dt){
    boundary.apply(*this); //Apply Boundary Conditions before every sweep
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
        auto __B = DRAGONWING::requestPrimitiveArrays(1, ny, ghosts);
    ExtendedArray1D<PrimitiveState>& _w = *__B[0];

    for(int k=-ghosts; k<nz+ghosts; k++){
        for(int i=-ghosts; i<nx+ghosts; i++) {
            for(int j=-ghosts; j<ny+ghosts; j++) _w[j] = w[i,j,k].swappedXY(); //Dimension swap + copy to a 1D array

            Godunov::sweep(_w, dt/dy); //Sweep through the 1D array

            for(int j=-ghosts; j<ny+ghosts; j++) w[i,j,k] = _w[j].swappedXY();  //Dimension swap back + copy back to grid
       }
    }
}
void Grid3D::advanceZ(double dt){
    boundary.apply(*this); //Apply Boundary Conditions before every sweep
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
        auto __B = DRAGONWING::requestPrimitiveArrays(1, nz, ghosts);
    ExtendedArray1D<PrimitiveState>& _w = *__B[0];

    
    for(int i=-ghosts; i<nx+ghosts; i++) {
        for(int j=-ghosts; j<ny+ghosts; j++) {
            for(int k=-ghosts; k<nz+ghosts; k++) _w[k] = w[i,j,k].swappedXZ(); //Dimension swap + copy to a 1D array
            
            Godunov::sweep(_w, dt/dz); //Sweep through the 1D array

            for(int k=-ghosts; k<nz+ghosts; k++) w[i,j,k] = _w[k].swappedXZ(); //Dimension swap back + copy back to grid
       }
    }
}
