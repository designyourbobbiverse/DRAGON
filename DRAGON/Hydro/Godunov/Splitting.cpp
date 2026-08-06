//
//  Splitting.cpp
//  DRAGON/Hydro/Godunov
//
//  Created by Bobbie Markwick on 16/05/2026.
//  Implementation based mostly on Toro (2009). https://doi.org/10.1007/b79761
//

#include "Grid.hpp"
#include "Godunov.hpp"

#include "Riemann.hpp"  //For Riemann Solvers
#include "TVD.hpp"      //For MUSCL reconstruction
#include "Boundary.hpp" //For boundary.apply

#include "Config.h"
#include <algorithm>      //For std::max
#include "DragonWing.hpp" //For memory management & synchornization
#include <utility>        //For std::move
#include <stdexcept>      //For exception throwing
#include <format>         //For formatting exception messages


//MARK: Godunov Sweep
void Godunov::sweep(FluidArray1D& w, double dt_dx){
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

//MARK: 1D Advance
void Grid1D::split_step(double dt){ Grid1D::unsplit_step(dt); }
void Grid1D::unsplit_step(double dt){
    //Copy steps to a clone
        auto __w = DRAGONWING::requestPrimitiveArrays(1, w.getSize(), w.getGhosts());
    FluidArray1D& _w = *__w[0];
    _w.clone(w);
    //Compute the updated states
    Godunov::sweep(_w, dt/dx);
    //Check physicality before comitting
    for(int i=0; i<w.getSize(); i++){
        if(!_w[i].isPhysical()) throw std::runtime_error(std::format("Unphysical state would be produced at ({})",i));
    }
    DRAGONWING::reportCheckpoint1();
    if(!DRAGONWING::waitForCheckpoint1()) return;
    //Commit updates
    w.clone(_w);
}

//MARK: 2D Split
void Grid2D::split_step(double dt){
    //Clone in case of failure
    Grid2D _w(w.getSizeX(),getSizeY(), dx, dy, getGhosts());
    _w.w.clone(w);
    _w.boundary = std::move(boundary);
    _w.sweep_step = sweep_step;
    
    try{//Advance (Strang Split), alternating which step comes first
        if (_w.sweep_step++ % 2 == 0) {
            _w.advanceX(dt/2);
            _w.advanceY(dt);
            _w.advanceX(dt/2);
        } else {
            _w.advanceY(dt/2);
            _w.advanceX(dt);
            _w.advanceY(dt/2);
        }
    } catch (...){
        boundary = std::move(_w.boundary);
        throw;
    }
    DRAGONWING::reportCheckpoint1();
    if(!DRAGONWING::waitForCheckpoint1()){
        boundary = std::move(_w.boundary);
        return;
    }
    w.clone(_w.w);
    sweep_step = _w.sweep_step;
    boundary = std::move(_w.boundary);
}
//MARK: 3D Split
void Grid3D::split_step(double dt){
    //Clone in case of failure
    Grid3D _w(w.getSizeX(),getSizeY(), getSizeZ(), dx, dy, dz, getGhosts());
    _w.w.clone(w);
    _w.boundary = std::move(boundary);
    _w.sweep_step = sweep_step;

    try{//Advance (Strang Split), rotating step orders
        switch (_w.sweep_step++ % 6) {
        case 0: //Cyclic XYZ
            _w.advanceX(dt/2);
            _w.advanceY(dt/2);
            _w.advanceZ(dt);
            _w.advanceY(dt/2);
            _w.advanceX(dt/2);
            break;
        case 1: //Cyclic ZXY
            _w.advanceZ(dt/2);
            _w.advanceX(dt/2);
            _w.advanceY(dt);
            _w.advanceX(dt/2);
            _w.advanceZ(dt/2);
            break;
        case 2: //Cyclic YZX
            _w.advanceY(dt/2);
            _w.advanceZ(dt/2);
            _w.advanceX(dt);
            _w.advanceZ(dt/2);
            _w.advanceY(dt/2);
            break;
        case 3: //Anticyclic ZYX
            _w.advanceZ(dt/2);
            _w.advanceY(dt/2);
            _w.advanceX(dt);
            _w.advanceY(dt/2);
            _w.advanceZ(dt/2);
            break;
        case 4: //Anticyclic XZY
            _w.advanceX(dt/2);
            _w.advanceZ(dt/2);
            _w.advanceY(dt);
            _w.advanceZ(dt/2);
            _w.advanceX(dt/2);
            break;
        case 5: //Anticyclic YXZ
            _w.advanceY(dt/2);
            _w.advanceX(dt/2);
            _w.advanceZ(dt);
            _w.advanceX(dt/2);
            _w.advanceY(dt/2);
            break;
        }
    } catch (...){
        boundary = std::move(_w.boundary);
        throw;
    }
    DRAGONWING::reportCheckpoint1();
    if(!DRAGONWING::waitForCheckpoint1()){
        boundary = std::move(_w.boundary);
        return;
    }
    //Commit updates
    w.clone(_w.w);
    sweep_step = _w.sweep_step;
    boundary = std::move(_w.boundary);
}


//MARK: 2D Component Sweeps
void Grid2D::advanceX(double dt){
    boundary.apply(*this); //Apply Boundary Conditions before every sweep
    const int nx = w.getSizeX(), ny = w.getSizeY(), ghosts = w.getGhosts();
        auto __B = DRAGONWING::requestPrimitiveArrays(1, nx, ghosts);
    FluidArray1D& _w = *__B[0];

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
    FluidArray1D& _w = *__B[0];

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
    FluidArray1D& _w = *__B[0];
    
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
    FluidArray1D& _w = *__B[0];

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
    FluidArray1D& _w = *__B[0];

    
    for(int i=-ghosts; i<nx+ghosts; i++) {
        for(int j=-ghosts; j<ny+ghosts; j++) {
            for(int k=-ghosts; k<nz+ghosts; k++) _w[k] = w[i,j,k].swappedXZ(); //Dimension swap + copy to a 1D array
            
            Godunov::sweep(_w, dt/dz); //Sweep through the 1D array

            for(int k=-ghosts; k<nz+ghosts; k++) w[i,j,k] = _w[k].swappedXZ(); //Dimension swap back + copy back to grid
       }
    }
}
