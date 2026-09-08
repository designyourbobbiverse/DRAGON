//
//  Splitting.cpp
//  DRAGON/Hydro/Godunov
//
//  Created by Bobbie Markwick on 16/05/2026.
//  Implementation based mostly on Toro (2009). https://doi.org/10.1007/b79761
//

#include "Hydro/Grid.hpp"
#include "Godunov.hpp"

#include "Hydro/Riemann/Riemann.hpp"    //For Riemann Solvers
#include "Hydro/Reconstruction/TVD.hpp" //For MUSCL reconstruction
#include "Boundary/Boundary.hpp"        //For boundary.apply

#include "Config.h"
#include <algorithm>      //For std::max
#include "DragonWing.hpp" //For memory management & synchornization
#include <utility>        //For std::move
#include <stdexcept>      //For exception throwing
#include <format>         //For formatting exception messages
using namespace DRAGON;


//MARK: Godunov Sweep
void Godunov::sweep(FluidArray1D& w, double dt_dx, PassiveArray1D& q){
    const int size = w.getSize(), ghosts = w.getGhosts();
    //Buffers for passives
        auto __w = DRAGONWING::requestPrimitiveArrays(1, w.getSize(), w.getGhosts());
        auto __f = DRAGONWING::requestFluxArrays(1, w.getSize(), w.getGhosts());
    FluxArray1D& f = *__f[0];
    FluidArray1D& _w = *__w[0];
    _w.clone(w);
    
    //Compute Fluxes
    PrimitiveState _LR = w[-ghosts], _RL, _RR;
    //Compute the leftmost flux
    TVD::MUSCL(w[-ghosts], _RL, w[-ghosts+1], _RR, w[-ghosts+2], dt_dx);
    f[-ghosts+1] = Riemann(_LR, _RL).flux_X(dt_dx);
    //Sweep through the array
    for (int i=-ghosts+2; i<size+ghosts; i++) {
        _LR = _RR;//Move Right state from previous cycle
        //Reconstruct Half-States (if applicable)
        if (i+1 < size+ghosts) {
            TVD::MUSCL(w[i-1], _RL, w[i], _RR, w[i+1], dt_dx);
        } else {
            _RL = w[i];
        }
        //ComputeFlux
        f[i] = Riemann(_LR, _RL).flux_X(dt_dx);
    }

    //Apply Flux
    for (int i=-ghosts+1; i<size+ghosts-1; i++){
        w[i] += (f[i] - f[i+1]) * dt_dx; //Apply flux to cell
    }
    //Advect passives
    auto _q = q.advected(f, _w, w, dt_dx);
    q.clone(*_q);
}

//MARK: 1D Advance
void Grid1D::split_step(double dt){ Grid1D::unsplit_step(dt); }
void Grid1D::unsplit_step(double dt){
    //Do everything on a clone in case we need to restart the step
        auto __w = DRAGONWING::requestPrimitiveArrays(1, w.getSize(), w.getGhosts());
    FluidArray1D& _w = *__w[0];
    _w.clone(w);
    PassiveArray1D _q(w.getSize(), 1);
    _q.clone(q);
    
    //Compute the updated states
    Godunov::sweep(_w, dt/dx, _q);
    //Check physicality before comitting
    for (int i=0; i<w.getSize(); i++) {
        if (!_w[i].isPhysical()) throw std::runtime_error(std::format("Unphysical state would be produced at ({})",i));
    }
    //If in a domain-composed group, wait for the other grids to finish before committing
    DRAGONWING::reportCheckpoint1();
    if (!DRAGONWING::waitForCheckpoint1()) return; //Only proceed once everyone is done and if nobody had an error
    //Commit updates
    w.clone(_w);
    q.clone(_q);
}

//MARK: 2D Split
void Grid2D::split_step(double dt){
    //Clone in case of failure
    Grid2D _w(w.getSizeX(),getSizeY(), dx, dy, getGhosts());
    _w.w.clone(w);
    _w.boundary = std::move(boundary);
    _w.sweep_step = sweep_step;
    
    try{//Advance (Strang Split), alternating which step comes first to reduce directional bias
        if (_w.sweep_step++ % 2 == 0) {
            _w.advanceX(dt/2);
            _w.advanceY(dt);
            _w.advanceX(dt/2);
        } else {
            _w.advanceY(dt/2);
            _w.advanceX(dt);
            _w.advanceY(dt/2);
        }
    } catch (...) {
        boundary = std::move(_w.boundary);
        throw;
    }
    //If in a domain-composed group, wait for the other grids to finish before committing
    DRAGONWING::reportCheckpoint1();
    if (!DRAGONWING::waitForCheckpoint1()) {
        boundary = std::move(_w.boundary);
        return; //Somebody had an error, have to restart
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

    try{//Advance (Strang Split), rotating step orders to reduce directional bias
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
    } catch (...) {
        boundary = std::move(_w.boundary);
        throw;
    }
    //If in a domain-composed group, wait for the other grids to finish before committing
    DRAGONWING::reportCheckpoint1();
    if (!DRAGONWING::waitForCheckpoint1()) {
        boundary = std::move(_w.boundary);
        return; //Somebody had an error, have to restart
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
    PassiveArray1D _q {nx, 1};
    
    for (int j=-ghosts; j<ny+ghosts; j++) {
        for (int i=-ghosts; i<nx+ghosts; i++){ //Copy to a 1D array
            _w[i] = w[i,j];
            if(i >= -1 && i <= nx && j >= -1 && j <= ny) _q[i] = q[i,j];
        }

        Godunov::sweep(_w, dt/dx, _q);
        
        for (int i=-ghosts; i<nx+ghosts; i++) { //Copy 1D array back to grid
            w[i,j] = _w[i];
            if(i >= -1 && i <= nx && j >= -1 && j <= ny) q[i,j] = _q[i];
        }
    }
}
void Grid2D::advanceY(double dt){
    boundary.apply(*this); //Apply Boundary Conditions before every sweep
    const int nx = w.getSizeX(), ny = w.getSizeY(), ghosts = w.getGhosts();
        auto __B = DRAGONWING::requestPrimitiveArrays(1, ny, ghosts);
    FluidArray1D& _w = *__B[0];
    PassiveArray1D _q {ny, 1};

    for (int i=-ghosts; i<nx+ghosts; i++) {
        for (int j=-ghosts; j<ny+ghosts; j++){ //Dimension swap + copy to a 1D array
            _w[j] = w[i,j].swappedXY();
            if(i >= -1 && i <= nx && j >= -1 && j <= ny) _q[j] = q[i,j];
        }

        Godunov::sweep(_w, dt/dx, _q);
        
        for (int j=-ghosts; j<ny+ghosts; j++) { //Dimension swap + copy to a 1D array
            w[i,j] = _w[j].swappedXY();
            if(i >= -1 && i <= nx && j >= -1 && j <= ny) q[i,j] = _q[j];
        }
    }
}


//MARK: 3D Component Sweeps
void Grid3D::advanceX(double dt){
    boundary.apply(*this); //Apply Boundary Conditions before every sweep
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
        auto __B = DRAGONWING::requestPrimitiveArrays(1, nx, ghosts);
    FluidArray1D& _w = *__B[0];
    PassiveArray1D _q {nx, 1};
    
    for (int k=-ghosts; k<nz+ghosts; k++) {
        for (int j=-ghosts; j<ny+ghosts; j++) {
            for (int i=-ghosts; i<nx+ghosts; i++){ //Copy to a 1D array
                _w[i] = w[i,j,k];
                if(i >= -1 && i <= nx && j >= -1 && j <= ny & k >= -1 && k <= nz) _q[i] = q[i,j,k];
            }
            
            Godunov::sweep(_w, dt/dx, _q); //Sweep through the 1D array
            
            for (int i=-ghosts; i<nx+ghosts; i++) { //Copy 1D array back to grid
                w[i,j,k] = _w[i];
                if(i >= -1 && i <= nx && j >= -1 && j <= ny & k >= -1 && k <= nz) q[i,j,k] = _q[i];
            }
        }
    }
}

void Grid3D::advanceY(double dt){
    boundary.apply(*this); //Apply Boundary Conditions before every sweep
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
        auto __B = DRAGONWING::requestPrimitiveArrays(1, ny, ghosts);
    FluidArray1D& _w = *__B[0];
    PassiveArray1D _q {ny, 1};

    for (int k=-ghosts; k<nz+ghosts; k++) {
        for (int i=-ghosts; i<nx+ghosts; i++) {
            for (int j=-ghosts; j<ny+ghosts; j++){ //Dimension swap + copy to a 1D array
                _w[j] = w[i,j,k].swappedXY();
                if(i >= -1 && i <= nx && j >= -1 && j <= ny & k >= -1 && k <= nz) _q[j] = q[i,j,k];
            }
            Godunov::sweep(_w, dt/dy, _q); //Sweep through the 1D array

            for (int j=-ghosts; j<ny+ghosts; j++) { //Dimension swap + copy to a 1D array
                w[i,j,k] = _w[j].swappedXY();
                if(i >= -1 && i <= nx && j >= -1 && j <= ny & k >= -1 && k <= nz) q[i,j,k] = _q[j];
            }
       }
    }
}
void Grid3D::advanceZ(double dt){
    boundary.apply(*this); //Apply Boundary Conditions before every sweep
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
        auto __B = DRAGONWING::requestPrimitiveArrays(1, nz, ghosts);
    FluidArray1D& _w = *__B[0];
    PassiveArray1D _q {nz, 1};
    
    for (int i=-ghosts; i<nx+ghosts; i++) {
        for (int j=-ghosts; j<ny+ghosts; j++) {
            for (int k=-ghosts; k<nz+ghosts; k++) { //Dimension swap + copy to a 1D array
                _w[k] = w[i,j,k].swappedXZ();
                if(i >= -1 && i <= nx && j >= -1 && j <= ny & k >= -1 && k <= nz) _q[k] = q[i,j,k];
            }
            
            Godunov::sweep(_w, dt/dz, _q); //Sweep through the 1D array

            for (int k=-ghosts; k<nz+ghosts; k++) { //Dimension swap + copy to a 1D array
                w[i,j,k] = _w[k].swappedXZ();
                if(i >= -1 && i <= nx && j >= -1 && j <= ny & k >= -1 && k <= nz) q[i,j,k] = _q[k];
            }
       }
    }
}
