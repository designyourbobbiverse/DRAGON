//
//  Unsplit.cpp
//  DRAGON/Hydro/Godunov
//
//  Created by Bobbie Markwick on 10/06/2026.
//  Implementation based in part on
//      Toro (2009) https://doi.org/10.1007/b79761
//      Gardiner and Stone (2008) https://arxiv.org/abs/0712.2634



#include "Hydro/Grid.hpp"
#include "Godunov.hpp"

#include "Config.h"
#include "MHD/CT.hpp"      //For MHD

#include "DragonWing.hpp" //For memory management & synchronization

#include <stdexcept>      //For error handling
#include <format>         //For error message formatting
using namespace DRAGON;
using namespace Godunov;

//General Procedure
// 1) MUSCL (preliminary interface states)
// 2) CTU (transverse corrections to interface states)
// 3) Final fluxes
// 4) Update hydro states
// 5) Constrained Transport (E fields then B fields)
// 6) Check that the solution is acceptable (if not, throw exception to trigger step restart)
// 7) If domain-decomposed, wait until everyone finishes (and is error-free)
// 8) Commit the update


//MARK: 2D Unsplit Step

void Grid2D::unsplit_step(double dt){
    const int nx = w.getSizeX(), ny = w.getSizeY(), ghosts = w.getGhosts();
    const double dt_dx = dt/dx, dt_dy = dt/dy;
    
    if (!DRAGONWING::waitForRelease()) return; //On memory-constrained systems, we might have to wait until it's our turn
    
    //Compute Interface States
        auto __half_states = DRAGONWING::requestPrimitiveArrays(4, nx, ny, ghosts);
    FluidArray2D& _xL = *__half_states[0];
    FluidArray2D& _xR = *__half_states[1];
    FluidArray2D& _yL = *__half_states[2];
    FluidArray2D& _yR = *__half_states[3];
    computeHalfStates_X(_xL, (*this), _xR, dt);
    computeHalfStates_Y(_yL, (*this), _yR, dt);
    
    #ifdef CTU //Gardiner and Stone (2008) https://arxiv.org/abs/0712.2634
    #ifdef MHD
        auto __E_half = DRAGONWING::requestVec3Arrays(1, nx+1, ny+1, ghosts);
    MagneticArray2D& _E_half = *__E_half[0];
    ctu_sweep_MHD(_xL, _xR, _yL, _yR, B, w, _E_half, dt_dx, dt_dy);
    #else
    ctu_sweep_hydro(_xL, _xR, _yL, _yR, dt_dx, dt_dy);
    #endif
    #endif

    //Compute Fluxes
        auto __fluxes = DRAGONWING::requestFluxArrays(2, nx, ny, ghosts);
    FluxArray2D& F_X = *__fluxes[0];
    FluxArray2D& F_Y = *__fluxes[1];
    #ifdef MHD //MHD also needs the transverse fluxes in the first ghost layer to calculate E
    computeFlux_X(_xL, _xR, F_X, 0, nx, -1, ny+1, dt_dx);
    computeFlux_Y(_yL, _yR, F_Y, -1, nx+1, 0, ny, dt_dy);
    #else //Hydro doesn't need transverse fluxes in the first ghost layer
    computeFlux_X(_xL, _xR, F_X, 0, nx, 0, ny, dt_dx);
    computeFlux_Y(_yL, _yR, F_Y, 0, nx, 0, ny, dt_dy);
    #endif
        __half_states.release();
    
    //Preliminarily apply all fluxes
        auto __w = DRAGONWING::requestPrimitiveArrays(1,nx, ny, ghosts); //Auto-releases when the function terminates
    FluidArray2D& _w = *__w[0];
    applyFluxes(w, _w, F_X, F_Y, dt_dx, dt_dy);

    //Preliminary CT Update
    #ifdef MHD
        auto __Elec = DRAGONWING::requestVec3Arrays(1, nx+1, ny+1, ghosts);
    //Compute Electric Fields
    MagneticArray2D& E = *__Elec[0];
    CT::computeElectric(E, F_X, F_Y);
    #ifdef CTU
    CT::upwindElectric(E, F_X, F_Y, _E_half);
        __E_half.release();
    #endif
    //Update B
        auto __B = DRAGONWING::requestVec3Arrays(1, nx+1, ny+1, ghosts); //Auto-releases when the function terminates
    MagneticArray2D& _B = *__B[0]; //Done with _E_half, reuse it
    _B.clone(B);
    CT::Faraday(E, _B, dt_dx, dt_dy);
        __Elec.release();
    CT::computeBodyFields(_B, _w);
    #endif

    //Verify Physicality of solution
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            if (!_w[i,j].isPhysical())
                throw std::runtime_error(std::format("Unphysical state would be produced at ({},{})",i,j));
        }
    }
    
    //Passive Scalar Advection
    auto _q = q.advected(F_X, F_Y, w, _w, dt_dx, dt_dy);
        __fluxes.release();
        
    //Wait for any parallel grids to finish
    DRAGONWING::reportCheckpoint1();
    if (!DRAGONWING::waitForCheckpoint1()) return;
    
    //Commit flux updates
    w.clone(_w, false);
    #ifdef MHD
    B.clone(_B, false);
    #endif
    q.clone(*_q);
    
}

//MARK: 3D Unsplit Step

void Grid3D::unsplit_step(double dt){
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
    const double dt_dx = dt/dx, dt_dy = dt/dy, dt_dz = dt/dz;
    
    if (!DRAGONWING::waitForRelease()) return; //On memory-constrained systems, we might have to wait until it's our turn.

    //Compute Half States
        auto __half_states = DRAGONWING::requestPrimitiveArrays(6, nx, ny, nz, ghosts);
    FluidArray3D& _xL = *__half_states[0];
    FluidArray3D& _xR = *__half_states[1];
    FluidArray3D& _yL = *__half_states[2];
    FluidArray3D& _yR = *__half_states[3];
    FluidArray3D& _zL = *__half_states[4];
    FluidArray3D& _zR = *__half_states[5];
    computeHalfStates_X(_xL, (*this), _xR, dt);
    computeHalfStates_Y(_yL, (*this), _yR, dt);
    computeHalfStates_Z(_zL, (*this), _zR, dt);
    
    #ifdef CTU //Gardiner and Stone (2008) https://arxiv.org/abs/0712.2634
    #ifdef MHD
        auto __E_half = DRAGONWING::requestVec3Arrays(1, nx+1, ny+1, nz+1, ghosts);
    MagneticArray3D& _E_half = *__E_half[0];
    ctu_sweep_MHD(_xL, _xR, _yL, _yR, _zL, _zR, B, w, _E_half, dt_dx, dt_dy, dt_dz);
    #else
    ctu_sweep_hydro(_xL, _xR, _yL, _yR, _zL, _zR, dt_dx, dt_dy, dt_dz);
    #endif
    #endif

    //Compute Fluxes
        auto __fluxes = DRAGONWING::requestFluxArrays(3, nx, ny, nz, ghosts);
    FluxArray3D& F_X = *__fluxes[0];
    FluxArray3D& F_Y = *__fluxes[1];
    FluxArray3D& F_Z = *__fluxes[2];
    #ifdef MHD //MHD also needs the transverse fluxes in the first ghost layer to calculate E
    computeFlux_X(_xL, _xR, F_X, 0, nx, -1, ny+1, -1, nz+1, dt_dx);
    computeFlux_Y(_yL, _yR, F_Y, -1, nx+1, 0, ny, -1, nz+1, dt_dy);
    computeFlux_Z(_zL, _zR, F_Z, -1, nx+1, -1, ny+1, 0, nz, dt_dz);
    #else //Hydro doesn't need transverse fluxes in the first ghost layer
    computeFlux_X(_xL, _xR, F_X, 0, nx, 0, ny, 0, nz, dt_dx);
    computeFlux_Y(_yL, _yR, F_Y, 0, nx, 0, ny, 0, nz, dt_dy);
    computeFlux_Z(_zL, _zR, F_Z, 0, nx, 0, ny, 0, nz, dt_dz);
    #endif
        __half_states.release();
    
    //Preliminarily apply all fluxes
        auto __w = DRAGONWING::requestPrimitiveArrays(1, nx, ny, nz, ghosts); //Auto-releases when the function terminates
    FluidArray3D& _w = *__w[0];
    applyFluxes(w, _w, F_X, F_Y, F_Z, dt_dx, dt_dy, dt_dz);
    
    //Preliminary CT Update
    #ifdef MHD
        auto __Elec = DRAGONWING::requestVec3Arrays(1, nx+1, ny+1, nz+1, ghosts);
    //Compute Electric Fields
    MagneticArray3D& E = *__Elec[0];
    CT::computeElectric(E, F_X, F_Y, F_Z);
    #ifdef CTU
    CT::upwindElectric(E, F_X, F_Y, F_Z, _E_half);
        __E_half.release();
    #endif
    //Update B
        auto __B = DRAGONWING::requestVec3Arrays(1, nx+1, ny+1, nz+1, ghosts); //Auto-releases when the function terminates
    MagneticArray3D& _B = *__B[0]; //Done with _E_half, reuse it
    _B.clone(B);
    CT::Faraday(E, _B, dt_dx, dt_dy, dt_dz);
        __Elec.release();
    CT::computeBodyFields(_B, _w);
    #endif
        __fluxes.release();

    
    //Check Physicality
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            for (int k=0; k<nz; k++) {
                if (!_w[i,j,k].isPhysical())
                    throw std::runtime_error(std::format("Unphysical state would be produced at ({},{},{})",i,j,k));
            }
        }
    }
    
    //Passive Scalar Advection
    auto _q = q.advected(F_X, F_Y, F_Z, w, _w, dt_dx, dt_dy, dt_dz);
        __fluxes.release();
    
    //Wait for any parallel grids to finish
    DRAGONWING::reportCheckpoint1();
    if (!DRAGONWING::waitForCheckpoint1()) return;
    
    //Commit Flux updates
    w.clone(_w, false);
    #ifdef MHD
    B.clone(_B, false);
    #endif
    q.clone(*_q);
}
