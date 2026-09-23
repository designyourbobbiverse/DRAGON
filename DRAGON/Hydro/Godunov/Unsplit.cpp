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
// 1) MUSCL (preliminary interface states. Includes half a source step)
// 2) Preliminary fluxes
// 3) Compute a half-step predictor value (needed for MHD CTU, non-stiff Source Terms)
// 4) CTU (transverse corrections to preliminary interface states)
// 5) Final fluxes
// 6) Update hydro states (apply fluxes & sources)
// 7) Constrained Transport (E fields then B fields)
// 8) Check that the solution is acceptable (if not, throw exception to trigger step restart)
// 9) If domain-decomposed, wait until everyone finishes (and is error-free)
// 10) advect passive scalars
// 11) Commit the update


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
    
    //Compute Initial Fluxes
        auto __fluxes = DRAGONWING::requestFluxArrays(2, nx, ny, ghosts);
    FluxArray2D& F_X = *__fluxes[0];
    FluxArray2D& F_Y = *__fluxes[1];
    #ifdef MHD //MHD needs an extra layer to calculate E
    computeFlux_X(_xL, _xR, F_X, -1, nx+1, -2, ny+2,  dt_dx);
    computeFlux_Y(_yL, _yR, F_Y, -2, nx+2, -1, ny+1,  dt_dy);
    #else //Hydro doesn't have electric fields
    computeFlux_X(_xL, _xR, F_X, 0, nx, -1, ny+1, dt_dx);
    computeFlux_Y(_yL, _yR, F_Y, -1, nx+1, 0, ny, dt_dy);
    #endif
        #ifndef CTU
        __half_states.release();
        #endif
    
    //Calculate half-time state (needed for source terms + constrained transport)
        auto __w_half = DRAGONWING::requestPrimitiveArrays(1, nx, ny,  ghosts);
    FluidArray2D& _w_half = *__w_half[0];
    applyFluxes(w, _w_half, F_X, F_Y, sources.non_stiff_terms(), w, dt_dx*0.5, dt_dy*0.5, dt*0.5, 1);

    #ifdef MHD //Preliminary CT update
        auto __mag_half = DRAGONWING::requestVec3Arrays(2, nx+1, ny+1, ghosts);
    MagneticArray2D &_B_half = *__mag_half[0], &_E = *__mag_half[1];
    //Compute Electric fields
    CT::computeElectric(_E, F_X, F_Y, 1);
    CT::bodyElectric(w, _B_half, 1); //Use _B_half buffer for Eref since we don't need it yet
    CT::upwindElectric(_E, F_X, F_Y, _B_half,1);
    //Compute Magnetic Fields
    _B_half.clone(B);
    CT::Faraday(_E, _B_half, 0.5*dt_dx, 0.5*dt_dy, 1);
    CT::computeBodyFields(_B_half, _w_half);
    #endif
    
    
    //CTU corrections & Final Fluxes
    #ifdef CTU //Gardiner and Stone (2008) https://arxiv.org/abs/0712.2634
    #ifdef MHD
    ctu_sweep_MHD(_xL, _xR, _yL, _yR, F_X, F_Y, _E, _B_half, w, B, dt_dx, dt_dy);
    //MHD also needs the transverse fluxes in the first ghost layer to calculate E
    computeFlux_X(_xL, _xR, F_X, 0, nx, -1, ny+1, dt_dx);
    computeFlux_Y(_yL, _yR, F_Y, -1, nx+1, 0, ny, dt_dy);
    #else
    ctu_sweep_hydro(_xL, _xR, _yL, _yR, _zL, _zR, dt_dx, dt_dy, dt_dz);
    //Hydro doesn't need transverse fluxes in the first ghost layer
    computeFlux_X(_xL, _xR, F_X, 0, nx, 0, ny, dt_dx);
    computeFlux_Y(_yL, _yR, F_Y, 0, nx, 0, ny, dt_dy);
    #endif
        __half_states.release();
    #endif
        #ifdef MHD
        __mag_half.release();
        #endif
    
    //Preliminarily apply all fluxes + sources
        auto __w = DRAGONWING::requestPrimitiveArrays(1, nx, ny, ghosts); //Auto-releases when the function terminates
    FluidArray2D& _w = *__w[0];
    applyFluxes(w, _w, F_X, F_Y, sources.non_stiff_terms(), _w_half, dt_dx, dt_dy, dt);

    
    //Preliminary CT Update
    #ifdef MHD
        auto __Elec = DRAGONWING::requestVec3Arrays(2, nx+1, ny+1, ghosts);
    //Compute Electric Fields
    MagneticArray2D& E = *__Elec[0], &_E_half = *__Elec[1];
    CT::computeElectric(E, F_X, F_Y);
    CT::bodyElectric(_w_half, _E_half);
    CT::upwindElectric(E, F_X, F_Y, _E_half);
    //Update B
        auto __B = DRAGONWING::requestVec3Arrays(1, nx+1, ny+1,  ghosts); //Auto-releases when the function terminates
    MagneticArray2D& _B = *__B[0];
    _B.clone(B);
    CT::Faraday(E, _B, dt_dx, dt_dy);
        __Elec.release();
    CT::computeBodyFields(_B, _w);
    #endif
        __w_half.release();

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
    
    //Compute Preliminary Fluxes
        auto __fluxes = DRAGONWING::requestFluxArrays(3, nx, ny, nz, ghosts);
    FluxArray3D& F_X = *__fluxes[0];
    FluxArray3D& F_Y = *__fluxes[1];
    FluxArray3D& F_Z = *__fluxes[2];
    #ifdef MHD //MHD needs an extra layer to calculate E
    computeFlux_X(_xL, _xR, F_X, -1, nx+1, -2, ny+2, -2, nz+2, dt_dx);
    computeFlux_Y(_yL, _yR, F_Y, -2, nx+2, -1, ny+1, -2, nz+2, dt_dy);
    computeFlux_Z(_zL, _zR, F_Z, -2, nx+2, -2, ny+2, -1, nz+1, dt_dz);
    #else //Hydro doesn't have electric fields
    computeFlux_X(_xL, _xR, F_X, 0, nx, -1, ny+1, -1, nz+1, dt_dx);
    computeFlux_Y(_yL, _yR, F_Y, -1, nx+1, 0, ny, -1, nz+1, dt_dy);
    computeFlux_Z(_zL, _zR, F_Z, -1, nx+1, -1, ny+1, 0, nz, dt_dz);
    #endif
        #ifndef CTU
        __half_states.release();
        #endif
    
    //Calculate half-time state (needed for source terms + constrained transport)
        auto __w_half = DRAGONWING::requestPrimitiveArrays(1, nx, ny, nz, ghosts);
    FluidArray3D& _w_half = *__w_half[0];
    applyFluxes(w, _w_half, F_X, F_Y, F_Z, sources.non_stiff_terms(), w, dt_dx*0.5, dt_dy*0.5, dt_dz*0.5, dt*0.5, 1);
    #ifdef MHD //Preliminary CT update
        auto __mag_half = DRAGONWING::requestVec3Arrays(2, nx+1, ny+1, nz+1, ghosts);
    MagneticArray3D &_B_half = *__mag_half[0], &_E = *__mag_half[1];
    //Compute Electric fields
    CT::computeElectric(_E, F_X, F_Y, F_Z,1);
    CT::bodyElectric(w, _B_half, 1); //Use _B_half buffer for Eref since we don't need it yet
    CT::upwindElectric(_E, F_X, F_Y, F_Z, _B_half,1);
    //Compute Magnetic Fields
    _B_half.clone(B);
    CT::Faraday(_E, _B_half, 0.5*dt_dx, 0.5*dt_dy, 0.5*dt_dz, 1);
    CT::computeBodyFields(_B_half, _w_half);
    #endif
    
    
    //CTU corrections & Final Fluxes
    #ifdef CTU //Gardiner and Stone (2008) https://arxiv.org/abs/0712.2634
    #ifdef MHD
    ctu_sweep_MHD(_xL, _xR, _yL, _yR, _zL, _zR, F_X, F_Y, F_Z, _E, _B_half, w, B, dt_dx, dt_dy, dt_dz);
    //MHD needs the transverse fluxes in the first ghost layer to calculate E
    computeFlux_X(_xL, _xR, F_X, 0, nx, -1, ny+1, -1, nz+1, dt_dx);
    computeFlux_Y(_yL, _yR, F_Y, -1, nx+1, 0, ny, -1, nz+1, dt_dy);
    computeFlux_Z(_zL, _zR, F_Z, -1, nx+1, -1, ny+1, 0, nz, dt_dz);
    #else
    ctu_sweep_hydro(_xL, _xR, _yL, _yR, _zL, _zR, F_X, F_Y, F_Z, dt_dx, dt_dy, dt_dz);
    //Hydro doesn't need transverse fluxes in the first ghost layer
    computeFlux_X(_xL, _xR, F_X, 0, nx, 0, ny, 0, nz, dt_dx);
    computeFlux_Y(_yL, _yR, F_Y, 0, nx, 0, ny, 0, nz, dt_dy);
    computeFlux_Z(_zL, _zR, F_Z, 0, nx, 0, ny, 0, nz, dt_dz);
    #endif
        __half_states.release();
    #endif
        #ifdef MHD
        __mag_half.release();
        #endif
    
    //Preliminarily apply all fluxes + sources
        auto __w = DRAGONWING::requestPrimitiveArrays(1, nx, ny, nz, ghosts); //Auto-releases when the function terminates
    FluidArray3D& _w = *__w[0];
    applyFluxes(w, _w, F_X, F_Y, F_Z, sources.non_stiff_terms(), _w_half, dt_dx, dt_dy, dt_dz, dt);

    
    //Preliminary CT Update
    #ifdef MHD
        auto __Elec = DRAGONWING::requestVec3Arrays(2, nx+1, ny+1, nz+1, ghosts);
    //Compute Electric Fields
    MagneticArray3D& E = *__Elec[0], &_E_half = *__Elec[1];
    CT::computeElectric(E, F_X, F_Y, F_Z);
    CT::bodyElectric(_w_half, _E_half);
    CT::upwindElectric(E, F_X, F_Y, F_Z, _E_half);
    //Update B
        auto __B = DRAGONWING::requestVec3Arrays(1, nx+1, ny+1, nz+1, ghosts); //Auto-releases when the function terminates
    MagneticArray3D& _B = *__B[0];
    _B.clone(B);
    CT::Faraday(E, _B, dt_dx, dt_dy, dt_dz);
        __Elec.release();
    CT::computeBodyFields(_B, _w);
    #endif
        __w_half.release();
    
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
