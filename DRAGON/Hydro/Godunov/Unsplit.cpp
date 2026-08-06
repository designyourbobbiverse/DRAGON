//
//  Unsplit.cpp
//  DRAGON/Hydro/Godunov
//
//  Created by Bobbie Markwick on 10/06/2026.
//  Implementation based in part on  Toro (2009). https://doi.org/10.1007/b79761


#include "Grid.hpp"
#include "Godunov.hpp"

#include "Config.h"
#include "CT.hpp"      //For MHD

#include "DragonWing.hpp" //For memory management & synchronization
#include <stdexcept>      //For error handling
#include <format>         //For error message formatting

using namespace Godunov;

bool Grid::on_step_fail(const std::exception &e){
    return !DRAGONWING::requestRestart(e.what());
}


//MARK: 2D Unsplit Step

void Grid2D::unsplit_step(double dt){
    const int nx = w.getSizeX(), ny = w.getSizeY(), ghosts = w.getGhosts();
    const double dt_dx = dt/dx, dt_dy = dt/dy;
    
    if(!DRAGONWING::waitForRelease()) return;
    #ifndef MHD //Dummy B array
        auto B = MagneticArray2D(0,0);
    #endif
    
    //Compute Half States
        auto __half_states = DRAGONWING::requestPrimitiveArrays(4, nx, ny, ghosts);
    FluidArray2D& _xL = *__half_states[0];
    FluidArray2D& _xR = *__half_states[1];
    FluidArray2D& _yL = *__half_states[2];
    FluidArray2D& _yR = *__half_states[3];
    computeHalfStates_X(_xL, (*this), _xR, B, dt);
    computeHalfStates_Y(_yL, (*this), _yR, B, dt);
    
    #ifdef CTU
        #ifdef MHD //Gardiner and Stone (2005). https://doi.org/10.1016/j.jcp.2004.11.016
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
        auto __w = DRAGONWING::requestPrimitiveArrays(1,nx, ny, ghosts);
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
        auto __B = DRAGONWING::requestVec3Arrays(1, nx+1, ny+1, ghosts);
    MagneticArray2D& _B = *__B[0]; //Done with _E_half, reuse it
    _B.clone(B);
    CT::Faraday(E, _B, dt_dx, dt_dy);
        __Elec.release();
    CT::computeBodyFields(_B, _w);
    #endif
        __fluxes.release();

    //Verify Physicality of solution
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            if(!_w[i,j].isPhysical())
                throw std::runtime_error(std::format("Unphysical state would be produced at ({},{})",i,j));
        }
    }
        
    //Wait for any parallel grids to finish
    DRAGONWING::reportCheckpoint1();
    if(!DRAGONWING::waitForCheckpoint1()) return;
    
    //Commit flux updates
    w.clone(_w, false);
    #ifdef MHD
    B.clone(_B, false);
    #endif
    
}

//MARK: 3D Unsplit Step

void Grid3D::unsplit_step(double dt){
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
    const double dt_dx = dt/dx, dt_dy = dt/dy, dt_dz = dt/dz;
    
    if(!DRAGONWING::waitForRelease()) return;
    #ifndef MHD //Face Fields Dummy array
        auto B = MagneticArray3D(0,0,0);
    #endif
    
    //Compute Half States
        auto __half_states = DRAGONWING::requestPrimitiveArrays(6, nx, ny, nz, ghosts);
    FluidArray3D& _xL = *__half_states[0];
    FluidArray3D& _xR = *__half_states[1];
    FluidArray3D& _yL = *__half_states[2];
    FluidArray3D& _yR = *__half_states[3];
    FluidArray3D& _zL = *__half_states[4];
    FluidArray3D& _zR = *__half_states[5];
    computeHalfStates_X(_xL, (*this), _xR, B, dt);
    computeHalfStates_Y(_yL, (*this), _yR, B, dt);
    computeHalfStates_Z(_zL, (*this), _zR, B, dt);
    
    #ifdef CTU
        #ifdef MHD //Gardiner and Stone (2005). https://doi.org/10.1016/j.jcp.2004.11.016
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
        auto __w = DRAGONWING::requestPrimitiveArrays(1, nx, ny, nz, ghosts);
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
        auto __B = DRAGONWING::requestVec3Arrays(1, nx+1, ny+1, nz+1, ghosts);
    MagneticArray3D& _B = *__B[0]; //Done with _E_half, reuse it
    _B.clone(B);
    CT::Faraday(E, _B, dt_dx, dt_dy, dt_dz);
        __Elec.release();
    CT::computeBodyFields(_B, _w);
    #endif
        __fluxes.release();

    
    //Check Physicality
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                if(!_w[i,j,k].isPhysical())
                    throw std::runtime_error(std::format("Unphysical state would be produced at ({},{},{})",i,j,k));
            }
        }
    }
    
    //Wait for any parallel grids to finish
    DRAGONWING::reportCheckpoint1();
    if(!DRAGONWING::waitForCheckpoint1()) return;
    
    //Commit Flux updates
    w.clone(_w, false);
    #ifdef MHD
    B.clone(_B, false);
    #endif
}
