//
//  Unsplit.cpp
//  DRAGON/Hydro/Godunov
//
//  Created by Bobbie Markwick on 10/06/2026.
//  Implementation based in part on  Toro (2009). https://doi.org/10.1007/b79761
//      Colella (1990). https://doi.org/10.1016/0021-9991(90)90233-Q


#include "Grid.hpp"
#include "Riemann.hpp"
#include "Config.h"
#include "CFL.hpp"
#include "TVD.hpp"
#include "CT.hpp"
#include <cassert>
#include <iostream>
#include <stdexcept>
#include "DragonWing.hpp"

bool Grid::on_step_fail(const std::exception &e){
    return !DRAGONWING::requestRestart(e.what());
}


//MARK: 2D Unsplit Step
void computeFlux_X(const FluidArray2D& _L, const FluidArray2D& _R, FluxArray2D& F, int xL, int xR, int yL, int yR, double dt_dx);
void computeFlux_Y(const FluidArray2D& _L, const FluidArray2D& _R, FluxArray2D& F, int xL, int xR, int yL, int yR, double dt_dy);
void applyFluxes(const FluidArray2D& w, FluidArray2D& _w, const FluxArray2D& F_X, const FluxArray2D& F_Y,  double dt_dx, double dt_dy,  int g = 0);
void computeHalfStates_X(FluidArray2D& _L, const Grid2D& _W,  FluidArray2D& _R, double dt);
void computeHalfStates_Y( FluidArray2D& _L,const Grid2D& _W, FluidArray2D& _R, double dt);
void correctState(FluidArray2D& _L, FluidArray2D& _R, const FluxArray2D& F, double dt_dL, int dim);

void Grid2D::unsplit_step(double dt){
    int nx = w.getSizeX(), ny = w.getSizeY(), ghosts = w.getGhosts();
    
    if(!DRAGONWING::waitForRelease()) return;
    //Half States
        auto __half_states = DRAGONWING::requestPrimitiveArrays(4, nx, ny, ghosts);
    FluidArray2D& _xL = *__half_states[0];
    FluidArray2D& _xR = *__half_states[1];
    FluidArray2D& _yL = *__half_states[2];
    FluidArray2D& _yR = *__half_states[3];
    
    //Compute Half States
    computeHalfStates_X(_xL, (*this), _xR, dt);//_xLR needs (-1...nx, -1...ny), (-2...nx+1, -1...ny)for MHD
    computeHalfStates_Y(_yL, (*this), _yR, dt);//_yLR needs (-1...nx, -1...ny), (-1...nx, -2...ny+1) for MHD
    #ifdef MHD//Compute face-normal fields
        auto __B = DRAGONWING::requestVec3Arrays(1, nx+1, ny+1, ghosts);
    MagneticArray2D& B = *__B[0];
    CT::computeFaceFields(A, B, dx, dy);
    CT::copyFaceFields_X(_xL, B, _xR);
    CT::copyFaceFields_Y(_yL, B, _yR);
    #endif
    
    //Fluxes
        auto __fluxes = DRAGONWING::requestFluxArrays(2, nx, ny, ghosts);
    FluxArray2D& F_X = *__fluxes[0];
    FluxArray2D& F_Y = *__fluxes[1];
#ifdef CTU //Colella (1990). https://doi.org/10.1016/0021-9991(90)90233-Q
    //Compute preliminary Fluxes
    computeFlux_X(_xL, _xR, F_X, -1, nx+1, -2, ny+2, dt/dx);  //F_X needs (0...nx, -1...ny), (-1...nx+1, -1...ny)for MHD
    computeFlux_Y(_yL, _yR, F_Y, -2, nx+2, -1, ny+1, dt/dy);  //F_Y needs (-1...nx, 0...ny), (-1...nx, -1...ny+1) for MHD
    #ifdef MHD//Half-update A and compute the face-normal fields
        auto __CTU_A = DRAGONWING::requestVec3Arrays(2, nx+1, ny+1, ghosts);
    MagneticArray2D& _A_half = *__CTU_A[0];
    _A_half.clone(A);
    //Compute the new A
    MagneticArray2D& _E = *__CTU_A[1];
    CT::computeElectric(_E, F_X, F_Y, 1);
    CT::updatePotential(_A_half, _E, dt/2,1);
    //Update the B fields
    CT::computeFaceFields(_A_half, B, dx, dy);
        __CTU_A.release();
    #endif
    
    //Correct states
    correctState(_xL, _xR, F_Y, (0.5*dt/dy), 1); //_xLR needs (-1...nx, 0...ny-1), (-1...nx, -1...ny) for MHD
    correctState(_yL, _yR, F_X, (0.5*dt/dx), 0); //_yLR needs (0...nx-1, -1...ny), (-1...nx, -1...ny) for MHD
 
    #ifdef MHD //Apply half-step face-normal B
    CT::copyFaceFields_X(_xL, B, _xR);
    CT::copyFaceFields_Y(_yL, B, _yR);
    #endif
#endif
    //Compute Fluxes
    computeFlux_X(_xL, _xR, F_X, 0, nx, -1, ny+1, dt/dx); //F_X needs (0...nx, 0...ny-1), (0...nx, -1...ny) for MHD
    computeFlux_Y(_yL, _yR, F_Y, -1, nx+1, 0, ny, dt/dy); //F_Y needs (0...nx-1, 0...ny), (-1...nx, 0...ny) for MHD
        __half_states.release();
    
    //Preliminarily apply all fluxes
        auto __w = DRAGONWING::requestPrimitiveArrays(1,nx, ny, ghosts);
    FluidArray2D& _w = *__w[0];
    applyFluxes(w, _w, F_X, F_Y, dt/dx, dt/dy);

    //Preliminary CT Update
    #ifdef MHD
        auto __A = DRAGONWING::requestVec3Arrays(2, nx+1, ny+1, ghosts);
    MagneticArray2D& _A = *__A[0];
    _A.clone(A);
    //Compute the new A
    MagneticArray2D& E = *__A[1];
    CT::computeElectric(E, F_X, F_Y);
    CT::updatePotential(_A, E, dt);
    //Update the B fields
    CT::computeFaceFields(_A, B, dx, dy);
    CT::computeBodyFields(B, _w);
        __B.release();
    #endif
        __fluxes.release();


    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            if(!_w[i,j].isPhysical()) throw std::runtime_error(std::format("Unphysical state would be produced at ({},{})",i,j));
        }
    }
        
    //Wait for any parallel grids to finish
    DRAGONWING::reportCheckpoint1();
    if(!DRAGONWING::waitForCheckpoint1()) return;
    
    //Commit flux updates
    w.clone(_w, false);
    #ifdef MHD
    A.clone(_A, false);
    #endif
    
}

//MARK: 3D Unsplit Step
void computeFlux_X(const FluidArray3D& _L, const FluidArray3D& _R, FluxArray3D& F, int xL, int xR, int yL, int yR, int zL, int zR, double dt_dx);
void computeFlux_Y(const FluidArray3D& _L, const FluidArray3D& _R, FluxArray3D& F, int xL, int xR, int yL, int yR, int zL, int zR, double dt_dy);
void computeFlux_Z(const FluidArray3D& _L, const FluidArray3D& _R, FluxArray3D& F, int xL, int xR, int yL, int yR, int zL, int zR, double dt_dz);
void applyFluxes(const FluidArray3D& w, FluidArray3D& _w, const FluxArray3D& F_X, const FluxArray3D& F_Y, const FluxArray3D& F_Z, double dt_dx, double dt_dy, double dt_dz, int g=0);
void computeHalfStates_X(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, double dt);
void computeHalfStates_Y(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, double dt);
void computeHalfStates_Z(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, double dt);
void correctState(FluidArray3D& _L, FluidArray3D& _R, const FluxArray3D& F, double dt_dL, int dim);
void computeCTUFlux_X(const FluidArray3D& _L, const FluidArray3D& _R, const FluxArray3D& Ftrans, FluxArray3D& F, double dt_dx, double dt_dy, int dim);
void computeCTUFlux_Y(const FluidArray3D& _L, const FluidArray3D& _R, const FluxArray3D& Ftrans, FluxArray3D& F,  double dt_dy, double dt_dz, int dim);
void computeCTUFlux_Z(const FluidArray3D& _L, const FluidArray3D& _R, const FluxArray3D& Ftrans, FluxArray3D& F, double dt_dz, double dt_dy, int dim);


void Grid3D::unsplit_step(double dt){
    int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
    
    if(!DRAGONWING::waitForRelease()) return;
    //Half States
        auto __half_states = DRAGONWING::requestPrimitiveArrays(6, nx, ny, nz, ghosts);
    FluidArray3D& _xL = *__half_states[0];
    FluidArray3D& _xR = *__half_states[1];
    FluidArray3D& _yL = *__half_states[2];
    FluidArray3D& _yR = *__half_states[3];
    FluidArray3D& _zL = *__half_states[4];
    FluidArray3D& _zR = *__half_states[5];
           
    //Compute Half States
    computeHalfStates_X(_xL, (*this), _xR, dt);
    computeHalfStates_Y(_yL, (*this), _yR, dt);
    computeHalfStates_Z(_zL, (*this), _zR, dt);
    #ifdef MHD//Compute face-normal fields
        auto __B = DRAGONWING::requestVec3Arrays(1, nx+1, ny+1, nz+1, ghosts);
    MagneticArray3D& B = *__B[0];
    CT::computeFaceFields(A, B, dx, dy, dz);
    CT::copyFaceFields_X(_xL, B, _xR);
    CT::copyFaceFields_Y(_yL, B, _yR);
    CT::copyFaceFields_Z(_zL, B, _zR);
    #endif

    //Fluxes
        auto __fluxes = DRAGONWING::requestFluxArrays(3, nx, ny, nz, ghosts);
    FluxArray3D& F_X = *__fluxes[0];
    FluxArray3D& F_Y = *__fluxes[1];
    FluxArray3D& F_Z = *__fluxes[2];
#ifdef CTU //Colella (1990). https://doi.org/10.1016/0021-9991(90)90233-Q
    //Compute preliminary face fluxes
    computeFlux_X(_xL, _xR, F_X, -1, nx+1, -2, ny+2, -2, nz+2, dt/dx);
    computeFlux_Y(_yL, _yR, F_Y, -2, nx+2, -1, ny+1, -2, nz+2, dt/dy);
    computeFlux_Z(_zL, _zR, F_Z, -2, nx+2, -2, ny+2, -1, nz+1, dt/dz);
    
    
    #ifdef MHD//Half-update A and compute the face-normal fields
    //Do this now so CTU can borrow the flux aux grids
        auto __CTU_A = DRAGONWING::requestVec3Arrays(2, nx+1, ny+1, nz+1, ghosts);
    MagneticArray3D& _A_half = *__CTU_A[0];
    _A_half.clone(A);
    //Compute the new A
    MagneticArray3D& _E = *__CTU_A[1];
    CT::computeElectric(_E, F_X, F_Y, F_Z, 1);
    CT::updatePotential(_A_half, _E, dt/2,1);
    //Update the B fields
    CT::computeFaceFields(_A_half, B, dx, dy, dz);
        __CTU_A.release();
    #endif
        
    //Compute Edge-correct the fluxes
        auto __CTU_fluxes = DRAGONWING::requestFluxArrays(2, nx, ny, nz, ghosts);
    FluxArray3D& F_Xz = *__CTU_fluxes[0];
    computeCTUFlux_X(_xL, _xR, F_Z, F_Xz, dt/dx, (0.5*dt/dz), 2);
    FluxArray3D& F_Yz = *__CTU_fluxes[1];
    computeCTUFlux_Y(_yL, _yR, F_Z, F_Yz, dt/dy, (0.5*dt/dz), 2);

    FluxArray3D& F_Xy = F_Z; //F_Z isn't used again before it gets recomputed, so we can resue it here
    computeCTUFlux_X(_xL, _xR, F_Y, F_Xy, dt/dx, (0.5*dt/dy), 1);
    FluxArray3D& F_Zy = F_Y;  //In computeCTUFlux, the last read of [i,j,k] happens before [i,j,k] gets written, so we can likewise repurpose F_Y here
    computeCTUFlux_Z(_zL, _zR, F_Y, F_Zy, dt/dz, (0.5*dt/dy),1);

    //Update the X half states based on the YZ corner fluxes
    correctState(_xL, _xR, F_Yz, (0.5*dt/dy), 1);
    correctState(_xL, _xR, F_Zy, (0.5*dt/dz), 2);
    //Doing this early means we can reuse the F_Yz and F_Zy grids for F_Yx and F_Zx
    
    FluxArray3D& F_Yx = F_Yz; //Reuse existing grid that has already served its purpose
    computeCTUFlux_Y(_yL, _yR, F_X, F_Yx, dt/dy, (0.5*dt/dx), 0);
    FluxArray3D& F_Zx = F_Zy; //Reuse existing grid that has already served its purpose
    computeCTUFlux_Z(_zL, _zR, F_X, F_Zx, dt/dz, (0.5*dt/dx),0);

    //Update the Y half states based on the XZ corner fluxes
    correctState(_yL, _yR, F_Xz, (0.5*dt/dx), 0);
    correctState(_yL, _yR, F_Zx, (0.5*dt/dz), 2);
    //Update the Z half states based on the XY corner fluxes
    correctState(_zL, _zR, F_Xy, (0.5*dt/dx), 0);
    correctState(_zL, _zR, F_Yx, (0.5*dt/dy), 1);
        __CTU_fluxes.release();
    
    #ifdef MHD//Apply the half-step face-normal fields
    CT::copyFaceFields_X(_xL, B, _xR);
    CT::copyFaceFields_Y(_yL, B, _yR);
    CT::copyFaceFields_Z(_zL, B, _zR);
    #endif
#endif
    //Compute Fluxes
    computeFlux_X(_xL, _xR, F_X, 0, nx, -1, ny+1, -1, nz+1, dt/dx); //F_X needs (0...nx, -1...ny, -1...nz)
    computeFlux_Y(_yL, _yR, F_Y, -1, nx+1, 0, ny, -1, nz+1, dt/dy); //F_Y needs (-1...nx, 0...ny, -1...nz)
    computeFlux_Z(_zL, _zR, F_Z, -1, nx+1, -1, ny+1, 0, nz, dt/dz); //F_Z needs (-1...nx, -1...ny, 0...nz)
        __half_states.release();
    
    //Preliminarily apply all fluxes
        auto __w = DRAGONWING::requestPrimitiveArrays(1, nx, ny, nz, ghosts);
    FluidArray3D& _w = *__w[0];
    applyFluxes(w, _w, F_X, F_Y, F_Z, dt/dx, dt/dy, dt/dz);
    
    //Preliminary CT Update
    #ifdef MHD
        auto __mags = DRAGONWING::requestVec3Arrays(2, nx+1, ny+1, nz+1, ghosts);
    MagneticArray3D& _A = *__mags[0];
    _A.clone(A);
    //Compute the new A
    MagneticArray3D& E = *__mags[1];
    CT::computeElectric(E, F_X, F_Y, F_Z);
    CT::updatePotential(_A, E, dt);
    //Update the B fields
    CT::computeFaceFields(_A, B, dx, dy, dz);
    CT::computeBodyFields(B, _w);
        __B.release();
    #endif
        __fluxes.release();
    
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                if(!_w[i,j,k].isPhysical()) throw std::runtime_error(std::format("Unphysical state would be produced at ({},{},{})",i,j,k));
            }
        }
    }
    
    //Wait for any parallel grids to finish
    DRAGONWING::reportCheckpoint1();
    if(!DRAGONWING::waitForCheckpoint1()) return;
    
    //Commit Flux updates
    w.clone(_w, false);
    #ifdef MHD
    A.clone(_A, false);
    #endif
}

//MARK: Flux Calculation
//Compute X fluxes between Right(_R) and Left (_L) half-states over the entire grid
void computeFlux_X(const FluidArray2D& _L, const FluidArray2D& _R, FluxArray2D& F, int xL, int xR, int yL, int yR, double dt_dx ){
    for(int i=xL; i<=xR; i++){
        for(int j=yL; j<yR; j++){
            F[i,j] = Riemann(_R[i-1,j], _L[i,j]).flux_X(dt_dx);
        }
    }
}
//Compute Y fluxes between Right(_R) and Left (_L) half-states over the entire grid
void computeFlux_Y(const FluidArray2D& _L, const FluidArray2D& _R, FluxArray2D& F, int xL, int xR, int yL, int yR, double dt_dy){
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<=yR; j++){
            F[i,j] = Riemann(_R[i,j-1], _L[i,j]).flux_Y(dt_dy);
        }
    }
}
//Compute X fluxes between Right(_R) and Left (_L) half-states over the entire grid
void computeFlux_X(const FluidArray3D& _L, const FluidArray3D& _R, FluxArray3D& F, int xL, int xR, int yL, int yR, int zL, int zR, double dt_dx){
    for(int i=xL; i<=xR; i++){
        for(int j=yL; j<yR; j++){
            for(int k=zL; k<zR; k++){
                F[i,j,k] = Riemann(_R[i-1,j,k], _L[i,j,k]).flux_X(dt_dx);
            }
        }
    }
}
//Compute Y fluxes between Right(_R) and Left (_L) half-states over the entire grid
void computeFlux_Y(const FluidArray3D& _L, const FluidArray3D& _R, FluxArray3D& F, int xL, int xR, int yL, int yR, int zL, int zR, double dt_dy){
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<=yR; j++){
            for(int k=zL; k<zR; k++){
                F[i,j,k] = Riemann(_R[i,j-1,k], _L[i,j,k]).flux_Y(dt_dy);
            }
        }
    }
}
//Compute Z fluxes between Right(_R) and Left (_L) half-states over the entire grid
void computeFlux_Z(const FluidArray3D& _L, const FluidArray3D& _R, FluxArray3D& F, int xL, int xR, int yL, int yR, int zL, int zR, double dt_dz){
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<yR; j++){
            for(int k=zL; k<=zR; k++){
                F[i,j,k] = Riemann(_R[i,j,k-1], _L[i,j,k]).flux_Z(dt_dz);
            }
        }
    }
}
//MARK: Flux Application
void applyFluxes(const FluidArray2D& w, FluidArray2D& _w, const FluxArray2D& F_X, const FluxArray2D& F_Y,  double dt_dx, double dt_dy,  int g){
    const int nx = w.getSizeX(), ny = w.getSizeY();
    
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            ConservativeState U(w[i,j]);
            U += dt_dx * (F_X[i,j] - F_X[i+1,j]);
            U += dt_dy * (F_Y[i,j] - F_Y[i,j+1]);
            _w[i,j] = U;

            if(!U.isFinite()) throw std::runtime_error(std::format("\tNaN state would be produced at ({},{})\n",i,j));
        }
    }
}
void applyFluxes(const FluidArray3D& w, FluidArray3D& _w, const FluxArray3D& F_X, const FluxArray3D& F_Y, const FluxArray3D& F_Z, double dt_dx, double dt_dy, double dt_dz, int g){
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ();
    
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g; k<nz+g; k++){
                ConservativeState U(w[i,j,k]);
                U += dt_dx * (F_X[i,j,k] - F_X[i+1,j,k]);
                U += dt_dy * (F_Y[i,j,k] - F_Y[i,j+1,k]);
                U += dt_dz * (F_Z[i,j,k] - F_Z[i,j,k+1]);
                _w[i,j,k] = U;

                if(!U.isFinite()) throw std::runtime_error(std::format("\tNaN state would be produced at ({},{},{})\n",i,j,k));
            }
        }
    }
}

//MARK: MUSCL
//Apply MUSCL over the entire grid
void computeHalfStates_X(FluidArray2D& _L, const Grid2D& _W, FluidArray2D& _R, double dt){
    const double dt_dL = dt/_W.dx;//Compute once
    const int nx = _W.getSizeX(), ny = _W.getSizeY(), g = _W.getGhosts();
    //MUSCL Reconstruction
    for(int i=-g+1; i<nx+g - 1; i++){
        for(int j=-g; j<ny+g; j++){
            auto wL =_W[i-1,j], wR = _W[i+1,j];
            TVD::MUSCL(wL, _L[i,j], _W[i,j], _R[i,j], wR, dt_dL);
        }
    }
    //Leftmost and rightmost ghosts
    for(int j=-g; j<ny+g; j++){
        _L[-g,j] = _W[-g,j]; _R[-g,j] = _W[-g,j];
        _L[nx-1+g,j] = _W[nx-1+g,j]; _R[nx-1+g,j] = _W[nx-1+g,j];
    }
}
void computeHalfStates_Y(FluidArray2D& _L, const Grid2D& _W, FluidArray2D& _R, double dt){
    const double dt_dL = dt/_W.dy;//Compute once
    const int nx = _W.getSizeX(), ny = _W.getSizeY(), g = _W.getGhosts();

    //MUSCL Reconstruction
    for(int i=-g; i<nx+g; i++){
        for(int j=-g + 1; j<ny+g - 1; j++){
            //Swap XY inputs to MUSCL, then swap output back
            auto wL =_W[i,j-1].swappedXY(), wR = _W[i,j+1].swappedXY();
            TVD::MUSCL(wL, _L[i,j], _W[i,j].swappedXY(), _R[i,j],wR, dt_dL);
            _L[i,j].swapXY(); _R[i,j].swapXY();
        }
    }
    //Leftmost and rightmost ghosts
    for(int i=-g; i<nx+g; i++){
        _L[i,-g] = _W[i,-g]; _R[i,-g] = _W[i,-g];
        _L[i,ny-1+g] = _W[i,ny-1+g]; _R[i,ny-1+g] = _W[i,ny-1+g];
    }
}
void computeHalfStates_X(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, double dt){
    const double dt_dL = dt/_W.dx;//Compute once
    const int nx = _W.getSizeX(), ny = _W.getSizeY(),nz = _W.getSizeZ(), g = _W.getGhosts();

    //MUSCL Reconstruction
    for(int i=-g + 1; i<nx+g - 1; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g; k<nz+g; k++){
                auto wL = _W[i-1,j,k], wR = _W[i+1,j,k];
                TVD::MUSCL(wL, _L[i,j,k], _W[i,j,k], _R[i,j,k], wR, dt_dL);
            }
        }
    }
    //Leftmost and rightmost ghosts
    for(int j=-g; j<ny+g; j++){
        for(int k=-g; k<nz+g; k++){
            _L[-g,j,k] = _W[-g,j,k]; _R[-g,j,k] = _W[-g,j,k];
            _L[nx-1+g,j,k] = _W[nx-1+g,j,k]; _R[nx-1+g,j,k] = _W[nx-1+g,j,k];
        }
    }
}
void computeHalfStates_Y(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, double dt){
    const double dt_dL = dt/_W.dy;//Compute once
    const int nx = _W.getSizeX(), ny = _W.getSizeY(),nz = _W.getSizeZ(), g = _W.getGhosts();
    
    //MUSCL Reconstruction
    for(int i=-g; i<nx+g; i++){
        for(int j=-g + 1; j<ny+g - 1; j++){
            for(int k=-g; k<nz+g; k++){
                //Swap XY inputs to MUSCL, then swap output back
                auto wL =_W[i,j-1,k].swappedXY(), wR = _W[i,j+1,k].swappedXY();
                TVD::MUSCL(wL, _L[i,j,k], _W[i,j,k].swappedXY(), _R[i,j,k],wR, dt_dL);
                _L[i,j,k].swapXY(); _R[i,j,k].swapXY();
            }
        }
    }
    //Leftmost and rightmost ghosts
    for(int i=-g; i<nx+g; i++){
        for(int k=-g; k<nz+g; k++){
            _L[i,-g,k] = _W[i,-g,k]; _R[i,-g,k] = _W[i,-g,k];
            _L[i,ny-1+g,k] = _W[i,ny-1+g,k]; _R[i,ny-1+g,k] = _W[i,ny-1+g,k];
        }
    }
}
void computeHalfStates_Z(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, double dt){
    const double dt_dL = dt/_W.dz;//Compute once
    const int nx = _W.getSizeX(), ny = _W.getSizeY(),nz = _W.getSizeZ(), g = _W.getGhosts();
    
    //MUSCL Reconstruction
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g + 1; k<nz+g - 1; k++){
                //Swap XZ inputs to MUSCL, then swap output back
                auto wL = _W[i,j,k-1].swappedXZ(), wR = _W[i,j,k+1].swappedXZ();
                TVD::MUSCL(wL, _L[i,j,k], _W[i,j,k].swappedXZ(), _R[i,j,k],wR, dt_dL);
                _L[i,j,k].swapXZ(); _R[i,j,k].swapXZ();
            }
        }
    }
    //Leftmost and rightmost ghosts
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            _L[i,j,-g] = _W[i,j,-g]; _R[i,j,-g] = _W[i,j,-g];
            _L[i,j,nz-1+g] = _W[i,j,nz-1+g]; _R[i,j,nz-1+g] = _W[i,j,nz-1+g];
        }
    }
}

//MARK: CTU Corrections
#ifdef CTU
#ifndef MUSCL_Hancock
#error CTU requires MUSCL reconstruction. Please enable MUSCL_Hancock in Config.h
#endif
#endif

void correctState(FluidArray2D& _L, FluidArray2D& _R, const FluxArray2D& F, double dt_dL, int dim){
    const int xL = -1, xR = F.getSizeX()+1, yL = -1, yR = F.getSizeY()+1;
    //Extract direction encoding
    int isX = dim%2 == 0 ? 1 : 0;
    int isY = dim%2 == 1 ? 1 : 0;
    //Cycle
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<yR; j++){
            auto trans = (F[i+isX,j+isY] - F[i,j]) * dt_dL;
            _L[i,j] = _L[i,j] -  trans;
            _R[i,j] = _R[i,j] -  trans;
        }
    }
}
void correctState(FluidArray3D& _L, FluidArray3D& _R, const FluxArray3D& F, double dt_dL, int dim){
    const int xL = -1, xR = F.getSizeX()+1, yL = -1, yR = F.getSizeY()+1, zL = -1, zR = F.getSizeZ()+1;
    //Extract direction encoding
    int isX = dim%3 == 0 ? 1 : 0;
    int isY = dim%3 == 1 ? 1 : 0;
    int isZ = dim%3 == 2 ? 1 : 0;
    //cycle
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<yR; j++){
            for(int k=zL; k<zR; k++){
                auto trans = (F[i+isX, j+isY, k+isZ] - F[i,j,k]) * dt_dL;
                _L[i,j,k] = _L[i,j,k] -  trans;
                _R[i,j,k] = _R[i,j,k] -  trans;
            }
        }
    }
}

//Compute X fluxes between Right(_R) and Left (_L) half-states, as corrected by transverse fluxes FYZ
//Equivalent to correctState -> computeFlux_X but without needing intermediate arrays
void computeCTUFlux_X(const FluidArray3D& _L, const FluidArray3D& _R, const FluxArray3D& FYZ, FluxArray3D& F,  double dt_dx, double dt_dy, int dim){
    const int xL = -1, xR = F.getSizeX()+1, yL = -1, yR = F.getSizeY()+1, zL = -1, zR = F.getSizeZ()+1;
    //Extract direction encoding
    int isY = dim%3 == 1 ? 1 : 0;
    int isZ = dim%3 == 2 ? 1 : 0;
    //cycle
    for(int j=yL; j<yR; j++){
        for(int k=zL; k<zR; k++){
            //Compute Left half-state on the first cell
            PrimitiveState  __R, __L =  _R[xL-1, j,k] - (FYZ[xL-1, j+isY, k+isZ] - FYZ[xL-1,j,k]) * dt_dy;
            for(int i=xL; i<=xR; i++){
                auto trans = (FYZ[i, j+isY, k+isZ] - FYZ[i,j,k]) * dt_dy;
                __R = _L[i,j,k] -  trans;//Compute Right half-state for this flux
                F[i,j,k] = Riemann(__L, __R).flux_X(dt_dx);
                __L = _R[i,j,k] -  trans;//Compute Left half-state for the next flux
            }
        }
    }
}
//Compute Y fluxes between Right(_R) and Left (_L) half-states, as corrected by transverse fluxes FXZ
//Equivalent to correctState -> computeFlux_Y but without needing intermediate arrays
void computeCTUFlux_Y(const FluidArray3D& _L, const FluidArray3D& _R, const FluxArray3D& FXZ, FluxArray3D& F, double dt_dy, double dt_dz, int dim){
    const int xL = -1, xR = F.getSizeX()+1, yL = -1, yR = F.getSizeY()+1, zL = -1, zR = F.getSizeZ()+1;

    //Extract direction encoding
    int isX = dim%3 == 0 ? 1 : 0;
    int isZ = dim%3 == 2 ? 1 : 0;
    //cycle
    for(int i=xL; i<xR; i++){
        for(int k=zL; k<zR; k++){
            //Compute Left half-state on the first cell
            PrimitiveState  __R, __L =  _R[i,yL-1,k] - (FXZ[i+isX,yL-1, k+isZ] - FXZ[i,yL-1,k]) * dt_dz;
            for(int j=yL; j<=yR; j++){
                auto trans = (FXZ[i+isX, j, k+isZ] - FXZ[i,j,k]) * dt_dz;
                __R = _L[i,j,k] -  trans;//Compute Right half-state for this flux
                F[i,j,k] = Riemann(__L, __R).flux_Y(dt_dy);
                __L = _R[i,j,k] -  trans;//Compute Left half-state for the next flux
            }
        }
    }
}
//Compute Z fluxes between Right(_R) and Left (_L) half-states, as corrected by transverse fluxes FXY
//Equivalent to correctState -> computeFlux_Z but without needing intermediate arrays
void computeCTUFlux_Z(const FluidArray3D& _L, const FluidArray3D& _R, const FluxArray3D& FXY, FluxArray3D& F, double dt_dz, double dt_dy, int dim){
    const int xL = -1, xR = F.getSizeX()+1, yL = -1, yR = F.getSizeY()+1, zL = -1, zR = F.getSizeZ()+1;
    //Extract direction encoding
    int isX = dim%3 == 0 ? 1 : 0;
    int isY = dim%3 == 1 ? 1 : 0;
    //cycle
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<yR; j++){
            //Compute Left half-state on the first cell
            PrimitiveState  __R, __L =  _R[i,j,zL-1] - (FXY[i+isX,j+isY,zL-1] - FXY[i,j,zL-1]) * dt_dy;
            for(int k=zL; k<=zR; k++){
                auto trans = (FXY[i+isX, j+isY, k] - FXY[i,j,k]) * dt_dy;
                __R = _L[i,j,k] -  trans;//Compute Right half-state for this flux
                F[i,j,k] = Riemann(__L, __R).flux_Z(dt_dz);
                __L = _R[i,j,k] -  trans;//Compute Left half-state for the next flux
            }
        }
    }
}
