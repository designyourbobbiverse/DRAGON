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


//MARK: Split vs Unsplit
//Use split or unsplit depending on whether DIMENSION_UNSPLIT is defined in Config.h
void Grid2D::advance(double dt,bool check_cfl){
#ifdef DIMENSION_UNSPLIT
    advance_unsplit(dt,check_cfl);
#else
    advance_split(dt,check_cfl);
#endif
}
void Grid3D::advance(double dt, bool check_cfl){
#ifdef DIMENSION_UNSPLIT
    advance_unsplit(dt,check_cfl);
#else
    advance_split(dt,check_cfl);
#endif
}


//MARK: CFL
#ifdef MHD
typedef ExtendedArray2D<vec3> MagneticArray2D;
#endif
void Grid2D::advance_unsplit(double dt, bool check_cfl){
    #ifdef MHD
    initialize_B_fields();
    #endif
    while(dt > CONFIG::Timestep_Tolerance){
        //CFL Time Constraint
        double t1 = check_cfl ? std::min(dt,CFL::cfl_time(*this)) : dt;
        //Advance
        bool success = false;
        do{
            try{
                advanceXY(t1);
                success = true;
            } catch(const std::exception &exc){
                if(DRAGONWING::requestRestart(exc.what())){ return; }
                else std::cout<<"\t"<<exc.what()<<"\n";
            }
            if(!success) t1 *= 0.5;
        } while(!success);
        
        dt -= t1;
    }
}
void Grid3D::advance_unsplit(double dt, bool check_cfl){
    #ifdef MHD //Ensure B is initialised
    initialize_B_fields();
    #endif
    while(dt > CONFIG::Timestep_Tolerance){
        //CFL Time Constraint
        double t1 = check_cfl ? std::min(dt,CFL::cfl_time(*this)) : dt;
        //Advance
        bool success = false;
        do{
            try{
                advanceXYZ(t1);
                success = true;
            } catch(const std::exception &exc){
                if(DRAGONWING::requestRestart(exc.what())){ return; }
                else std::cout<<"\t"<<exc.what()<<"\n";
            }
            if(!success) t1 *= 0.5;
        } while(!success);
        
        dt -= t1;
    }
}

//MARK: 2D Unsplit Step
void computeFlux_X(const FluidArray2D& _L, const FluidArray2D& _R, FluxArray2D& F, int xL, int xR, int yL, int yR, double dt_dx);
void computeFlux_Y(const FluidArray2D& _L, const FluidArray2D& _R, FluxArray2D& F, int xL, int xR, int yL, int yR, double dt_dy);
void computeHalfStates_X(FluidArray2D& _L, const Grid2D& _W,  FluidArray2D& _R, double dt);
void computeHalfStates_Y( FluidArray2D& _L,const Grid2D& _W, FluidArray2D& _R, double dt);
void correctState(FluidArray2D& _L, FluidArray2D& _R, const FluxArray2D& F, double dt_dL, int xL, int xR, int yL, int yR, int dim);

void Grid2D::advanceXY(double dt){
    int nx = w.getSizeX(), ny = w.getSizeY(), ghosts = w.getGhosts();
    boundary.apply(*this);
    
    if(!DRAGONWING::waitForRelease()) return;
    //Half States
    FluidArray2D& _xL = *DRAGONWING::requestPrimitiveArray(nx, ny, ghosts);
    FluidArray2D& _xR = *DRAGONWING::requestPrimitiveArray(nx, ny, ghosts);
    FluidArray2D& _yL = *DRAGONWING::requestPrimitiveArray(nx, ny, ghosts);
    FluidArray2D& _yR = *DRAGONWING::requestPrimitiveArray(nx, ny, ghosts);
    //Fluxes
    FluxArray2D& F_X = *DRAGONWING::requestFluxArray(nx, ny, ghosts);
    FluxArray2D& F_Y = *DRAGONWING::requestFluxArray(nx, ny, ghosts);
    
    //Compute Half States
    computeHalfStates_X(_xL, (*this), _xR, dt);//_xLR needs (-1...nx, -1...ny), (-2...nx+1, -1...ny)for MHD
    computeHalfStates_Y(_yL, (*this), _yR, dt);//_yLR needs (-1...nx, -1...ny), (-1...nx, -2...ny+1) for MHD
    
    #ifdef MHD//Compute face-normal fields
    MagneticArray2D& B = *DRAGONWING::requestVec3Array(nx+1, ny+1, ghosts);
    CT::computeFaceFields(A, B, dx, dy);
    CT::copyFaceFields_X(_xL, B, _xR);
    CT::copyFaceFields_Y(_yL, B, _yR);
    #endif
    
#ifdef CTU //Colella (1990). https://doi.org/10.1016/0021-9991(90)90233-Q
    //Compute preliminary Fluxes
    computeFlux_X(_xL, _xR, F_X, -1, nx+1, -1, ny+1, dt/dx);  //F_X needs (0...nx, -1...ny), (-1...nx+1, -1...ny)for MHD
    computeFlux_Y(_yL, _yR, F_Y, -1, nx+1, -1, ny+1, dt/dy);  //F_Y needs (-1...nx, 0...ny), (-1...nx, -1...ny+1) for MHD
    //Correct states
    correctState(_xL, _xR, F_Y, (0.5*dt/dy), -1, nx+1, -1, ny+1, 1); //_xLR needs (-1...nx, 0...ny-1), (-1...nx, -1...ny) for MHD
    correctState(_yL, _yR, F_X, (0.5*dt/dx), -1, nx+1, -1, ny+1, 0); //_yLR needs (0...nx-1, -1...ny), (-1...nx, -1...ny) for MHD
 
    #ifdef MHD//Restore face-normal fields
    CT::copyFaceFields_X(_xL, B, _xR);
    CT::copyFaceFields_Y(_yL, B, _yR);
    #endif
#endif

    
    
    //Compute Fluxes
    computeFlux_X(_xL, _xR, F_X, 0, nx, -1, ny+1, dt/dx); //F_X needs (0...nx, 0...ny-1), (0...nx, -1...ny) for MHD
    computeFlux_Y(_yL, _yR, F_Y, -1, nx+1, 0, ny, dt/dy); //F_Y needs (0...nx-1, 0...ny), (-1...nx, 0...ny) for MHD
    
    FluidArray2D& _w = _xL; //Repurpose grid that isn't being used anymore
        DRAGONWING::releaseArray(&_xR); //Release the others
        DRAGONWING::releaseArray(&_yL);
        DRAGONWING::releaseArray(&_yR);
    
    //Preliminarily apply all fluxes
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            ConservativeState  U(w[i,j]);
            U += (dt/dx) * (F_X[i,j] - F_X[i+1,j]);
            U += (dt/dy) * (F_Y[i,j] - F_Y[i,j+1]);
            _w[i,j] = U;
            
            if(!U.isFinite()){
                    DRAGONWING::releaseArray(&F_X);
                    DRAGONWING::releaseArray(&F_Y);
                    DRAGONWING::releaseArray(&_w);
                    #ifdef MHD
                    DRAGONWING::releaseArray(&B);
                    #endif
                throw std::runtime_error(std::format("NaN state would be produced at ({},{})",i,j));
            }
        }
    }

    //Preliminary CT Update
    #ifdef MHD
    MagneticArray2D& _A = *DRAGONWING::requestVec3Array(nx+1, ny+1, ghosts);
    _A.clone(A);
    CT::updatePotential(_A, F_X, F_Y, dt);
    CT::computeFaceFields(_A, B, dx, dy);
    CT::computeBodyFields(B, _w);
        DRAGONWING::releaseArray(&B);
    #endif
        DRAGONWING::releaseArray(&F_X);
        DRAGONWING::releaseArray(&F_Y);


    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            if(!_w[i,j].isPhysical()){
                    DRAGONWING::releaseArray(&_w);
                    #ifdef MHD
                    DRAGONWING::releaseArray(&_A);
                    #endif
                throw std::runtime_error(std::format("Unphysical state would be produced by CT at ({},{})",i,j));
            }
        }
    }
        
    //Wait for any parallel grids to finish
    DRAGONWING::reportCheckpoint1();
    if(!DRAGONWING::waitForCheckpoint1()) return;
    
    //Commit flux updates
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            w[i,j] = _w[i,j];
        }
    }
        DRAGONWING::releaseArray(&_w);
    #ifdef MHD
    A.clone(_A);
        DRAGONWING::releaseArray(&_A);
    #endif
    
}

//MARK: 3D Unsplit Step
typedef ExtendedArray3D<PrimitiveState> FluidArray3D;
typedef ExtendedArray3D<ConservativeState> FluxArray3D;
void computeFlux_X(const FluidArray3D& _L, const FluidArray3D& _R, FluxArray3D& F, int xL, int xR, int yL, int yR, int zL, int zR, double dt_dx);
void computeFlux_Y(const FluidArray3D& _L, const FluidArray3D& _R, FluxArray3D& F, int xL, int xR, int yL, int yR, int zL, int zR, double dt_dy);
void computeFlux_Z(const FluidArray3D& _L, const FluidArray3D& _R, FluxArray3D& F, int xL, int xR, int yL, int yR, int zL, int zR, double dt_dz);
void computeHalfStates_X(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, double dt);
void computeHalfStates_Y(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, double dt);
void computeHalfStates_Z(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, double dt);
void correctState(const FluidArray3D& _L0, const FluidArray3D& _R0, FluidArray3D& _L, FluidArray3D& _R, const FluxArray3D& F, double dt_dL, int xL, int xR, int yL, int yR, int zL, int zR, int dim);


void Grid3D::advanceXYZ(double dt){
    int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
    boundary.apply(*this);

    if(!DRAGONWING::waitForRelease()) return;
    //Half States
    FluidArray3D& _xL = *DRAGONWING::requestPrimitiveArray(nx, ny, nz, ghosts);
    FluidArray3D& _xR = *DRAGONWING::requestPrimitiveArray(nx, ny, nz, ghosts);
    FluidArray3D& _yL = *DRAGONWING::requestPrimitiveArray(nx, ny, nz, ghosts);
    FluidArray3D& _yR = *DRAGONWING::requestPrimitiveArray(nx, ny, nz, ghosts);
    FluidArray3D& _zL = *DRAGONWING::requestPrimitiveArray(nx, ny, nz, ghosts);
    FluidArray3D& _zR = *DRAGONWING::requestPrimitiveArray(nx, ny, nz, ghosts);
    //Fluxes
    FluxArray3D& F_X = *DRAGONWING::requestFluxArray(nx, ny, nz, ghosts);
    FluxArray3D& F_Y = *DRAGONWING::requestFluxArray(nx, ny, nz, ghosts);
    FluxArray3D& F_Z = *DRAGONWING::requestFluxArray(nx, ny, nz, ghosts);
           
    //Compute Half States
    computeHalfStates_X(_xL, (*this), _xR, dt);
    computeHalfStates_Y(_yL, (*this), _yR, dt);
    computeHalfStates_Z(_zL, (*this), _zR, dt);
    
    #ifdef MHD//Compute face-normal fields
    MagneticArray3D& B = *DRAGONWING::requestVec3Array(nx+1, ny+1, nz+1, ghosts);
    CT::computeFaceFields(A, B, dx, dy, dz);
    CT::copyFaceFields_X(_xL, B, _xR);
    CT::copyFaceFields_Y(_yL, B, _yR);
    CT::copyFaceFields_Z(_zL, B, _zR);
    #endif

#ifdef CTU //Colella (1990). https://doi.org/10.1016/0021-9991(90)90233-Q
    //Compute preliminary face fluxes
    computeFlux_X(_xL, _xR, F_X, -1, nx+1, -2, ny+2, -2, nz+2, dt/dx); //F_X needs (0...nx, -1...ny, -1...nz)
    computeFlux_Y(_yL, _yR, F_Y, -2, nx+2, -1, ny+1, -2, nz+2, dt/dy); //F_Y needs (-1...nx, 0...ny, -1...nz)
    computeFlux_Z(_zL, _zR, F_Z, -2, nx+2, -2, ny+2, -1, nz+1, dt/dz); //F_Z needs (-1...nx, -1...ny, 0...nz)
    
    FluidArray3D& __L = *DRAGONWING::requestPrimitiveArray(nx, ny, nz, ghosts);
    FluidArray3D& __R = *DRAGONWING::requestPrimitiveArray(nx, ny, nz, ghosts);
    
    //Compute Edge-correct the fluxes
    //F_Xy (needs 0...nx, 0...ny-1, -1...nz) from _xyLR (needs -1...nx, 0...ny-1, -1...nz)
    FluxArray3D& F_Xy = *DRAGONWING::requestFluxArray(nx, ny, nz, ghosts);
    correctState(_xL, _xR, __L, __R, F_Y, (0.5*dt/dy) , -2, nx+2, -1, ny+1, -1, nz+1, 1);
    computeFlux_X(__L, __R, F_Xy, -1, nx+1, -1, ny+1, -1, nz+1, dt/dx);
    //F_Xz (needs 0...nx, -1...ny, 0...nz-1) from _xzLR (needs -1...nx, -1...ny, 0...nz-1)
    FluxArray3D& F_Xz = *DRAGONWING::requestFluxArray(nx, ny, nz, ghosts);
    correctState(_xL, _xR, __L, __R, F_Z, (0.5*dt/dz) , -2, nx+2, -1, ny+1, -1, nz+1, 2);
    computeFlux_X(__L, __R, F_Xz, -1, nx+1, -1, ny+1, -1, nz+1, dt/dx);
    //F_Yx (needs 0...nx-1, 0...ny, -1...nz) from _yxLR (needs 0...nx-1, -1...ny, -1...nz)
    FluxArray3D& F_Yx = *DRAGONWING::requestFluxArray(nx, ny, nz, ghosts);
    correctState(_yL, _yR, __L, __R, F_X, (0.5*dt/dx), -1, nx+1, -2, ny+2, -1, nz+1, 0);
    computeFlux_Y(__L, __R, F_Yx, -1, nx+1, -1, ny+1, -1, nz+1, dt/dy);
    //F_Yz (needs -1...nx, 0...ny, 0...nz-1) from _yzLR (needs -1...nx, -1...ny, 0...nz-1)
    FluxArray3D& F_Yz = *DRAGONWING::requestFluxArray(nx, ny, nz, ghosts);
    correctState(_yL, _yR, __L, __R, F_Z, (0.5*dt/dz), -1, nx+1, -2, ny+2, -1, nz+1, 2);
    computeFlux_Y(__L, __R, F_Yz, -1, nx+1, -1, ny+1, -1, nz+1, dt/dy);
    //F_Zx (needs 0...nx-1, -1...ny, 0...nz) from _zxLR (needs 0...nx-1, -1...ny, -1...nz)
    FluxArray3D& F_Zx = *DRAGONWING::requestFluxArray(nx, ny, nz, ghosts);
    correctState(_zL, _zR, __L, __R, F_X, (0.5*dt/dx), -1, nx+1, -1, ny+1, -2, nz+2, 0);
    computeFlux_Z(__L, __R, F_Zx, -1, nx+1, -1, ny+1, -1, nz+1, dt/dz);
    //F_Zy (needs -1...nx, 0...ny-1, 0...nz) from _zyLR (needs -1...nx, 0...ny-1, -1...nz)
    FluxArray3D& F_Zy = *DRAGONWING::requestFluxArray(nx, ny, nz, ghosts);
    correctState(_zL, _zR, __L, __R, F_Y, (0.5*dt/dy), -1, nx+1, -1, ny+1, -2, nz+2, 1);
    computeFlux_Z(__L, __R, F_Zy, -1, nx+1, -1, ny+1, -1, nz+1, dt/dz);
    
        DRAGONWING::releaseArray(&__L);
        DRAGONWING::releaseArray(&__R);
    
    //Apply second round of transverse corrections
    //_xLR needs (-1...nx, -1...ny, -1...nz)
    correctState(_xL, _xR, _xL,_xR, F_Yz, (0.5*dt/dy), -1, nx+1, -1, ny+1, -1, nz+1, 1);
        DRAGONWING::releaseArray(&F_Yz);
    correctState(_xL, _xR, _xL,_xR, F_Zy, (0.5*dt/dz), -1, nx+1, -1, ny+1, -1, nz+1, 2);
        DRAGONWING::releaseArray(&F_Zy);
    //_yLR needs (-1...nx, -1...ny, -1...nz)
    correctState(_yL, _yR, _yL,_yR, F_Xz, (0.5*dt/dx), -1, nx+1, -1, ny+1, -1, nz+1, 0);
        DRAGONWING::releaseArray(&F_Xz);
    correctState(_yL, _yR, _yL,_yR, F_Zx, (0.5*dt/dz), -1, nx+1, -1, ny+1, -1, nz+1, 2);
        DRAGONWING::releaseArray(&F_Zx);
    //_zLR needs (-1...nx, -1...ny, -1...nz)
    correctState(_zL, _zR, _zL,_zR, F_Xy, (0.5*dt/dx), -1, nx+1, -1, ny+1, -1, nz+1, 0);
        DRAGONWING::releaseArray(&F_Xy);
    correctState(_zL, _zR, _zL,_zR, F_Yx, (0.5*dt/dy), -1, nx+1, -1, ny+1, -1, nz+1, 1);
        DRAGONWING::releaseArray(&F_Yx);
    
    #ifdef MHD//Restore face-normal fields
    CT::copyFaceFields_X(_xL, B, _xR);
    CT::copyFaceFields_Y(_yL, B, _yR);
    CT::copyFaceFields_Z(_zL, B, _zR);
    #endif
#endif

    //Compute Fluxes
    computeFlux_X(_xL, _xR, F_X, 0, nx, -1, ny+1, -1, nz+1, dt/dx); //F_X needs (0...nx, -1...ny, -1...nz)
    computeFlux_Y(_yL, _yR, F_Y, -1, nx+1, 0, ny, -1, nz+1, dt/dy); //F_Y needs (-1...nx, 0...ny, -1...nz)
    computeFlux_Z(_zL, _zR, F_Z, -1, nx+1, -1, ny+1, 0, nz, dt/dz); //F_Z needs (-1...nx, -1...ny, 0...nz)
    
    FluidArray3D& _w = _xL; //Repurpose grid that isn't being used anymore
        DRAGONWING::releaseArray(&_xR);//Release the others
        DRAGONWING::releaseArray(&_yL);
        DRAGONWING::releaseArray(&_yR);
        DRAGONWING::releaseArray(&_zL);
        DRAGONWING::releaseArray(&_zR);
    
    
    //Preliminarily apply all fluxes
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                ConservativeState U(w[i,j,k]);
                U += (dt/dx) * (F_X[i,j,k] - F_X[i+1,j,k]);
                U += (dt/dy) * (F_Y[i,j,k] - F_Y[i,j+1,k]);
                U += (dt/dz) * (F_Z[i,j,k] - F_Z[i,j,k+1]);
                _w[i,j,k] = U;

                if(!U.isFinite()){
                        DRAGONWING::releaseArray(&F_X);
                        DRAGONWING::releaseArray(&F_Y);
                        DRAGONWING::releaseArray(&F_Z);
                        DRAGONWING::releaseArray(&_w);
                        #ifdef MHD
                        DRAGONWING::releaseArray(&B);
                        #endif
                    throw std::runtime_error(std::format("\tNaN state would be produced at ({},{},{})\n",i,j,k));
                }
            }
        }
    }
    //Preliminary CT Update
    #ifdef MHD
    MagneticArray3D& _A = *DRAGONWING::requestVec3Array(nx+1, ny+1, nz+1, ghosts);
    _A.clone(A);
    CT::updatePotential(_A, F_X, F_Y, F_Z, dt);
    CT::computeFaceFields(_A, B, dx, dy, dz);
    CT::computeBodyFields(B, _w);
        DRAGONWING::releaseArray(&B);
    #endif
        DRAGONWING::releaseArray(&F_X);
        DRAGONWING::releaseArray(&F_Y);
        DRAGONWING::releaseArray(&F_Z);
    
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                if(!_w[i,j,k].isPhysical()){
                        DRAGONWING::releaseArray(&_w);
                        #ifdef MHD
                        DRAGONWING::releaseArray(&_A);
                        #endif
                    throw std::runtime_error(std::format("Unphysical state would be produced by CT at ({},{},{})",i,j,k));
                }
                   
            }
        }
    }
    
    //Wait for any parallel grids to finish
    DRAGONWING::reportCheckpoint1();
    if(!DRAGONWING::waitForCheckpoint1()) return;
    
    //Commit Flux updates
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                w[i,j,k] = _w[i,j,k];
            }
        }
    }
        DRAGONWING::releaseArray(&_w);
    #ifdef MHD
    A.clone(_A);
        DRAGONWING::releaseArray(&_A);
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

void correctState(FluidArray2D& _L, FluidArray2D& _R, const FluxArray2D& F, double dt_dL, int xL, int xR, int yL, int yR, int dim){
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
void correctState(const FluidArray3D& _L0, const FluidArray3D& _R0, FluidArray3D& _L, FluidArray3D& _R, const FluxArray3D& F, double dt_dL, int xL, int xR, int yL, int yR, int zL, int zR, int dim){
    //Extract direction encoding
    int isX = dim%3 == 0 ? 1 : 0;
    int isY = dim%3 == 1 ? 1 : 0;
    int isZ = dim%3 == 2 ? 1 : 0;
    //cycle
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<yR; j++){
            for(int k=zL; k<zR; k++){
                auto trans = (F[i+isX, j+isY, k+isZ] - F[i,j,k]) * dt_dL;
                _L[i,j,k] = _L0[i,j,k] -  trans;
                _R[i,j,k] = _R0[i,j,k] -  trans;
            }
        }
    }
}
