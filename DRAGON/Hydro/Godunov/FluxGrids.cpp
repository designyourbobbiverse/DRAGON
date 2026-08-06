//
//  FluxGrids.cpp
//  DRAGON/Hydro/Godunov
//
//  Created by Bobbie Markwick on 5/07/2026.


#include "Grid.hpp"
#include "Godunov.hpp"

#include "Config.h"
#include "Riemann.hpp" //For Riemann Solvers
#include "TVD.hpp"     //For MUSCL Reconstruction
#include "CT.hpp"      //For CT::copyFaceFields

#include <stdexcept>      //For error handling
#include <format>         //For error message formatting
using namespace DRAGON;



//MARK: Flux Application
void Godunov::applyFluxes(const FluidArray2D& w, FluidArray2D& _w, const FluxArray2D& F_X, const FluxArray2D& F_Y,  double dt_dx, double dt_dy,  int g){
    const int nx = w.getSizeX(), ny = w.getSizeY();
    _w.clone(w);
    
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
void Godunov::applyFluxes(const FluidArray3D& w, FluidArray3D& _w, const FluxArray3D& F_X, const FluxArray3D& F_Y, const FluxArray3D& F_Z, double dt_dx, double dt_dy, double dt_dz, int g){
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ();
    _w.clone(w);
    
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


//MARK: Riemann Solver loops
//Compute X fluxes between Right(_R) and Left (_L) half-states over the entire grid
void Godunov::computeFlux_X(const FluidArray2D& _L, const FluidArray2D& _R, FluxArray2D& F, int xL, int xR, int yL, int yR, double dt_dx ){
    for(int i=xL; i<=xR; i++){
        for(int j=yL; j<yR; j++){
            F[i,j] = Riemann(_R[i-1,j], _L[i,j]).flux_X(dt_dx);
        }
    }
}
//Compute Y fluxes between Right(_R) and Left (_L) half-states over the entire grid
void Godunov::computeFlux_Y(const FluidArray2D& _L, const FluidArray2D& _R, FluxArray2D& F, int xL, int xR, int yL, int yR, double dt_dy){
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<=yR; j++){
            F[i,j] = Riemann(_R[i,j-1], _L[i,j]).flux_Y(dt_dy);
        }
    }
}
//Compute X fluxes between Right(_R) and Left (_L) half-states over the entire grid
void Godunov::computeFlux_X(const FluidArray3D& _L, const FluidArray3D& _R, FluxArray3D& F, int xL, int xR, int yL, int yR, int zL, int zR, double dt_dx){
    for(int i=xL; i<=xR; i++){
        for(int j=yL; j<yR; j++){
            for(int k=zL; k<zR; k++){
                F[i,j,k] = Riemann(_R[i-1,j,k], _L[i,j,k]).flux_X(dt_dx);
            }
        }
    }
}
//Compute Y fluxes between Right(_R) and Left (_L) half-states over the entire grid
void Godunov::computeFlux_Y(const FluidArray3D& _L, const FluidArray3D& _R, FluxArray3D& F, int xL, int xR, int yL, int yR, int zL, int zR, double dt_dy){
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<=yR; j++){
            for(int k=zL; k<zR; k++){
                F[i,j,k] = Riemann(_R[i,j-1,k], _L[i,j,k]).flux_Y(dt_dy);
            }
        }
    }
}
//Compute Z fluxes between Right(_R) and Left (_L) half-states over the entire grid
void Godunov::computeFlux_Z(const FluidArray3D& _L, const FluidArray3D& _R, FluxArray3D& F, int xL, int xR, int yL, int yR, int zL, int zR, double dt_dz){
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<yR; j++){
            for(int k=zL; k<=zR; k++){
                F[i,j,k] = Riemann(_R[i,j,k-1], _L[i,j,k]).flux_Z(dt_dz);
            }
        }
    }
}

//MARK: MUSCL Loops
//Apply MUSCL over the entire grid
void Godunov::computeHalfStates_X(FluidArray2D& _L, const Grid2D& _W, FluidArray2D& _R, const MagneticArray2D& B,  double dt){
    const double dt_dx = dt/_W.dx, dt_dy = dt/_W.dy;//Compute once
    const int nx = _W.getSizeX(), ny = _W.getSizeY(), g = _W.getGhosts();
    //MUSCL Reconstruction
    for(int i=-g+1; i<nx+g - 1; i++){
        for(int j=-g; j<ny+g; j++){
            vec3 dB = {0,0,0};
            #ifdef MHD
            dB.x = (B[i+1,j].x - B[i,j].x) * dt_dx;
            dB.y = (B[i,j+1].y - B[i,j].y) * dt_dy;
            #endif
            auto wL =_W[i-1,j], wR = _W[i+1,j];
            TVD::MUSCL(wL, _L[i,j], _W[i,j], _R[i,j], wR, dt_dx, dB);
        }
    }
    //Leftmost and rightmost ghosts
    for(int j=-g; j<ny+g; j++){
        _L[-g,j] = _W[-g,j]; _R[-g,j] = _W[-g,j];
        _L[nx-1+g,j] = _W[nx-1+g,j]; _R[nx-1+g,j] = _W[nx-1+g,j];
    }
    #ifdef MHD
    CT::copyFaceFields_X(_L, B, _R);
    #endif
}
void Godunov::computeHalfStates_Y(FluidArray2D& _L, const Grid2D& _W, FluidArray2D& _R, const MagneticArray2D& B, double dt){
    const double dt_dx = dt/_W.dx, dt_dy = dt/_W.dy;//Compute once
    const int nx = _W.getSizeX(), ny = _W.getSizeY(), g = _W.getGhosts();

    //MUSCL Reconstruction
    for(int i=-g; i<nx+g; i++){
        for(int j=-g + 1; j<ny+g - 1; j++){
            vec3 dB = {0,0,0};
            #ifdef MHD //Calculate with swapped XY
            dB.y = (B[i+1,j].x - B[i,j].x) * dt_dx;
            dB.x = (B[i,j+1].y - B[i,j].y) * dt_dy;
            #endif
            //Swap XY inputs to MUSCL, then swap output back
            auto wL =_W[i,j-1].swappedXY(), wR = _W[i,j+1].swappedXY();
            TVD::MUSCL(wL, _L[i,j], _W[i,j].swappedXY(), _R[i,j],wR, dt_dy, dB);
            _L[i,j].swapXY(); _R[i,j].swapXY();
        }
    }
    //Leftmost and rightmost ghosts
    for(int i=-g; i<nx+g; i++){
        _L[i,-g] = _W[i,-g]; _R[i,-g] = _W[i,-g];
        _L[i,ny-1+g] = _W[i,ny-1+g]; _R[i,ny-1+g] = _W[i,ny-1+g];
    }
    #ifdef MHD
    CT::copyFaceFields_Y(_L, B, _R);
    #endif
}
void Godunov::computeHalfStates_X(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, const MagneticArray3D& B, double dt){
    const double dt_dx = dt/_W.dx, dt_dy = dt/_W.dy, dt_dz = dt/_W.dz;//Compute once
    const int nx = _W.getSizeX(), ny = _W.getSizeY(),nz = _W.getSizeZ(), g = _W.getGhosts();

    //MUSCL Reconstruction
    for(int i=-g + 1; i<nx+g - 1; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g; k<nz+g; k++){
                vec3 dB = {0,0,0};
                #ifdef MHD
                dB.x = (B[i+1,j,k].x - B[i,j,k].x) * dt_dx;
                dB.y = (B[i,j+1,k].y - B[i,j,k].y) * dt_dy;
                dB.z = (B[i,j,k+1].z - B[i,j,k].z) * dt_dz;
                #endif
                auto wL = _W[i-1,j,k], wR = _W[i+1,j,k];
                TVD::MUSCL(wL, _L[i,j,k], _W[i,j,k], _R[i,j,k], wR, dt_dx, dB);
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
    #ifdef MHD
    CT::copyFaceFields_X(_L, B, _R);
    #endif
}
void Godunov::computeHalfStates_Y(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, const MagneticArray3D& B, double dt){
    const double dt_dx = dt/_W.dx, dt_dy = dt/_W.dy, dt_dz = dt/_W.dz;//Compute once
    const int nx = _W.getSizeX(), ny = _W.getSizeY(),nz = _W.getSizeZ(), g = _W.getGhosts();
    
    //MUSCL Reconstruction
    for(int i=-g; i<nx+g; i++){
        for(int j=-g + 1; j<ny+g - 1; j++){
            for(int k=-g; k<nz+g; k++){
                vec3 dB = {0,0,0};
                #ifdef MHD //Calculate with swapped XY
                dB.y = (B[i+1,j,k].x - B[i,j,k].x) * dt_dx;
                dB.x = (B[i,j+1,k].y - B[i,j,k].y) * dt_dy;
                dB.z = (B[i,j,k+1].z - B[i,j,k].z) * dt_dz;
                #endif
                //Swap XY inputs to MUSCL, then swap output back
                auto wL =_W[i,j-1,k].swappedXY(), wR = _W[i,j+1,k].swappedXY();
                TVD::MUSCL(wL, _L[i,j,k], _W[i,j,k].swappedXY(), _R[i,j,k],wR, dt_dy,dB);
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
    #ifdef MHD
    CT::copyFaceFields_Y(_L, B, _R);
    #endif
}
void Godunov::computeHalfStates_Z(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, const MagneticArray3D& B,  double dt){
    const double dt_dx = dt/_W.dx, dt_dy = dt/_W.dy, dt_dz = dt/_W.dz;//Compute once
    const int nx = _W.getSizeX(), ny = _W.getSizeY(),nz = _W.getSizeZ(), g = _W.getGhosts();
    
    //MUSCL Reconstruction
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g + 1; k<nz+g - 1; k++){
                vec3 dB = {0,0,0};
                #ifdef MHD //Calculate with swapped XZ
                dB.z = (B[i+1,j,k].x - B[i,j,k].x) * dt_dx;
                dB.y = (B[i,j+1,k].y - B[i,j,k].y) * dt_dy;
                dB.x = (B[i,j,k+1].z - B[i,j,k].z) * dt_dz;
                #endif
                //Swap XZ inputs to MUSCL, then swap output back
                auto wL = _W[i,j,k-1].swappedXZ(), wR = _W[i,j,k+1].swappedXZ();
                TVD::MUSCL(wL, _L[i,j,k], _W[i,j,k].swappedXZ(), _R[i,j,k],wR, dt_dz,dB);
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
    #ifdef MHD
    CT::copyFaceFields_Z(_L, B, _R);
    #endif
}
