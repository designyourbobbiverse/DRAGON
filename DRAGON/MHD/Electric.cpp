
//
//  Electric.cpp
//  DRAGON/MHD
//
//  Created by Bobbie Markwick on 24/06/2026.
//  CT Implementaiton based in part on https://doi.org/10.1088/0067-0049/182/2/519
//      Evans and Hawley (1988). https://doi.org/10.1086/166684
//      Gardiner and Stone (2005). https://doi.org/10.1016/j.jcp.2004.11.016


#include "CT.hpp"
#include "Config.h"
#include "DragonWing.hpp"
#include <cmath>


#ifdef MHD
#ifndef DIMENSION_UNSPLIT
#error MHD requires Unsplit dimensions. Please enable DIMENSION_UNSPLIT in Config.h
#endif



//MARK: Vector Potential Update
void CT::updatePotential(MagneticArray2D& _A, const MagneticArray2D& E, double dt, int g){
    const int nx = _A.getSizeX()-1, ny = _A.getSizeY()-1;
    
    for(int i=-g; i<=nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            _A[i,j] -= E[i,j] * dt;
        }
    }
}
void CT::updatePotential(MagneticArray3D& _A, const MagneticArray3D& E, double dt, int g){
    const int nx = _A.getSizeX()-1, ny = _A.getSizeY()-1, nz = _A.getSizeZ()-1;
    
    for(int i=-g; i<=nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            for(int k=-g; k<=nz+g; k++){
                _A[i,j,k] -= E[i,j,k] * dt;
            }
        }
    }
}



//MARK: E field from Fluxes
void CT::computeElectric(MagneticArray2D& E, const FluxArray2D& F_X,const FluxArray2D& F_Y,  const int g){
    const int nx = E.getSizeX()-1, ny = E.getSizeY()-1;
    for(int i=-g; i<=nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            E[i,j].x = - F_Y[i,j].B.z;
            E[i,j].y = F_X[i,j].B.z;
            E[i,j].z = (F_Y[i-1,j].B.x -  F_X[i,j-1].B.y + F_Y[i,j].B.x  - F_X[i,j].B.y) * 0.25;
        }
    }
}


void CT::computeElectric(MagneticArray3D& E, const FluxArray3D& F_X,const FluxArray3D& F_Y, const FluxArray3D& F_Z, const int g){
    const int nx = E.getSizeX()-1, ny = E.getSizeY()-1, nz = E.getSizeZ()-1;
    for(int i=-g; i<=nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            for(int k=-g; k<=nz+g; k++){
                E[i,j,k].x = (F_Z[i,j-1,k].B.y -  F_Y[i,j,k-1].B.z + F_Z[i,j,k].B.y  - F_Y[i,j,k].B.z) * 0.25;
                E[i,j,k].y = (F_X[i,j,k-1].B.z -  F_Z[i-1,j,k].B.x + F_X[i,j,k].B.z  - F_Z[i,j,k].B.x) * 0.25;
                E[i,j,k].z = (F_Y[i-1,j,k].B.x -  F_X[i,j-1,k].B.y + F_Y[i,j,k].B.x  - F_X[i,j,k].B.y) * 0.25;
            }
        }
    }
}

    
//MARK: CTU Upwinding

constexpr double CTU_TOL = 1e-18;
inline double upwindCorr(double massFlux, double faceEminus, double faceEplus, double refMinus, double refPlus) {
    if (massFlux > CTU_TOL) return 0.25 * (faceEminus - refMinus);
    if (massFlux < -CTU_TOL) return 0.25 * (faceEplus - refPlus);
    return 0.125 * (faceEminus + faceEplus - refMinus - refPlus);
}

void CT::upwindElectric(MagneticArray2D& E, const FluxArray2D& F_X,const FluxArray2D& F_Y,const MagneticArray2D& Eref, int g){
    const int nx = E.getSizeX()-1, ny = E.getSizeY()-1;
    for(int i=-g; i<=nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            // X upwinds
            E[i,j].z += upwindCorr(F_X[i,j].rho, F_Y[i-1,j].B.x, F_Y[i,j].B.x, Eref[i-1,j].z, Eref[i,j].z);
            E[i,j].z += upwindCorr(F_X[i,j-1].rho, F_Y[i-1,j].B.x, F_Y[i,j].B.x, Eref[i-1,j-1].z, Eref[i,j-1].z);
            // Y upwinds
            E[i,j].z += upwindCorr(F_Y[i,j].rho, -F_X[i,j-1].B.y, -F_X[i,j].B.y, Eref[i,j-1].z, Eref[i,j].z);
            E[i,j].z += upwindCorr(F_Y[i-1,j].rho, -F_X[i,j-1].B.y, -F_X[i,j].B.y, Eref[i-1,j-1].z, Eref[i-1,j].z);
            // Z upwinds
            E[i,j].x += 0.25 * (-2*F_Y[i,j].B.z - Eref[i,j].x - Eref[i,j-1].x);
            E[i,j].y += 0.25 * (2*F_X[i,j].B.z - Eref[i,j].y - Eref[i-1,j].y);
        }
    }
}

    
void CT::upwindElectric(MagneticArray3D& E, const FluxArray3D& F_X,const FluxArray3D& F_Y, const FluxArray3D& F_Z, const MagneticArray3D& Eref, int g){
    const int nx = E.getSizeX()-1, ny = E.getSizeY()-1, nz = E.getSizeZ()-1;
    for(int i=-g; i<=nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            for(int k=-g; k<=nz+g; k++){
                // X upwinds
                E[i,j,k].y += upwindCorr(F_X[i,j,k].rho, -F_Z[i-1,j,k].B.x, -F_Z[i,j,k].B.x, Eref[i-1,j,k].y, Eref[i,j,k].y);
                E[i,j,k].y += upwindCorr(F_X[i,j,k-1].rho, -F_Z[i-1,j,k].B.x, -F_Z[i,j,k].B.x, Eref[i-1,j,k-1].y, Eref[i,j,k-1].y);
                E[i,j,k].z += upwindCorr(F_X[i,j,k].rho, F_Y[i-1,j,k].B.x, F_Y[i,j,k].B.x, Eref[i-1,j,k].z, Eref[i,j,k].z);
                E[i,j,k].z += upwindCorr(F_X[i,j-1,k].rho, F_Y[i-1,j,k].B.x, F_Y[i,j,k].B.x, Eref[i-1,j-1,k].z, Eref[i,j-1,k].z);
                // Y upwinds
                E[i,j,k].x += upwindCorr(F_Y[i,j,k].rho, F_Z[i,j-1,k].B.y, F_Z[i,j,k].B.y, Eref[i,j-1,k].x, Eref[i,j,k].x);
                E[i,j,k].x += upwindCorr(F_Y[i,j,k-1].rho, F_Z[i,j-1,k].B.y, F_Z[i,j,k].B.y, Eref[i,j-1,k-1].x, Eref[i,j,k-1].x);
                E[i,j,k].z += upwindCorr(F_Y[i,j,k].rho, -F_X[i,j-1,k].B.y, -F_X[i,j,k].B.y, Eref[i,j-1,k].z, Eref[i,j,k].z);
                E[i,j,k].z += upwindCorr(F_Y[i-1,j,k].rho, -F_X[i,j-1,k].B.y, -F_X[i,j,k].B.y, Eref[i-1,j-1,k].z, Eref[i-1,j,k].z);
                // Z upwinds
                E[i,j,k].x += upwindCorr(F_Z[i,j,k].rho, -F_Y[i,j,k-1].B.z, -F_Y[i,j,k].B.z, Eref[i,j,k-1].x, Eref[i,j,k].x);
                E[i,j,k].x += upwindCorr(F_Z[i,j-1,k].rho, -F_Y[i,j,k-1].B.z, -F_Y[i,j,k].B.z, Eref[i,j-1,k-1].x, Eref[i,j-1,k].x);
                E[i,j,k].y += upwindCorr(F_Z[i,j,k].rho, F_X[i,j,k-1].B.z, F_X[i,j,k].B.z, Eref[i,j,k-1].y, Eref[i,j,k].y);
                E[i,j,k].y += upwindCorr(F_Z[i-1,j,k].rho, F_X[i,j,k-1].B.z, F_X[i,j,k].B.z, Eref[i-1,j,k-1].y, Eref[i-1,j,k].y);
            }
        }
    }
}






//MARK: Body Electric Fields
void CT::bodyElectric(const FluidArray2D &w, MagneticArray2D &E, int g){
    const int nx = w.getSizeX(), ny = w.getSizeY();
    for(int i=-g-1; i <= nx+g; i++){
        for(int j=-g-1; j <= ny+g; j++){
            E[i,j] = cross(w[i,j].B, w[i,j].v);
        }
    }
}
void CT::bodyElectric(const FluidArray3D &w, MagneticArray3D &E, int g){
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ();
    for(int i=-g-1; i <= nx+g; i++){
        for(int j=-g-1; j <= ny+g; j++){
            for(int k=-g-1; k <= nz+g; k++){
                E[i,j,k] = cross(w[i,j,k].B, w[i,j,k].v);
            }
        }
    }
}

#endif


