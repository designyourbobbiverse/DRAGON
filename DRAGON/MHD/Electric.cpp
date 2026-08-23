
//
//  Electric.cpp
//  DRAGON/MHD
//
//  Created by Bobbie Markwick on 24/06/2026.
//  CT Implementaiton based in part on https://arxiv.org/abs/0710.0424
//      Evans and Hawley (1988). https://doi.org/10.1086/166684
//      Gardiner and Stone (2008). https://arxiv.org/abs/0712.2634


#include "MHD/CT.hpp"

#include "Config.h"
using namespace DRAGON;

#ifdef MHD
#ifndef DIMENSION_UNSPLIT
#error MHD requires Unsplit dimensions. Please enable DIMENSION_UNSPLIT in Config.h
#endif



//MARK: E field from Fluxes
//Ex is the By component of the flux through z-faces. Same for any cyclic permuation of xyz, or negative for anticyclic.
//For the elctric field on an edge, average the fluxes through the four adjacent faces

void CT::computeElectric(MagneticArray2D& E, const FluxArray2D& F_X,const FluxArray2D& F_Y,  const int g){
    const int nx = E.getSizeX()-1, ny = E.getSizeY()-1;
    for(int i=-g; i<=nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            E[i,j].x = -F_Y[i,j].B.z;
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
//See section 4 of Gardiner and Stone (2008). https://arxiv.org/abs/0712.2634
//Each edge-centered electric field component gets upwinded four times: once for each of the four connected faces
//      The full upwinding is an average of the partial-upwindings of these four faces, hence the x0.25
//Each of these partial-upwinding calculations needs (upwindCorr receives the possibilities for 2 and 3, then decides which one to use)
//      1) the velocity (or mass flux) through the corresponding face
//      2) the body-centered E-field cell upwind of that face, and
//      3) the E-field (extracted from flux, see above) on the other face that touches both the above body-cell and this edge


constexpr double CTU_TOL = 1e-18;
inline double upwindCorr(double massFlux, double faceEminus, double faceEplus, double bodyEminus, double bodyEplus) {
    //massFlux should have the same sign as velocity, so we can use it to detect the flow direction
    if (massFlux > CTU_TOL) return 0.25 * (faceEminus - bodyEminus); //Positive velocity -> upwind from the left
    if (massFlux < -CTU_TOL) return 0.25 * (faceEplus - bodyEplus); //Negative velocity -> upwind from the right
    return 0.125 * (faceEminus + faceEplus - bodyEminus - bodyEplus); //If close enough to zero, average from both directions.
}

void CT::upwindElectric(MagneticArray2D& E, const FluxArray2D& F_X,const FluxArray2D& F_Y,const MagneticArray2D& Ebody, int g){
    const int nx = E.getSizeX()-1, ny = E.getSizeY()-1;
    for(int i=-g; i<=nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            // X upwinds
            E[i,j].z += upwindCorr(F_X[i,j].rho, F_Y[i-1,j].B.x, F_Y[i,j].B.x, Ebody[i-1,j].z, Ebody[i,j].z);
            E[i,j].z += upwindCorr(F_X[i,j-1].rho, F_Y[i-1,j].B.x, F_Y[i,j].B.x, Ebody[i-1,j-1].z, Ebody[i,j-1].z);
            // Y upwinds
            E[i,j].z += upwindCorr(F_Y[i,j].rho, -F_X[i,j-1].B.y, -F_X[i,j].B.y, Ebody[i,j-1].z, Ebody[i,j].z);
            E[i,j].z += upwindCorr(F_Y[i-1,j].rho, -F_X[i,j-1].B.y, -F_X[i,j].B.y, Ebody[i-1,j-1].z, Ebody[i-1,j].z);
            // Z upwinds
            E[i,j].x += 0.25 * (-2*F_Y[i,j].B.z - Ebody[i,j].x - Ebody[i,j-1].x);
            E[i,j].y += 0.25 * (2*F_X[i,j].B.z - Ebody[i,j].y - Ebody[i-1,j].y);
        }
    }
}

    
void CT::upwindElectric(MagneticArray3D& E, const FluxArray3D& F_X,const FluxArray3D& F_Y, const FluxArray3D& F_Z, const MagneticArray3D& Ebody, int g){
    const int nx = E.getSizeX()-1, ny = E.getSizeY()-1, nz = E.getSizeZ()-1;
    for(int i=-g; i<=nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            for(int k=-g; k<=nz+g; k++){
                // X upwinds
                E[i,j,k].y += upwindCorr(F_X[i,j,k].rho, -F_Z[i-1,j,k].B.x, -F_Z[i,j,k].B.x, Ebody[i-1,j,k].y, Ebody[i,j,k].y);
                E[i,j,k].y += upwindCorr(F_X[i,j,k-1].rho, -F_Z[i-1,j,k].B.x, -F_Z[i,j,k].B.x, Ebody[i-1,j,k-1].y, Ebody[i,j,k-1].y);
                E[i,j,k].z += upwindCorr(F_X[i,j,k].rho, F_Y[i-1,j,k].B.x, F_Y[i,j,k].B.x, Ebody[i-1,j,k].z, Ebody[i,j,k].z);
                E[i,j,k].z += upwindCorr(F_X[i,j-1,k].rho, F_Y[i-1,j,k].B.x, F_Y[i,j,k].B.x, Ebody[i-1,j-1,k].z, Ebody[i,j-1,k].z);
                // Y upwinds
                E[i,j,k].x += upwindCorr(F_Y[i,j,k].rho, F_Z[i,j-1,k].B.y, F_Z[i,j,k].B.y, Ebody[i,j-1,k].x, Ebody[i,j,k].x);
                E[i,j,k].x += upwindCorr(F_Y[i,j,k-1].rho, F_Z[i,j-1,k].B.y, F_Z[i,j,k].B.y, Ebody[i,j-1,k-1].x, Ebody[i,j,k-1].x);
                E[i,j,k].z += upwindCorr(F_Y[i,j,k].rho, -F_X[i,j-1,k].B.y, -F_X[i,j,k].B.y, Ebody[i,j-1,k].z, Ebody[i,j,k].z);
                E[i,j,k].z += upwindCorr(F_Y[i-1,j,k].rho, -F_X[i,j-1,k].B.y, -F_X[i,j,k].B.y, Ebody[i-1,j-1,k].z, Ebody[i-1,j,k].z);
                // Z upwinds
                E[i,j,k].x += upwindCorr(F_Z[i,j,k].rho, -F_Y[i,j,k-1].B.z, -F_Y[i,j,k].B.z, Ebody[i,j,k-1].x, Ebody[i,j,k].x);
                E[i,j,k].x += upwindCorr(F_Z[i,j-1,k].rho, -F_Y[i,j,k-1].B.z, -F_Y[i,j,k].B.z, Ebody[i,j-1,k-1].x, Ebody[i,j-1,k].x);
                E[i,j,k].y += upwindCorr(F_Z[i,j,k].rho, F_X[i,j,k-1].B.z, F_X[i,j,k].B.z, Ebody[i,j,k-1].y, Ebody[i,j,k].y);
                E[i,j,k].y += upwindCorr(F_Z[i-1,j,k].rho, F_X[i,j,k-1].B.z, F_X[i,j,k].B.z, Ebody[i-1,j,k-1].y, Ebody[i-1,j,k].y);
            }
        }
    }
}






//MARK: Body Electric Fields
//E = - v x B = B x v
void CT::bodyElectric(const FluidArray2D &w, MagneticArray2D &E, int g){ //1001001 in distress (iykyk)
    const int nx = w.getSizeX(), ny = w.getSizeY();
    for(int i=-g-1; i <= nx+g; i++){
        for(int j=-g-1; j <= ny+g; j++){
            E[i,j] = cross(w[i,j].B, w[i,j].v);
        }
    }
}
void CT::bodyElectric(const FluidArray3D &w, MagneticArray3D &E, int g){ //1001001 in distress (iykyk)
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


