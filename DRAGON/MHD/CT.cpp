
//
//  CT.cpp
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

constexpr double CTU_TOL = 1e-18;
inline double upwindCorr(double massFlux, double faceEminus, double faceEplus, double refMinus, double refPlus);

#ifdef MHD

//MARK: Electric Field Calculation
void CT::computeElectric(MagneticArray2D& E, const FluxArray2D& F_X,const FluxArray2D& F_Y, int g){
    const int nx = E.getSizeX()-1, ny = E.getSizeY()-1;
    
    for(int i=-g; i<=nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            E[i,j].z = (F_Y[i-1,j].B.x -  F_X[i,j-1].B.y + F_Y[i,j].B.x  - F_X[i,j].B.y) * 0.25;
        }
    }
}


void CT::computeElectric(MagneticArray3D& E, const FluxArray3D& F_X,const FluxArray3D& F_Y, const FluxArray3D& F_Z, const FluidArray3D& wref, int g){
    const int nx = E.getSizeX()-1, ny = E.getSizeY()-1, nz = E.getSizeZ()-1;
    
    //Base Terms
    for(int i=-g; i<=nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            for(int k=-g; k<=nz+g; k++){
                E[i,j,k].x = (F_Z[i,j-1,k].B.y -  F_Y[i,j,k-1].B.z + F_Z[i,j,k].B.y  - F_Y[i,j,k].B.z) * 0.25;
                E[i,j,k].y = (F_X[i,j,k-1].B.z -  F_Z[i-1,j,k].B.x + F_X[i,j,k].B.z  - F_Z[i,j,k].B.x) * 0.25;
                E[i,j,k].z = (F_Y[i-1,j,k].B.x -  F_X[i,j-1,k].B.y + F_Y[i,j,k].B.x  - F_X[i,j,k].B.y) * 0.25;
            }
        }
    }
    
    #ifdef CTU //Upwind terms
    //Calculate Body-centred E
        auto __E = DRAGONWING::requestVec3Arrays(1, nx, ny, nz, E.getGhosts());
    MagneticArray3D& Eref = *__E[0];
    for(int i=-g-1; i <= nx+g; i++){
        for(int j=-g-1; j <= ny+g; j++){
            for(int k=-g-1; k <= nz+g; k++){
                Eref[i,j,k] = cross(wref[i,j,k].B, wref[i,j,k].v);
            }
        }
    }
    
    for(int i=-g; i<=nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            for(int k=-g; k<=nz+g; k++){
                // X upwinds
                E[i,j,k].y += upwindCorr(F_X[i,j,k].rho, -F_Z[i-1,j,k].B.x, -F_Z[i,j,k].B.x, Eref[i-1,j,k].y, Eref[i,j,k].y);
                E[i,j,k].y -= upwindCorr(F_X[i,j,k-1].rho, -F_Z[i-1,j,k].B.x, -F_Z[i,j,k].B.x, Eref[i-1,j,k-1].y, Eref[i,j,k-1].y);
                E[i,j,k].z += upwindCorr(F_X[i,j,k].rho, F_Y[i-1,j,k].B.x, F_Y[i,j,k].B.x, Eref[i-1,j,k].z, Eref[i,j,k].z);
                E[i,j,k].z -= upwindCorr(F_X[i,j-1,k].rho, F_Y[i-1,j,k].B.x, F_Y[i,j,k].B.x, Eref[i-1,j-1,k].z, Eref[i,j-1,k].z);
                // Y upwinds
                E[i,j,k].x += upwindCorr(F_Y[i,j,k].rho, F_Z[i,j-1,k].B.y, F_Z[i,j,k].B.y, Eref[i,j-1,k].x, Eref[i,j,k].x);
                E[i,j,k].x -= upwindCorr(F_Y[i,j,k-1].rho, F_Z[i,j-1,k].B.y, F_Z[i,j,k].B.y, Eref[i,j-1,k-1].x, Eref[i,j,k-1].x);
                E[i,j,k].z += upwindCorr(F_Y[i,j,k].rho, -F_X[i,j-1,k].B.y, -F_X[i,j,k].B.y, Eref[i,j-1,k].z, Eref[i,j,k].z);
                E[i,j,k].z -= upwindCorr(F_Y[i-1,j,k].rho, -F_X[i,j-1,k].B.y, -F_X[i,j,k].B.y, Eref[i-1,j-1,k].z, Eref[i-1,j,k].z);
                // Z upwinds
                E[i,j,k].x += upwindCorr(F_Z[i,j,k].rho, -F_Y[i,j,k-1].B.z, -F_Y[i,j,k].B.z, Eref[i,j,k-1].x, Eref[i,j,k].x);
                E[i,j,k].x -= upwindCorr(F_Z[i,j-1,k].rho, -F_Y[i,j,k-1].B.z, -F_Y[i,j,k].B.z, Eref[i,j-1,k-1].x, Eref[i,j-1,k].x);
                E[i,j,k].y += upwindCorr(F_Z[i,j,k].rho, F_X[i,j,k-1].B.z, F_X[i,j,k].B.z, Eref[i,j,k-1].y, Eref[i,j,k].y);
                E[i,j,k].y -= upwindCorr(F_Z[i-1,j,k].rho, F_X[i,j,k-1].B.z, F_X[i,j,k].B.z, Eref[i-1,j,k-1].y, Eref[i-1,j,k].y);
            }
        }
    }
    #endif
}




//MARK: CTU Corrections

inline double upwindCorr(double massFlux, double faceEminus, double faceEplus, double refMinus, double refPlus) {
    if (massFlux > CTU_TOL) return 0.25 * (faceEminus - refMinus);
    if (massFlux < -CTU_TOL) return 0.25 * (faceEplus - refPlus);
    return 0.125 * (faceEminus + faceEplus - refMinus - refPlus);
}






//MARK: Vector Potential Update
void CT::updatePotential(MagneticArray2D& _A, const MagneticArray2D& E, double dt, int g){
    const int nx = _A.getSizeX()-1, ny = _A.getSizeY()-1;
    
    for(int i=-g; i<=nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            _A[i,j].z -= E[i,j].z * dt;
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
//MARK: Face Fields
void CT::computeFaceFields(const MagneticArray2D& _A, MagneticArray2D& _B, double dx, double dy){
    const int nx = _A.getSizeX()-1, ny = _A.getSizeY()-1, g = _B.getGhosts();
    const double _dx = 1/dx, _dy = 1/dy; //One division + n^2 multiplications > n^2 divisions
    
    for(int i=-g; i<=nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            _B[i,j].x = (_A[i,j+1].z - _A[i,j].z) * _dy;
        }
    }
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            _B[i,j].y = (_A[i,j].z - _A[i+1,j].z) * _dx;
        }
    }
}
void CT::computeFaceFields(const MagneticArray3D& _A, MagneticArray3D& _B, double dx, double dy, double dz){
    const int nx = _A.getSizeX()-1, ny = _A.getSizeY()-1, nz = _A.getSizeZ()-1, g = _B.getGhosts();
    const double _dx = 1/dx, _dy = 1/dy, _dz = 1/dz; //One division + n^3 multiplications > n^3 divisions
    
    for(int i=-g; i<=nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g; k<nz+g; k++){
                _B[i,j,k].x = (_A[i,j+1,k].z - _A[i,j,k].z)*_dy - (_A[i,j,k+1].y - _A[i,j,k].y) * _dz;
            }
        }
    }
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            for(int k=-g; k<nz+g; k++){
                _B[i,j,k].y = (_A[i,j,k+1].x - _A[i,j,k].x)*_dz - (_A[i+1,j,k].z - _A[i,j,k].z) * _dx;
            }
        }
    }
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g; k<=nz+g; k++){
                _B[i,j,k].z = (_A[i+1,j,k].y - _A[i,j,k].y)*_dx -  (_A[i,j+1,k].x - _A[i,j,k].x) * _dy;
            }
        }
    }
}

//MARK: Body Fields
void Grid2D::initialize_B_fields(){
    const int nx = w.getSizeX(), ny = w.getSizeY(), ng = w.getGhosts();
    boundary.apply(*this);

    MagneticArray2D B(nx+1,ny+1,w.getGhosts());
    CT::computeFaceFields(A, B, dx, dy);
    
    for(int i=-ng; i<nx+ng; i++){
        for(int j=-ng; j<ny+ng; j++){
            w[i,j].B.x = (B[i,j].x + B[i+1,j].x)/2;
            w[i,j].B.y = (B[i,j].y + B[i,j+1].y)/2;
        }
    }
}
void CT::computeBodyFields(const MagneticArray2D& B, FluidArray2D& w){
    const int nx = w.getSizeX(), ny = w.getSizeY(), ng = w.getGhosts();

    for(int i=-ng; i<nx+ng; i++){
        for(int j=-ng; j<ny+ng; j++){
            ConservativeState U(w[i,j]);//Update the conservative element to keep quantities conserved
            U.B.x = (B[i,j].x + B[i+1,j].x)/2;
            U.B.y = (B[i,j].y + B[i,j+1].y)/2;
            w[i,j] = U;
        }
    }
}
void Grid2D::computeBodyAveragedFields(const MagneticArray2D& B){ CT::computeBodyFields(B, w); }

void Grid3D::initialize_B_fields(){
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ng = w.getGhosts();
    boundary.apply(*this);

    MagneticArray3D B(nx+1,ny+1,nz+1,w.getGhosts());
    CT::computeFaceFields(A, B, dx, dy, dz);
    
    for(int i=-ng; i<nx+ng; i++){
        for(int j=-ng; j<ny+ng; j++){
            for(int k=-ng; k<nz+ng; k++){
                w[i,j,k].B.x = (B[i,j,k].x + B[i+1,j,k].x)/2;
                w[i,j,k].B.y = (B[i,j,k].y + B[i,j+1,k].y)/2;
                w[i,j,k].B.z = (B[i,j,k].z + B[i,j,k+1].z)/2;
            }
        }
    }

}

void Grid3D::computeBodyAveragedFields(const MagneticArray3D& B){ CT::computeBodyFields(B, w); }
void CT::computeBodyFields(const MagneticArray3D& B, FluidArray3D& w){
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ng = w.getGhosts();

        for(int i=-ng; i<nx+ng; i++){
            for(int j=-ng; j<ny+ng; j++){
                for(int k=-ng; k<nz+ng; k++){
                    ConservativeState U(w[i,j,k]);//Update the conservative element to keep quantities conserved
                    U.B.x = (B[i,j,k].x + B[i+1,j,k].x)/2;
                    U.B.y = (B[i,j,k].y + B[i,j+1,k].y)/2;
                    U.B.z = (B[i,j,k].z + B[i,j,k+1].z)/2;
                    w[i,j,k] = U;
                }
            }
        }
}



//MARK: Face-Centred Field Copy
void CT::copyFaceFields_X( FluidArray2D& _L,const MagneticArray2D& _B, FluidArray2D& _R){
    const int nx = _L.getSizeX(), ny = _L.getSizeY(), g = _B.getGhosts();
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            _L[i,j].B.x = _B[i,j].x;
            _R[i,j].B.x = _B[i+1,j].x;
        }
    }
}
void CT::copyFaceFields_Y( FluidArray2D& _L,const MagneticArray2D& _B, FluidArray2D& _R){
    const int nx = _L.getSizeX(), ny = _L.getSizeY(), g = _B.getGhosts();
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            _L[i,j].B.y = _B[i,j].y;
            _R[i,j].B.y = _B[i,j+1].y;
        }
    }
}
void CT::copyFaceFields_X( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R){
    const int nx = _L.getSizeX(), ny = _L.getSizeY(), nz = _L.getSizeZ(), g = _B.getGhosts();
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g; k<nz+g; k++){
                _L[i,j,k].B.x = _B[i,j,k].x;
                _R[i,j,k].B.x = _B[i+1,j,k].x;
            }
        }
    }
}
void CT::copyFaceFields_Y( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R){
    const int nx = _L.getSizeX(), ny = _L.getSizeY(), nz = _L.getSizeZ(), g = _B.getGhosts();
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g; k<nz+g; k++){
                _L[i,j,k].B.y = _B[i,j,k].y;
                _R[i,j,k].B.y = _B[i,j+1,k].y;
            }
        }
    }
}
void CT::copyFaceFields_Z( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R){
    const int nx = _L.getSizeX(), ny = _L.getSizeY(), nz = _L.getSizeZ(), g = _B.getGhosts();
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g; k<nz+g; k++){
                _L[i,j,k].B.z = _B[i,j,k].z;
                _R[i,j,k].B.z = _B[i,j,k+1].z;
            }
        }
    }
}

#ifndef DIMENSION_UNSPLIT
#error MHD requires Unsplit dimensions. Please enable DIMENSION_UNSPLIT in Config.h
#endif
#endif


