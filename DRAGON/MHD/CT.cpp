
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

#ifdef MHD

//MARK: CT Updates

void CT::updatePotential(MagneticArray2D& _A, const FluxArray2D& F_X,const FluxArray2D& F_Y, double dt){
    const int nx = _A.getSizeX()-1, ny = _A.getSizeY()-1;
    
    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++){
            _A[i,j].z -= (F_Y[i-1,j].B.x -  F_X[i,j-1].B.y + F_Y[i,j].B.x  - F_X[i,j].B.y) * 0.25 * dt;
        }
    }
}
void CT::updatePotential(MagneticArray3D& _A, const FluxArray3D& F_X,const FluxArray3D& F_Y, const FluxArray3D& F_Z,  double dt){
    const int nx = _A.getSizeX()-1, ny = _A.getSizeY()-1, nz = _A.getSizeZ()-1;
    
    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++){
            for(int k=0; k<=nz; k++){
                _A[i,j,k].x -= (F_Z[i,j-1,k].B.y -  F_Y[i,j,k-1].B.z + F_Z[i,j,k].B.y  - F_Y[i,j,k].B.z) * 0.25 * dt;
                _A[i,j,k].y -= (F_X[i,j,k-1].B.z -  F_Z[i-1,j,k].B.x + F_X[i,j,k].B.z  - F_Z[i,j,k].B.x) * 0.25* dt;
                _A[i,j,k].z -= (F_Y[i-1,j,k].B.x -  F_X[i,j-1,k].B.y + F_Y[i,j,k].B.x  - F_X[i,j,k].B.y) * 0.25 * dt;
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
            _B[i,j].x = (_A[i,j+1].z - _A[i,j].z)*_dy;
        }
    }
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            _B[i,j].y = (_A[i,j].z - _A[i+1,j].z)*_dx;
        }
    }
}
void CT::computeFaceFields(const MagneticArray3D& _A, MagneticArray3D& _B, double dx, double dy, double dz){
    const int nx = _A.getSizeX()-1, ny = _A.getSizeY()-1, nz = _A.getSizeZ()-1, g = _B.getGhosts();
    const double _dx = 1/dx, _dy = 1/dy, _dz = 1/dz; //One division + n^3 multiplications > n^3 divisions
    
    for(int i=-g; i<=nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g; k<nz+g; k++){
                _B[i,j,k].x = (_A[i,j+1,k].z - _A[i,j,k].z)*_dy - (_A[i,j,k+1].y - _A[i,j,k].y)*_dz;
            }
        }
    }
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<=ny+g; j++){
            for(int k=-g; k<nz+g; k++){
                _B[i,j,k].y = (_A[i,j,k+1].x - _A[i,j,k].x)*_dz - (_A[i+1,j,k].z - _A[i,j,k].z)*_dx;
            }
        }
    }
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g; k<=nz+g; k++){
                _B[i,j,k].z = (_A[i+1,j,k].y - _A[i,j,k].y)*_dx -  (_A[i,j+1,k].x - _A[i,j,k].x)*_dy;
            }
        }
    }
}

//MARK: Body Fields
void Grid2D::initialize_B_fields(){
    const int nx = w.getSizeX(), ny = w.getSizeY();
    MagneticArray2D B(nx+1,ny+1,w.getGhosts());
    CT::computeFaceFields(A, B, dx, dy);
    
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            w[i,j].B.x = (B[i,j].x + B[i+1,j].x)/2;
            w[i,j].B.y = (B[i,j].y + B[i,j+1].y)/2;
        }
    }
}
void Grid2D::computeBodyAveragedFields(const MagneticArray2D& B){
    const int nx = w.getSizeX(), ny = w.getSizeY();

    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            ConservativeState U(w[i,j]);//Update the conservative element to keep quantities conserved
            U.B.x = (B[i,j].x + B[i+1,j].x)/2;
            U.B.y = (B[i,j].y + B[i,j+1].y)/2;
            w[i,j] = U;
        }
    }
}
void Grid3D::initialize_B_fields(){
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ();
    MagneticArray3D B(nx+1,ny+1,nz+1,w.getGhosts());
    CT::computeFaceFields(A, B, dx, dy, dz);
    
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                w[i,j,k].B.x = (B[i,j,k].x + B[i+1,j,k].x)/2;
                w[i,j,k].B.y = (B[i,j,k].y + B[i,j+1,k].y)/2;
                w[i,j,k].B.z = (B[i,j,k].z + B[i,j,k+1].z)/2;
            }
        }
    }

}
void Grid3D::computeBodyAveragedFields(const MagneticArray3D& B){
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ();

    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
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
#endif
