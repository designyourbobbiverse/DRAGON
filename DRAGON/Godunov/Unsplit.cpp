//
//  Unsplit.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 10/06/2026.
//

#include "Grid.hpp"
#include "Riemann.hpp"
#include "Config.h"
#include "CFL.hpp"
#include "TVD.hpp"
#include <cassert>
#include <iostream>

//MARK: Split vs Unsplit
void Grid2D::advance(double dt){
#ifdef DimensionUnsplit
    advance_unsplit(dt);
#else
    advance_split(dt);
#endif
}

void Grid3D::advance(double dt){
#ifdef DimensionUnsplit
    advance_unsplit(dt);
#else
    advance_split(dt);
#endif
}

//MARK: 2D Flux Array

struct FluxGrid2D{
    
    FluxGrid2D(int nx_, int ny_, int ghosts=0): nx(nx_), ny(ny_), ghosts(ghosts) {
        f = new ConservativeState[(nx+2*ghosts)*(ny+2*ghosts)];
    }
    ~FluxGrid2D(){ delete[] f; }
    FluxGrid2D(const FluxGrid2D&) = delete; //No copying
    FluxGrid2D& operator=(const FluxGrid2D&) = delete;
    
    //Grid access
    ConservativeState& operator[](int i,int j){
    #ifdef TESTMODE
        assert(i + ghosts >= 0 && i < nx+ghosts);
        assert(j + ghosts >= 0 && j < ny+ghosts);
    #endif
        int m = (i+ghosts)*(ny+2*ghosts) + (j+ghosts);
        return f[m];
    }
    const ConservativeState& operator[](int i,int j) const{
    #ifdef TESTMODE
        assert(i + ghosts >= 0 && i < nx+ghosts);
        assert(j + ghosts >= 0 && j < ny+ghosts);
    #endif
        int m = (i+ghosts)*(ny+2*ghosts) + (j+ghosts);
        return f[m];
    }
    int getSizeX() const{ return nx; }
    int getSizeY() const{ return ny; }
    int getGhosts() const{ return ghosts; }
private:
    ConservativeState* f;
    int ghosts, nx, ny;
};

//MARK: 3D Flux Array

struct FluxGrid3D{
    
    FluxGrid3D(int nx_, int ny_, int nz_, int ghosts=0): nx(nx_), ny(ny_), nz(nz_), ghosts(ghosts) {
        f = new ConservativeState[(nx+2*ghosts)*(ny+2*ghosts)*(nz+2*ghosts)];
    }
    ~FluxGrid3D(){ delete[] f; }
    FluxGrid3D(const FluxGrid3D&) = delete; //No copying
    FluxGrid3D& operator=(const FluxGrid3D&) = delete;
    
    //Grid access
    ConservativeState& operator[](int i,int j, int k){
    #ifdef TESTMODE
        assert(i + ghosts >= 0 && i < nx+ghosts);
        assert(j + ghosts >= 0 && j < ny+ghosts);
        assert(k + ghosts >= 0 && k < nz+ghosts);
    #endif
        int m = ((i+ghosts)*(ny+2*ghosts) + (j+ghosts)) * (nz+2*ghosts) + (k+ghosts);
        return f[m];
    }
    const ConservativeState& operator[](int i,int j,int k) const{
    #ifdef TESTMODE
        assert(i + ghosts >= 0 && i < nx+ghosts);
        assert(j + ghosts >= 0 && j < ny+ghosts);
        assert(k + ghosts >= 0 && k < nz+ghosts);
    #endif
        int m = ((i+ghosts)*(ny+2*ghosts) + (j+ghosts)) * (nz+2*ghosts) + (k+ghosts);
        return f[m];
    }
    int getSizeX() const{ return nx; }
    int getSizeY() const{ return ny; }
    int getSizeZ() const{ return nz; }
    int getGhosts() const{ return ghosts; }
private:
    ConservativeState* f;
    int ghosts, nx, ny, nz;
};
    

//MARK: Unsplit CFL

void Grid2D::advance_unsplit(double dt){
    while(dt > Timestep_Tolerance){
        //CFL Time Constraint
        double t1 = std::min(dt, CFL::cfl_time(*this));
        dt -= t1;
        //Advance
        advanceXY(t1);
    }
}
void Grid3D::advance_unsplit(double dt){
    while(dt > Timestep_Tolerance){
        //CFL Time Constraint
        double t1 = std::min(dt, CFL::cfl_time(*this));
        dt -= t1;
        //Advance
        advanceXYZ(t1);
    }
}

//MARK: Unsplit TVD

void Grid2D::computeHalfStates(Grid2D& _L, Grid2D& _R, double dt, int dim){
    //Dimension
    int isX = 0, isY = 0;
    switch(dim % 2){
        case 0: isX = 1; break;
        case 1: isY = 1; break;
    }
    double dt_dL = dt/(isX*dx + isY*dy);
    //Cycle
    for(int i=-ghosts + isX; i<nx+ghosts - isX; i++){
        for(int j=-ghosts + isY; j<ny+ghosts - isY; j++){
            auto wL =(*this)[i-isX,j-isY], wR = (*this)[i+isX,j+isY];
            TVD::MUSCL(wL, _L[i,j], (*this)[i,j], _R[i,j],wR, dt_dL);
        }
    }
}
void Grid3D::computeHalfStates(Grid3D& _L, Grid3D& _R, double dt, int dim){
    //Dimension
    int isX = 0, isY = 0, isZ = 0;
    switch(dim % 3){
        case 0: isX = 1; break;
        case 1: isY = 1; break;
        case 2: isZ = 1; break;
    }
    double dt_dL = dt/(isX*dx + isY*dy + isZ*dz);
    //Cycle
    for(int i=-ghosts + isX; i<nx+ghosts - isX; i++){
        for(int j=-ghosts + isY; j<ny+ghosts - isY; j++){
            for(int k=-ghosts + isZ; k<nz+ghosts - isZ; k++){
                auto wL =(*this)[i-isX,j-isY,k-isZ], wR = (*this)[i+isX,j+isY,k+isZ];
                TVD::MUSCL(wL, _L[i,j,k], (*this)[i,j,k], _R[i,j,k],wR, dt_dL);
            }
        }
    }
}


//MARK: 2D Unsplit Step

void Grid2D::advanceXY(double dt){
    Grid2D _L(nx,ny,dx,dy,ghosts), _R(nx,ny,dx,dy,ghosts);//Half States
    FluxGrid2D F_X(nx+1, ny), F_Y(nx, ny+1); //Fluxes

    boundary.apply(*this);

    //Compute X Fluxes
    computeHalfStates(_L, _R, dt, 0);
    for(int i=0; i<=nx; i++){
        for(int j=0; j<ny; j++){
            auto L = _R[i-1,j], R = _L[i,j];
            F_X[i,j] = Riemann(L,R).flux();
        }
    }
    //Compute Y fluxes
    computeHalfStates(_L, _R, dt, 1);
    for(int i=0; i<nx; i++){
        for(int j=0; j<=ny; j++){
            auto L = _R[i,j-1], R = _L[i,j];
            L.swapXY(); R.swapXY();
            F_Y[i,j] = Riemann(L,R).flux();
            F_Y[i,j].swapXY();
        }
    }
    //Apply all Fluxes
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            (*this)[i,j] += (dt/dx) * (F_X[i,j] - F_X[i+1,j]);
            (*this)[i,j] += (dt/dy) * (F_Y[i,j] - F_Y[i,j+1]);
        }
    }
}
//MARK: 3D Unsplit Advancement

void Grid3D::advanceXYZ(double dt){
    Grid3D _L(nx,ny,nz,dx,dy,dz,ghosts), _R(nx,ny,nz,dx,dy,dz,ghosts);
    FluxGrid3D F_X(nx+1, ny,nz), F_Y(nx,ny+1,nz), F_Z(nx,ny,nz+1);
    boundary.apply(*this);

    //Compute X fluxes
    computeHalfStates(_L, _R, dt, 0);
    for(int i=0; i<=nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                auto L = _R[i-1,j,k], R = _L[i,j,k];
                F_X[i,j,k] = Riemann(L,R).flux();
            }
        }
        
    }
    //Compute Y fluxes
    computeHalfStates(_L, _R, dt, 1);
    for(int i=0; i<nx; i++){
        for(int j=0; j<=ny; j++){
            for(int k=0; k<nz; k++){
                auto L = _R[i,j-1,k], R = _L[i,j,k];
                L.swapXY(); R.swapXY();
                F_Y[i,j,k] = Riemann(L,R).flux();
                F_Y[i,j,k].swapXY();
            }
        }
    }
    //Compute Z fluxes
    computeHalfStates(_L, _R, dt, 2);
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<=nz; k++){
                //Z
                auto L = _R[i,j,k-1], R = _L[i,j,k];
                L.swapXZ(); R.swapXZ();
                F_Z[i,j,k] = Riemann(L,R).flux();
                F_Z[i,j,k].swapXZ();
            }
        }
    }
    //Apply Fluxes
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                (*this)[i,j,k] += (dt/dx) * (F_X[i,j,k] - F_X[i+1,j,k]);
                (*this)[i,j,k] += (dt/dy) * (F_Y[i,j,k] - F_Y[i,j+1,k]);
                (*this)[i,j,k] += (dt/dz) * (F_Z[i,j,k] - F_Z[i,j,k+1]);
            }
        }
    }
}
