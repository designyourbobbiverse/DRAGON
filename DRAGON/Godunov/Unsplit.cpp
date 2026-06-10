//
//  Unsplit.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 10/06/2026.
//

#include "Grid.hpp"
#include "Boundary.hpp"
#include "Riemann.hpp"
#include "CFL.hpp"
#include <math.h>
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

//MARK: 2D Unsplit Step
void Grid2D::advanceXY(double dt){
    FluxGrid2D F_X(nx+1, ny),  F_Y(nx, ny+1);
    boundary.apply(*this);
    //Compute X fluxes
    for(int i=0; i<=nx; i++){
        for(int j=0; j<ny; j++){
            if(ghosts == 0 && (i == 0 || i == nx)) { //If no ghosts, fake an outflow boundary
                F_X[i,j] = (*this)[i == 0 ? i : i-1,j].flux();
            } else {
                auto L = (*this)[i-1,j], R = (*this)[i,j];
                F_X[i,j] = Riemann(L,R).flux();
            }
        }
    }
    //Compute Y fluxes
    for(int i=0; i<nx; i++){
        for(int j=0; j<=ny; j++){
            if(ghosts == 0 && (j == 0 || j == ny)) { //If no ghosts, fake an outflow boundary
                F_Y[i,j] = (*this)[i, j == 0 ? j : j-1].flux();
            } else {
                auto L = (*this)[i,j-1], R = (*this)[i,j];
                L.swapXY(); R.swapXY();
                F_Y[i,j] = Riemann(L,R).flux();
                F_Y[i,j].swapXY();
            }
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
    FluxGrid3D F_X(nx+1, ny,nz), F_Y(nx,ny+1,nz), F_Z(nx,ny,nz+1);
    boundary.apply(*this);

    //Compute X fluxes
    for(int i=0; i<=nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                if(ghosts == 0 && (i == 0 || i == nx)) { //If no ghosts, fake an outflow boundary
                    F_X[i,j,k] = (*this)[i == 0 ? i : i-1,j,k].flux();
                } else {
                    auto L = (*this)[i-1,j,k], R = (*this)[i,j,k];
                    F_X[i,j,k] = Riemann(L,R).flux();
                }
            }
        }
        
    }
    //Compute Y fluxes
    for(int i=0; i<nx; i++){
        for(int j=0; j<=ny; j++){
            for(int k=0; k<nz; k++){
                if(ghosts == 0 && (j == 0 || j == ny)) { //If no ghosts, fake an outflow boundary
                    F_Y[i,j,k] = (*this)[i, j == 0 ? j : j-1,k].flux();
                } else {
                    auto L = (*this)[i,j-1,k], R = (*this)[i,j,k];
                    L.swapXY(); R.swapXY();
                    F_Y[i,j,k] = Riemann(L,R).flux();
                    F_Y[i,j,k].swapXY();
                }
            }
        }
    }
    //Compute Z fluxes
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<=nz; k++){
                //Z
                if(ghosts == 0 && (k == 0 || k == nz)) { //If no ghosts, fake an outflow boundary
                    F_Z[i,j,k] = (*this)[i, j, k == 0 ? k : k-1].flux();
                } else {
                    auto L = (*this)[i,j,k-1], R = (*this)[i,j,k];
                    L.swapXZ(); R.swapXZ();
                    F_Z[i,j,k] = Riemann(L,R).flux();
                    F_Z[i,j,k].swapXZ();
                }
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
