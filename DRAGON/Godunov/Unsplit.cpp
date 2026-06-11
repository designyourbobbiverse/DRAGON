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
    

//MARK: CFL

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

//MARK: Flux Calculation
void computeFlux_X(const Grid2D& _L, const Grid2D& _R, FluxGrid2D& F, int xL, int xR, int yL, int yR){
    //Compute X Fluxes
    for(int i=xL; i<=xR; i++){
        for(int j=yL; j<yR; j++){
            auto L = _R[i-1,j], R = _L[i,j];
            F[i,j] = Riemann(L,R).flux();
        }
    }
}
void computeFlux_Y(const Grid2D& _L, const Grid2D& _R, FluxGrid2D& F, int xL, int xR, int yL, int yR){
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<=yR; j++){
            auto L = _R[i,j-1], R = _L[i,j];
            L.swapXY(); R.swapXY();
            F[i,j] = Riemann(L,R).flux();
            F[i,j].swapXY();
        }
    }
}
void computeFlux_X(const Grid3D& _L, const Grid3D& _R, FluxGrid3D& F, int xL, int xR, int yL, int yR, int zL, int zR){
    for(int i=xL; i<=xR; i++){
        for(int j=yL; j<yR; j++){
            for(int k=zL; k<zR; k++){
                auto L = _R[i-1,j,k], R = _L[i,j,k];
                F[i,j,k] = Riemann(L,R).flux();
            }
        }
    }
}
void computeFlux_Y(const Grid3D& _L, const Grid3D& _R, FluxGrid3D& F, int xL, int xR, int yL, int yR, int zL, int zR){
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<=yR; j++){
            for(int k=zL; k<zR; k++){
                auto L = _R[i,j-1,k], R = _L[i,j,k];
                L.swapXY(); R.swapXY();
                F[i,j,k] = Riemann(L,R).flux();
                F[i,j,k].swapXY();
            }
        }
    }
}
void computeFlux_Z(const Grid3D& _L, const Grid3D& _R, FluxGrid3D& F, int xL, int xR, int yL, int yR, int zL, int zR){
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<yR; j++){
            for(int k=zL; k<=zR; k++){
                auto L = _R[i,j,k-1], R = _L[i,j,k];
                L.swapXZ(); R.swapXZ();
                F[i,j,k] = Riemann(L,R).flux();
                F[i,j,k].swapXZ();
            }
        }
    }
}

//MARK: TVD
inline void dimensionFlags(int dim, int &isX, int &isY){
    isX = dim%2 == 0 ? 1 : 0;
    isY = dim%2 == 1 ? 1 : 0;
}
inline void dimensionFlags(int dim, int &isX, int &isY, int& isZ){
    isX = dim%3 == 0 ? 1 : 0;
    isY = dim%3 == 1 ? 1 : 0;
    isZ = dim%3 == 2 ? 1 : 0;
}

void Grid2D::computeHalfStates(Grid2D& _L, Grid2D& _R, double dt, int dim){
    //Dimension
    int isX, isY; dimensionFlags(dim, isX, isY);
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
    int isX, isY,isZ; dimensionFlags(dim, isX, isY,isZ);
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

//MARK: CTU Corrections
void correctState(Grid2D& _L, Grid2D& _R, const FluxGrid2D& F, double dt_dL, int xL, int xR, int yL, int yR, int dim){
    int isX, isY; dimensionFlags(dim, isX, isY);
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<yR; j++){
            auto trans = (F[i+isX,j+isY] - F[i,j]) * dt_dL;
            _L[i,j] = _L[i,j] -  trans;
            _R[i,j] = _R[i,j] -  trans;
        }
    }
}
void correctState(const Grid3D& _L0, const Grid3D& _R0, Grid3D& _L, Grid3D& _R, const FluxGrid3D& F, double dt_dL, int xL, int xR, int yL, int yR, int zL, int zR, int dim){
    int isX, isY, isZ; dimensionFlags(dim, isX, isY,isZ);
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<yR; j++){
            for(int k=zL; k<zR; k++){
                auto trans = (F[i+isX,j+isY,k+isZ] - F[i,j,k]) * dt_dL;
                _L[i,j,k] = _L0[i,j,k] -  trans;
                _R[i,j,k] = _R0[i,j,k] -  trans;
            }
        }
    }
}

//MARK: 2D Unsplit Step

void Grid2D::advanceXY(double dt){
    Grid2D _xL(nx,ny,dx,dy), _xR(nx,ny,dx,dy);//x Half States
    Grid2D _yL(nx,ny,dx,dy), _yR(nx,ny,dx,dy);//y Half States
    FluxGrid2D F_X(nx,ny,1), F_Y(nx,ny,1); //Fluxes
    
    boundary.apply(*this);

    //Compute Half States
    computeHalfStates(_xL, _xR, dt, 0);//_xLR needs (-1...nx, -1...ny)
    computeHalfStates(_yL, _yR, dt, 1);//_yLR needs (-1...nx, -1...ny)
#ifdef CTU
    //Compute preliminary Fluxes
    computeFlux_X(_xL, _xR, F_X, 0, nx, -1, ny+1);  //F_X needs (0...nx, -1...ny)
    computeFlux_Y(_yL, _yR, F_Y, -1, nx+1, 0, ny);  //F_Y needs (-1...nx, 0...ny)
    //Correct states
    correctState(_xL, _xR, F_Y, (0.5*dt/dy), -1, nx+1, 0, ny, 1); //_xLR needs (-1...nx, 0...ny-1)
    correctState(_yL, _yR, F_X, (0.5*dt/dx), 0, nx, -1, ny+1, 0); //_yLR needs (0...nx-1, -1...ny)
#endif
    //Compute Fluxes
    computeFlux_X(_xL, _xR, F_X, 0, nx, 0, ny); //F_X needs (0...nx, 0...ny-1)
    computeFlux_Y(_yL, _yR, F_Y, 0, nx, 0, ny); //F_Y needs (0...nx-1, 0...ny)
   
    //Apply all Fluxes
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            (*this)[i,j] += (dt/dx) * (F_X[i,j] - F_X[i+1,j]);
            (*this)[i,j] += (dt/dy) * (F_Y[i,j] - F_Y[i,j+1]);
        }
    }
}

//MARK: 3D Unsplit Step

void Grid3D::advanceXYZ(double dt){
    Grid3D _xL(nx,ny,nz,dx,dy,dz), _xR(nx,ny,nz,dx,dy,dz);
    Grid3D _yL(nx,ny,nz,dx,dy,dz), _yR(nx,ny,nz,dx,dy,dz);
    Grid3D _zL(nx,ny,nz,dx,dy,dz), _zR(nx,ny,nz,dx,dy,dz);
    FluxGrid3D F_X(nx, ny,nz,1), F_Y(nx,ny,nz,1), F_Z(nx,ny,nz,1);
    
    boundary.apply(*this);

    //Compute Half States
    computeHalfStates(_xL, _xR, dt, 0);
    computeHalfStates(_yL, _yR, dt, 1);
    computeHalfStates(_zL, _zR, dt, 2);
#ifdef CTU
    //Compute Face Fluxes
    computeFlux_X(_xL, _xR, F_X, 0, nx, -1, ny+1, -1, nz+1); //F_X needs (0...nx, -1...ny, -1...nz)
    computeFlux_Y(_yL, _yR, F_Y, -1, nx+1, 0, ny, -1, nz+1); //F_Y needs (-1...nx, 0...ny, -1...nz)
    computeFlux_Z(_zL, _zR, F_Z, -1, nx+1, -1, ny+1, 0, nz); //F_Z needs (-1...nx, -1...ny, 0...nz)
    //Compute Edge Corrections
    Grid3D _xyL(nx,ny,nz,dx,dy,dz), _xyR(nx,ny,nz,dx,dy,dz); //_xyLR needs (-1...nx, 0...ny-1, -1...nz)
    correctState(_xL, _xR, _xyL, _xyR, F_Y, (0.5*dt/dy) , -1, nx+1, 0, ny, -1, nz+1, 1);
    Grid3D _xzL(nx,ny,nz,dx,dy,dz), _xzR(nx,ny,nz,dx,dy,dz);//_xzLR needs (-1...nx, -1...ny, 0...nz-1)
    correctState(_xL, _xR, _xzL, _xzR, F_Z, (0.5*dt/dz) , -1, nx+1, -1, ny+1, 0, nz, 2);
    Grid3D _yxL(nx,ny,nz,dx,dy,dz), _yxR(nx,ny,nz,dx,dy,dz);//_yxLR needs (0...nx-1, -1...ny, -1...nz)
    correctState(_yL, _yR, _yxL, _yxR, F_X, (0.5*dt/dx), 0, nx, -1, ny+1, -1, nz+1, 0);
    Grid3D _yzL(nx,ny,nz,dx,dy,dz), _yzR(nx,ny,nz,dx,dy,dz);//_yzLR needs (-1...nx, -1...ny, 0...nz-1)
    correctState(_yL, _yR, _yzL, _yzR, F_Z, (0.5*dt/dz), -1, nx+1, -1, ny+1, 0, nz, 2);
    Grid3D _zxL(nx,ny,nz,dx,dy,dz), _zxR(nx,ny,nz,dx,dy,dz); //_zxLR needs (0...nx-1, -1...ny, -1...nz)
    correctState(_zL, _zR, _zxL, _zxR, F_X, (0.5*dt/dx), 0, nx, -1, ny+1, -1, nz+1, 0);
    Grid3D _zyL(nx,ny,nz,dx,dy,dz), _zyR(nx,ny,nz,dx,dy,dz);//_zyLR needs (-1...nx, 0...ny-1, -1...nz)
    correctState(_zL, _zR, _zyL, _zyR, F_Y, (0.5*dt/dx), -1, nx+1, 0, ny, -1, nz+1, 1);
    //Compute Edge Fluxes
    FluxGrid3D F_Xy(nx, ny,nz,1), F_Xz(nx, ny,nz,1);
    computeFlux_X(_xyL, _xyR, F_Xy, 0, nx, 0, ny, -1,nz+1); //F_Xy needs (0...nx, 0...ny-1, -1...nz)
    computeFlux_X(_xzL, _xzR, F_Xz, 0, nx, -1, ny+1, 0,nz); //F_Xz needs (0...nx, -1...ny, 0...nz-1)
    FluxGrid3D F_Yx(nx, ny,nz,1), F_Yz(nx, ny,nz,1);
    computeFlux_Y(_yxL, _yxR, F_Yx, 0, nx, 0, ny, -1, nz+1); //F_Yx needs (0...nx-1, 0...ny, -1...nz)
    computeFlux_Y(_yzL, _yzR, F_Yz, -1, nx+1, 0, ny, 0, nz); //F_Yz needs (-1...nx, 0...ny, 0...nz-1)
    FluxGrid3D F_Zx(nx,ny,nz,1), F_Zy(nx, ny,nz,1);
    computeFlux_Z(_zxL, _zxR, F_Zx, 0, nx, -1, ny+1, 0, nz); //F_Zx needs (0...nx-1, -1...ny, 0...nz)
    computeFlux_Z(_zyL, _zyR, F_Zy, -1, nx+1, 0, ny, 0, nz); //F_Zy needs (-1...nx, 0...ny-1, 0...nz)
    //Apply Corner Corrections
    correctState(_xL, _xR, _xL,_xR, F_Yz, (0.5*dt/dy), -1, nx+1, 0, ny, 0, nz, 1); //_xLR needs (-1...nx, 0...ny-1, 0...nz-1)
    correctState(_xL, _xR, _xL,_xR, F_Zy, (0.5*dt/dz), -1, nx+1, 0, ny, 0, nz, 2);
    correctState(_yL, _yR, _yL,_yR, F_Xz, (0.5*dt/dx), 0, nx, -1, ny+1, 0, nz, 0);//_yLR needs (0...nx-1, -1...ny, 0...nz-1)
    correctState(_yL, _yR, _yL,_yR, F_Zx, (0.5*dt/dz), 0, nx, -1, ny+1, 0, nz, 2);
    correctState(_zL, _zR, _zL,_zR, F_Xy, (0.5*dt/dx), 0, nx, 0, ny, -1, nz+1, 0);//_zLR needs (0...nx-1, 0...ny-1, -1...nz)
    correctState(_zL, _zR, _zL,_zR, F_Yx, (0.5*dt/dy), 0, nx, 0, ny, -1, nz+1, 1);
#endif
    //Compute Fluxes
    computeFlux_X(_xL, _xR, F_X, 0, nx, 0, ny, 0, nz); //F_X needs (0...nx, 0...ny-1, 0...nz-1)
    computeFlux_Y(_yL, _yR, F_Y, 0, nx, 0, ny, 0, nz); //F_Y needs (0...nx-1, 0...ny, 0...nz-1)
    computeFlux_Z(_zL, _zR, F_Z, 0, nx, 0, ny, 0, nz); //F_Z needs (0...nx-1, 0...ny-1, 0...nz)
    
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
