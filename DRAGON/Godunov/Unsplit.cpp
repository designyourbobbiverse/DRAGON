//
//  Unsplit.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 10/06/2026.
//  Implementation based in part on  Toro (2009). https://doi.org/10.1007/b79761
//      Colella (1990). https://doi.org/10.1016/0021-9991(90)90233-Q
//  CT Implementaiton based in part on https://doi.org/10.1088/0067-0049/182/2/519
//      Evans and Hawley (1988). https://doi.org/10.1086/166684
//      Gardiner and Stone (2005). https://doi.org/10.1016/j.jcp.2004.11.016


#include "Grid.hpp"
#include "Riemann.hpp"
#include "Config.h"
#include "CFL.hpp"
#include "TVD.hpp"
#include <cassert>
#include <iostream>


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
void Grid2D::advance_unsplit(double dt, bool check_cfl){
    while(dt > Timestep_Tolerance){
        //CFL Time Constraint
        double t1 = check_cfl ? std::min(dt,CFL::cfl_time(*this)) : dt;
        dt -= t1;
        //Advance
        advanceXY(t1);
    }
}
void Grid3D::advance_unsplit(double dt, bool check_cfl){
    while(dt > Timestep_Tolerance){
        //CFL Time Constraint
        double t1 = check_cfl ? std::min(dt,CFL::cfl_time(*this)) : dt;
        dt -= t1;
        //Advance
        advanceXYZ(t1);
    }
}

//MARK: 2D Unsplit Step
typedef ExtendedArray2D<PrimitiveState> FluidArray2D;
typedef ExtendedArray2D<ConservativeState> FluxArray2D;

void computeFlux_X(const FluidArray2D& _L, const FluidArray2D& _R, FluxArray2D& F, int xL, int xR, int yL, int yR, double dt_dx);
void computeFlux_Y(const FluidArray2D& _L, const FluidArray2D& _R, FluxArray2D& F, int xL, int xR, int yL, int yR, double dt_dy);
void computeHalfStates_X(FluidArray2D& _L, const Grid2D& _W,  FluidArray2D& _R, double dt);
void computeHalfStates_Y( FluidArray2D& _L,const Grid2D& _W, FluidArray2D& _R, double dt);
void correctState(FluidArray2D& _L, FluidArray2D& _R, const FluxArray2D& F, double dt_dL, int xL, int xR, int yL, int yR, int dim);
#ifdef MHD
typedef ExtendedArray2D<vec3> MagneticArray2D;
void computeFaceFields(const MagneticArray2D& _A, MagneticArray2D& _B, double dx, double dy);
void copyFaceFields_X( FluidArray2D& _L,const MagneticArray2D& _B, FluidArray2D& _R);
void copyFaceFields_Y( FluidArray2D& _L,const MagneticArray2D& _B, FluidArray2D& _R);
void updatePotential(MagneticArray2D& _B, const FluxArray2D& F_X,const FluxArray2D& F_Y, double dt);
#endif


void Grid2D::advanceXY(double dt){
    int nx = w.getSizeX(), ny = w.getSizeY(), ghosts = w.getGhosts();
    FluidArray2D _xL(nx,ny,ghosts), _xR(nx,ny,ghosts);//x Half States
    FluidArray2D _yL(nx,ny,ghosts), _yR(nx,ny,ghosts);//y Half States
    FluxArray2D F_X(nx,ny,ghosts), F_Y(nx,ny,ghosts); //Fluxes
    
    boundary.apply(*this);

    #ifdef MHD
    //Compute the face fields
    MagneticArray2D B(nx+1,ny+1,ghosts);
    computeFaceFields(A, B, dx, dy);
    //Set normal fields at faces to match face-centred fields
    copyFaceFields_X(_xL, B, _xR);
    copyFaceFields_Y(_yL, B, _yR);
    #endif
    
    //Compute Half States
    computeHalfStates_X(_xL, (*this), _xR, dt);//_xLR needs (-1...nx, -1...ny), (-2...nx+1, -1...ny)for MHD
    computeHalfStates_Y(_yL, (*this), _yR, dt);//_yLR needs (-1...nx, -1...ny), (-1...nx, -2...ny+1) for MHD
    
#ifdef CTU //Colella (1990). https://doi.org/10.1016/0021-9991(90)90233-Q
    //Compute preliminary Fluxes
    computeFlux_X(_xL, _xR, F_X, -1, nx+1, -1, ny+1, dt/dx);  //F_X needs (0...nx, -1...ny), (-1...nx+1, -1...ny)for MHD
    computeFlux_Y(_yL, _yR, F_Y, -1, nx+1, -1, ny+1, dt/dy);  //F_Y needs (-1...nx, 0...ny), (-1...nx, -1...ny+1) for MHD
    //Correct states
    correctState(_xL, _xR, F_Y, (0.5*dt/dy), -1, nx+1, -1, ny+1, 1); //_xLR needs (-1...nx, 0...ny-1), (-1...nx, -1...ny) for MHD
    correctState(_yL, _yR, F_X, (0.5*dt/dx), -1, nx+1, -1, ny+1, 0); //_yLR needs (0...nx-1, -1...ny), (-1...nx, -1...ny) for MHD
    
    #ifdef MHD//Restore face-normal fields
    copyFaceFields_X(_xL, B, _xR);
    copyFaceFields_Y(_yL, B, _yR);
    #endif
#endif
    //Compute Fluxes
    computeFlux_X(_xL, _xR, F_X, 0, nx, -1, ny+1, dt/dx); //F_X needs (0...nx, 0...ny-1), (0...nx, -1...ny) for MHD
    computeFlux_Y(_yL, _yR, F_Y, -1, nx+1, 0, ny, dt/dy); //F_Y needs (0...nx-1, 0...ny), (-1...nx, 0...ny) for MHD
    
    //Apply all Fluxes
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            w[i,j] += (dt/dx) * (F_X[i,j] - F_X[i+1,j]);
            w[i,j] += (dt/dy) * (F_Y[i,j] - F_Y[i,j+1]);
        }
    }
    
    #ifdef MHD //Do CT Update
    updatePotential(A, F_X, F_Y, dt);
    computeBodyAveragedFields();
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
#ifdef MHD
typedef ExtendedArray3D<vec3> MagneticArray3D;
void computeFaceFields(const MagneticArray3D& _A, MagneticArray3D& _B, double dx, double dy, double dz);
void copyFaceFields_X( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R);
void copyFaceFields_Y( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R);
void copyFaceFields_Z( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R);
void updatePotential(MagneticArray3D& _B, const FluxArray3D& F_X,const FluxArray3D& F_Y,const FluxArray3D& F_Z, double dt);
#endif

void Grid3D::advanceXYZ(double dt){
    int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ(), ghosts = w.getGhosts();
    FluidArray3D _xL(nx,ny,nz,ghosts), _xR(nx,ny,nz,ghosts); //X half States
    FluidArray3D _yL(nx,ny,nz,ghosts), _yR(nx,ny,nz,ghosts); //Y half states
    FluidArray3D _zL(nx,ny,nz,ghosts), _zR(nx,ny,nz,ghosts); //Z half states
    FluxArray3D F_X(nx, ny,nz,ghosts), F_Y(nx,ny,nz,ghosts), F_Z(nx,ny,nz,ghosts); //Fluxes
    
    boundary.apply(*this);
    
    
    #ifdef MHD
    computeBodyAveragedFields();
    //Compute the face fields
    MagneticArray3D B(nx+1,ny+1,nz+1,ghosts);
    computeFaceFields(A, B, dx, dy, dz);
    //Set normal fields at faces to match face-centred fields
    copyFaceFields_X(_xL, B, _xR);
    copyFaceFields_Y(_yL, B, _yR);
    copyFaceFields_Z(_zL, B, _zR);
    #endif
        
    //Compute Half States
    computeHalfStates_X(_xL, (*this), _xR, dt);
    computeHalfStates_Y(_yL, (*this), _yR, dt);
    computeHalfStates_Z(_zL, (*this), _zR, dt);

#ifdef CTU //Colella (1990). https://doi.org/10.1016/0021-9991(90)90233-Q
    //Compute preliminary face fluxes
    computeFlux_X(_xL, _xR, F_X, -1, nx+1, -2, ny+2, -2, nz+2, dt/dx); //F_X needs (0...nx, -1...ny, -1...nz)
    computeFlux_Y(_yL, _yR, F_Y, -2, nx+2, -1, ny+1, -2, nz+2, dt/dy); //F_Y needs (-1...nx, 0...ny, -1...nz)
    computeFlux_Z(_zL, _zR, F_Z, -2, nx+2, -2, ny+2, -1, nz+1, dt/dz); //F_Z needs (-1...nx, -1...ny, 0...nz)
    //Apply first round of transverse corrections
    FluidArray3D _xyL(nx,ny,nz,ghosts), _xyR(nx,ny,nz,ghosts); //_xyLR needs (-1...nx, 0...ny-1, -1...nz)
    correctState(_xL, _xR, _xyL, _xyR, F_Y, (0.5*dt/dy) , -2, nx+2, -1, ny+1, -1, nz+1, 1);
    FluidArray3D _xzL(nx,ny,nz,ghosts), _xzR(nx,ny,nz,ghosts);//_xzLR needs (-1...nx, -1...ny, 0...nz-1)
    correctState(_xL, _xR, _xzL, _xzR, F_Z, (0.5*dt/dz) , -2, nx+2, -1, ny+1, -1, nz+1, 2);
    FluidArray3D _yxL(nx,ny,nz,ghosts), _yxR(nx,ny,nz,ghosts);//_yxLR needs (0...nx-1, -1...ny, -1...nz)
    correctState(_yL, _yR, _yxL, _yxR, F_X, (0.5*dt/dx), -1, nx+1, -2, ny+2, -1, nz+1, 0);
    FluidArray3D _yzL(nx,ny,nz,ghosts), _yzR(nx,ny,nz,ghosts);//_yzLR needs (-1...nx, -1...ny, 0...nz-1)
    correctState(_yL, _yR, _yzL, _yzR, F_Z, (0.5*dt/dz), -1, nx+1, -2, ny+2, -1, nz+1, 2);
    FluidArray3D _zxL(nx,ny,nz,ghosts), _zxR(nx,ny,nz,ghosts); //_zxLR needs (0...nx-1, -1...ny, -1...nz)
    correctState(_zL, _zR, _zxL, _zxR, F_X, (0.5*dt/dx), -1, nx+1, -1, ny+1, -2, nz+2, 0);
    FluidArray3D _zyL(nx,ny,nz,ghosts), _zyR(nx,ny,nz,ghosts);//_zyLR needs (-1...nx, 0...ny-1, -1...nz)
    correctState(_zL, _zR, _zyL, _zyR, F_Y, (0.5*dt/dy), -1, nx+1, -1, ny+1, -2, nz+2, 1);
    //Recompute fluxes
    FluxArray3D F_Xy(nx, ny,nz,2), F_Xz(nx, ny,nz,2);
    computeFlux_X(_xyL, _xyR, F_Xy, -1, nx+1, -1, ny+1, -1, nz+1, dt/dx); //F_Xy needs (0...nx, 0...ny-1, -1...nz)
    computeFlux_X(_xzL, _xzR, F_Xz, -1, nx+1, -1, ny+1, -1, nz+1, dt/dx); //F_Xz needs (0...nx, -1...ny, 0...nz-1)
    FluxArray3D F_Yx(nx, ny,nz,2), F_Yz(nx, ny,nz,2);
    computeFlux_Y(_yxL, _yxR, F_Yx, -1, nx+1, -1, ny+1, -1, nz+1, dt/dy); //F_Yx needs (0...nx-1, 0...ny, -1...nz)
    computeFlux_Y(_yzL, _yzR, F_Yz, -1, nx+1, -1, ny+1, -1, nz+1, dt/dy); //F_Yz needs (-1...nx, 0...ny, 0...nz-1)
    FluxArray3D F_Zx(nx,ny,nz,2), F_Zy(nx, ny,nz,2);
    computeFlux_Z(_zxL, _zxR, F_Zx, -1, nx+1, -1, ny+1, -1, nz+1, dt/dz); //F_Zx needs (0...nx-1, -1...ny, 0...nz)
    computeFlux_Z(_zyL, _zyR, F_Zy, -1, nx+1, -1, ny+1, -1, nz+1, dt/dz); //F_Zy needs (-1...nx, 0...ny-1, 0...nz)
    //Apply second round of transverse corrections
    correctState(_xL, _xR, _xL,_xR, F_Yz, (0.5*dt/dy), -1, nx+1, -1, ny+1, -1, nz+1, 1); //_xLR needs (-1...nx, -1...ny, -1...nz)
    correctState(_xL, _xR, _xL,_xR, F_Zy, (0.5*dt/dz), -1, nx+1, -1, ny+1, -1, nz+1, 2);
    correctState(_yL, _yR, _yL,_yR, F_Xz, (0.5*dt/dx), -1, nx+1, -1, ny+1, -1, nz+1, 0);//_yLR needs (-1...nx, -1...ny, -1...nz)
    correctState(_yL, _yR, _yL,_yR, F_Zx, (0.5*dt/dz), -1, nx+1, -1, ny+1, -1, nz+1, 2);
    correctState(_zL, _zR, _zL,_zR, F_Xy, (0.5*dt/dx), -1, nx+1, -1, ny+1, -1, nz+1, 0);//_zLR needs (-1...nx, -1...ny, -1...nz)
    correctState(_zL, _zR, _zL,_zR, F_Yx, (0.5*dt/dy), -1, nx+1, -1, ny+1, -1, nz+1, 1);
    
    #ifdef MHD//Restore face-normal fields
    copyFaceFields_X(_xL, B, _xR);
    copyFaceFields_Y(_yL, B, _yR);
    copyFaceFields_Z(_zL, B, _zR);
    #endif
#endif
    
    //Compute Fluxes
    computeFlux_X(_xL, _xR, F_X, 0, nx, -1, ny+1, -1, nz+1, dt/dx); //F_X needs (0...nx, -1...ny, -1...nz)
    computeFlux_Y(_yL, _yR, F_Y, -1, nx+1, 0, ny, -1, nz+1, dt/dy); //F_Y needs (-1...nx, 0...ny, -1...nz)
    computeFlux_Z(_zL, _zR, F_Z, -1, nx+1, -1, ny+1, 0, nz, dt/dz); //F_Z needs (-1...nx, -1...ny, 0...nz)
    
    //Apply Fluxes
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){
                w[i,j,k] += (dt/dx) * (F_X[i,j,k] - F_X[i+1,j,k]);
                w[i,j,k] += (dt/dy) * (F_Y[i,j,k] - F_Y[i,j+1,k]);
                w[i,j,k] += (dt/dz) * (F_Z[i,j,k] - F_Z[i,j,k+1]);
            }
        }
    }
    
    #ifdef MHD
    updatePotential(A, F_X, F_Y, F_Z, dt);
    computeBodyAveragedFields();
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

//MARK: TVD
//Apply MUSCL over the entire grid
void computeHalfStates_X(FluidArray2D& _L, const Grid2D& _W, FluidArray2D& _R, double dt){
    const double dt_dL = dt/_W.dx;//Compute once
    const int nx = _W.getSizeX(), ny = _W.getSizeY(), g = _W.getGhosts();

    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
                #ifdef MHD
                double BL = _L[i,j].B.x, BR = _R[i,j].B.x;
                #endif
                _L[i,j] = _W[i,j]; _R[i,j] = _W[i,j];
                #ifdef MHD
                _L[i,j].B.x = BL; _R[i,j].B.x = BR;
                #endif
        }
    }
    for(int i=-g + 1; i<nx+g - 1; i++){
        for(int j=-g; j<ny+g; j++){
            auto wL =_W[i-1,j], wR = _W[i+1,j];
            TVD::MUSCL(wL, _L[i,j], _W[i,j], _R[i,j], wR, dt_dL);
        }
    }
}
void computeHalfStates_Y(FluidArray2D& _L, const Grid2D& _W, FluidArray2D& _R, double dt){
    const double dt_dL = dt/_W.dy;//Compute once
    const int nx = _W.getSizeX(), ny = _W.getSizeY(), g = _W.getGhosts();

    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
                #ifdef MHD
                double BL = _L[i,j].B.y, BR = _R[i,j].B.y;
                #endif
                _L[i,j] = _W[i,j]; _R[i,j] = _W[i,j];
                #ifdef MHD
                _L[i,j].B.y = BL; _R[i,j].B.y = BR;
                #endif
        }
    }
    for(int i=-g; i<nx+g; i++){
        for(int j=-g + 1; j<ny+g - 1; j++){
            #ifdef MHD//By has already set, but MUSCL needs to see it as Bx
            _L[i,j].swapXY(); _R[i,j].swapXY();
            #endif
            //Swap XY inputs to MUSCL, then swap output back
            auto wL =_W[i,j-1].swappedXY(), wR = _W[i,j+1].swappedXY();
            TVD::MUSCL(wL, _L[i,j], _W[i,j].swappedXY(), _R[i,j],wR, dt_dL);
            _L[i,j].swapXY(); _R[i,j].swapXY();
        }
    }
}
void computeHalfStates_X(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, double dt){
    const double dt_dL = dt/_W.dx;//Compute once
    const int nx = _W.getSizeX(), ny = _W.getSizeY(),nz = _W.getSizeZ(), g = _W.getGhosts();

    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g; k<nz+g; k++){
                #ifdef MHD
                double BL = _L[i,j,k].B.x, BR = _R[i,j,k].B.x;
                #endif
                _L[i,j,k] = _W[i,j,k]; _R[i,j,k] = _W[i,j,k];
                #ifdef MHD
                _L[i,j,k].B.x = BL; _R[i,j,k].B.x = BR;
                #endif
            }
        }
    }
    for(int i=-g + 1; i<nx+g - 1; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g; k<nz+g; k++){
                auto wL = _W[i-1,j,k], wR = _W[i+1,j,k];
                TVD::MUSCL(wL, _L[i,j,k], _W[i,j,k], _R[i,j,k], wR, dt_dL);
            }
        }
    }
}
void computeHalfStates_Y(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, double dt){
    const double dt_dL = dt/_W.dy;//Compute once
    const int nx = _W.getSizeX(), ny = _W.getSizeY(),nz = _W.getSizeZ(), g = _W.getGhosts();
    
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g; k<nz+g; k++){
                #ifdef MHD
                double BL = _L[i,j,k].B.y, BR = _R[i,j,k].B.y;
                #endif
                _L[i,j,k] = _W[i,j,k]; _R[i,j,k] = _W[i,j,k];
                #ifdef MHD
                _L[i,j,k].B.y = BL; _R[i,j,k].B.y = BR;
                #endif
            }
        }
    }
    for(int i=-g; i<nx+g; i++){
        for(int j=-g + 1; j<ny+g - 1; j++){
            for(int k=-g; k<nz+g; k++){
                #ifdef MHD//Bz has already set, but MUSCL needs to see it as Bx
                _L[i,j,k].swapXY(); _R[i,j,k].swapXY();
                #endif
                //Swap XY inputs to MUSCL, then swap output back
                auto wL =_W[i,j-1,k].swappedXY(), wR = _W[i,j+1,k].swappedXY();
                TVD::MUSCL(wL, _L[i,j,k], _W[i,j,k].swappedXY(), _R[i,j,k],wR, dt_dL);
                _L[i,j,k].swapXY(); _R[i,j,k].swapXY();
            }
        }
    }
}
void computeHalfStates_Z(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, double dt){
    const double dt_dL = dt/_W.dz;//Compute once
    const int nx = _W.getSizeX(), ny = _W.getSizeY(),nz = _W.getSizeZ(), g = _W.getGhosts();
    
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g; k<nz+g; k++){
                #ifdef MHD
                double BL = _L[i,j,k].B.z, BR = _R[i,j,k].B.z;
                #endif
                _L[i,j,k] = _W[i,j,k]; _R[i,j,k] = _W[i,j,k];
                #ifdef MHD
                _L[i,j,k].B.z = BL; _R[i,j,k].B.z = BR;
                #endif
            }
        }
    }
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            for(int k=-g + 1; k<nz+g - 1; k++){
                #ifdef MHD//Bz has already set, but MUSCL needs to see it as Bx
                _L[i,j,k].swapXZ(); _R[i,j,k].swapXZ();
                #endif
                //Swap XZ inputs to MUSCL, then swap output back
                auto wL = _W[i,j,k-1].swappedXZ(), wR = _W[i,j,k+1].swappedXZ();
                TVD::MUSCL(wL, _L[i,j,k], _W[i,j,k].swappedXZ(), _R[i,j,k],wR, dt_dL);
                _L[i,j,k].swapXZ(); _R[i,j,k].swapXZ();
            }
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
#ifdef MHD
//MARK: Face-Centred Field Copy
void copyFaceFields_X( FluidArray2D& _L,const MagneticArray2D& _B, FluidArray2D& _R){
    const int nx = _L.getSizeX(), ny = _L.getSizeY(), g = _B.getGhosts();
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            _L[i,j].B.x = _B[i,j].x;
            _R[i,j].B.x = _B[i+1,j].x;
        }
    }
}
void copyFaceFields_Y( FluidArray2D& _L,const MagneticArray2D& _B, FluidArray2D& _R){
    const int nx = _L.getSizeX(), ny = _L.getSizeY(), g = _B.getGhosts();
    for(int i=-g; i<nx+g; i++){
        for(int j=-g; j<ny+g; j++){
            _L[i,j].B.y = _B[i,j].y;
            _R[i,j].B.y = _B[i,j+1].y;
        }
    }
}
void copyFaceFields_X( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R){
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
void copyFaceFields_Y( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R){
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
void copyFaceFields_Z( FluidArray3D& _L,const MagneticArray3D& _B, FluidArray3D& _R){
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
//MARK: Constrained Transport

void updatePotential(MagneticArray2D& _A, const FluxArray2D& F_X,const FluxArray2D& F_Y, double dt){
    const int nx = _A.getSizeX()-1, ny = _A.getSizeY()-1;
    
    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++){
            _A[i,j].z -= (F_Y[i-1,j].B.x -  F_X[i,j-1].B.y + F_Y[i,j].B.x  - F_X[i,j].B.y) * 0.25 * dt;
        }
    }
}
void updatePotential(MagneticArray3D& _A, const FluxArray3D& F_X,const FluxArray3D& F_Y, const FluxArray3D& F_Z,  double dt){
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
void computeFaceFields(const MagneticArray2D& _A, MagneticArray2D& _B, double dx, double dy){
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
void computeFaceFields(const MagneticArray3D& _A, MagneticArray3D& _B, double dx, double dy, double dz){
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


void Grid2D::computeBodyAveragedFields(){
    const int nx = w.getSizeX(), ny = w.getSizeY();
    MagneticArray2D B(nx+1,ny+1,w.getGhosts());
    computeFaceFields(A, B, dx, dy);
    
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            ConservativeState U(w[i,j]);//Update the conservative element to keep quantities conserved
            U.B.x = (B[i,j].x + B[i+1,j].x)/2;
            U.B.y = (B[i,j].y + B[i,j+1].y)/2;
            w[i,j] = U;
        }
    }
}
void Grid3D::computeBodyAveragedFields(){
    const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ();
    MagneticArray3D B(nx+1,ny+1,nz+1,w.getGhosts());
    computeFaceFields(A, B, dx, dy, dz);

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


#endif
