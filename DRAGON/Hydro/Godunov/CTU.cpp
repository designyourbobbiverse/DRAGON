//
//  CTU.cpp
//  DRAGON/Hydro/Godunov
//
//  Created by Bobbie Markwick on 10/06/2026.
//  Implementation based in part on  Toro (2009). https://doi.org/10.1007/b79761
//      Colella (1990). https://doi.org/10.1016/0021-9991(90)90233-Q
//      Gardiner and Stone (2005). https://doi.org/10.1016/j.jcp.2004.11.016


#include "Grid.hpp"
#include "Unsplit.hpp"
#include "Riemann.hpp"
#include "Config.h"
#include "Constants.h"
#include "TVD.hpp"
#include "CT.hpp"
#include "DragonWing.hpp"

#include <iostream>

#ifdef CTU
#ifndef MUSCL_Hancock
#error CTU requires MUSCL reconstruction. Please enable MUSCL_Hancock in Config.h
#endif
#endif

//MARK: CTU 2D Hydro
void ctu_sweep_hydro(FluidArray2D& _xL, FluidArray2D& _xR, FluidArray2D& _yL, FluidArray2D& _yR, double dt_dx, double dt_dy){
    const int nx = _xL.getSizeX(), ny = _yL.getSizeY(), ghosts = _xL.getGhosts();
    //Compute preliminary Fluxes
        auto __fluxes = DRAGONWING::requestFluxArrays(2, nx, ny, ghosts);
    auto& F_X = *__fluxes[0];
    auto& F_Y = *__fluxes[1];
    computeFlux_X(_xL, _xR, F_X, -1, nx+1, -2, ny+2, dt_dx);
    computeFlux_Y(_yL, _yR, F_Y, -2, nx+2, -1, ny+1, dt_dy);
    
    //Correct states
    correctState(_xL, _xR, F_Y, (0.5*dt_dy), 1);
    correctState(_yL, _yR, F_X, (0.5*dt_dx), 0);
}


//MARK: CTU Hydro 12-Solve (3D)
void ctu_sweep_hydro(FluidArray3D& _xL, FluidArray3D& _xR, FluidArray3D& _yL, FluidArray3D& _yR, FluidArray3D& _zL, FluidArray3D& _zR, double dt_dx, double dt_dy, double dt_dz){
    const int nx = _xL.getSizeX(), ny = _xL.getSizeY(), nz = _xL.getSizeZ(), ghosts = _xL.getGhosts();
    
        auto __fluxes = DRAGONWING::requestFluxArrays(5, nx, ny, nz, ghosts);
    FluxArray3D& F_X = *__fluxes[0];
    FluxArray3D& F_Y = *__fluxes[1];
    FluxArray3D& F_Z = *__fluxes[2];
    //Compute preliminary face fluxes
    computeFlux_X(_xL, _xR, F_X, -2, nx+2, -2, ny+2, -2, nz+2, dt_dx);
    computeFlux_Y(_yL, _yR, F_Y, -2, nx+2, -2, ny+2, -2, nz+2, dt_dy);
    computeFlux_Z(_zL, _zR, F_Z, -2, nx+2, -2, ny+2, -2, nz+2, dt_dz);
        
    //Compute Edge-correct the fluxes
    FluxArray3D& F_Xz = *__fluxes[3];
    computeCTUFlux_X(_xL, _xR, F_Z, F_Xz, dt_dx, (0.5*dt_dz), 2);
    FluxArray3D& F_Yz = *__fluxes[4];
    computeCTUFlux_Y(_yL, _yR, F_Z, F_Yz, dt_dy, (0.5*dt_dz), 2);

    FluxArray3D& F_Xy = F_Z; //F_Z isn't used again before it gets recomputed, so we can resue it here
    computeCTUFlux_X(_xL, _xR, F_Y, F_Xy, dt_dx, (0.5*dt_dy), 1);
    FluxArray3D& F_Zy = F_Y;  //In computeCTUFlux, the last read of [i,j,k] happens before [i,j,k] gets written, so we can likewise repurpose F_Y here
    computeCTUFlux_Z(_zL, _zR, F_Y, F_Zy, dt_dz, (0.5*dt_dy),1);

    //Update the X half states based on the YZ corner fluxes
    correctState(_xL, _xR, F_Yz, (0.5*dt_dy), 1);
    correctState(_xL, _xR, F_Zy, (0.5*dt_dz), 2);
    //Doing this early means we can reuse the F_Yz and F_Zy grids for F_Yx and F_Zx
    
    FluxArray3D& F_Yx = F_Yz; //Reuse existing grid that has already served its purpose
    computeCTUFlux_Y(_yL, _yR, F_X, F_Yx, dt_dy, (0.5*dt_dx), 0);
    FluxArray3D& F_Zx = F_Zy; //Reuse existing grid that has already served its purpose
    computeCTUFlux_Z(_zL, _zR, F_X, F_Zx, dt_dz, (0.5*dt_dx),0);

    //Update the Y half states based on the XZ corner fluxes
    correctState(_yL, _yR, F_Xz, (0.5*dt_dx), 0);
    correctState(_yL, _yR, F_Zx, (0.5*dt_dz), 2);
    //Update the Z half states based on the XY corner fluxes
    correctState(_zL, _zR, F_Xy, (0.5*dt_dx), 0);
    correctState(_zL, _zR, F_Yx, (0.5*dt_dy), 1);
}

//MARK: CTU MHD 6-Solve (3D)
#ifdef MHD
void ctu_sweep_MHD(FluidArray3D& _xL, FluidArray3D& _xR, FluidArray3D& _yL, FluidArray3D& _yR, FluidArray3D& _zL, FluidArray3D& _zR, const MagneticArray3D &A, const MagneticArray3D &B, const FluidArray3D& w, MagneticArray3D& E, double dt, double dx, double dy, double dz){
    const int nx = _xL.getSizeX(), ny = _xL.getSizeY(), nz = _xL.getSizeZ(), g = _xL.getGhosts();
    const double dt_dx = dt/dx, dt_dy = dt/dy, dt_dz = dt/dz;

    //Preliminary Fluxes
        auto __fluxes = DRAGONWING::requestFluxArrays(3, nx, ny, nz, g);
    FluxArray3D& F_X = *__fluxes[0];
    FluxArray3D& F_Y = *__fluxes[1];
    FluxArray3D& F_Z = *__fluxes[2];
    computeFlux_X(_xL, _xR, F_X, -2, nx+2, -2, ny+2, -2, nz+2, dt_dx);
    computeFlux_Y(_yL, _yR, F_Y, -2, nx+2, -2, ny+2, -2, nz+2, dt_dy);
    computeFlux_Z(_zL, _zR, F_Z, -2, nx+2, -2, ny+2, -2, nz+2, dt_dz);
    
    //Preliminary CT Update
        auto __mags = DRAGONWING::requestVec3Arrays(3, nx+1, ny+1, nz+1, g);
    MagneticArray3D &_A = *__mags[0], &E0 = *__mags[1], &Bhalf = *__mags[2];
    //Compute E
    CT::computeElectric(E0, F_X, F_Y, F_Z,1);
    CT::bodyElectric(w, E,1); //Use Ehalf for Eref since we don't need it yet
    CT::upwindElectric(E0, F_X, F_Y, F_Z, E,1);
    //Compute A and B
    _A.clone(A);
    CT::updatePotential(_A, E0, dt/2,1);
    CT::computeFaceFields(_A, Bhalf, dx, dy, dz);
    //Construct Ehalf
        auto __whalf = DRAGONWING::requestPrimitiveArrays(1, nx, ny, nz, g);
    auto& whalf = *__whalf[0];
    applyFluxes(w, whalf, F_X, F_Y, F_Z, 0.5*dt_dx, 0.5*dt_dy, 0.5*dt_dz, 1);
    CT::computeBodyFields(Bhalf, whalf);
    CT::bodyElectric(whalf, E);
        __whalf.release();
    
    //CTU corrections (Fluid components)
    correctState(_xL, _xR, F_Y, 0.5*dt_dy, 1);
    correctState(_xL, _xR, F_Z, 0.5*dt_dz, 2);
    correctState(_yL, _yR, F_X, 0.5*dt_dx, 0);
    correctState(_yL, _yR, F_Z, 0.5*dt_dz, 2);
    correctState(_zL, _zR, F_X, 0.5*dt_dx, 0);
    correctState(_zL, _zR, F_Y, 0.5*dt_dy, 1);
        __fluxes.release();
    
    
    //Correct Magnetic fields
    for(int i=-1; i<nx+1; i++){
        for(int j=-1; j<ny+1; j++){
            for(int k=-1; k<nz+1; k++){
                double dBx = (B[i+1,j,k].x - B[i,j,k].x) * dt_dx;
                double dBy = (B[i,j+1,k].y - B[i,j,k].y) * dt_dy;
                double dBz = (B[i,j,k+1].z - B[i,j,k].z) * dt_dz;
                

                ConservativeState uL, uR;
                
                auto& W = w[i,j,k];
                
                
                auto corr = ConservativeState();
                double lim_y = TVD::minmod(-dBz, dBx);
                double lim_z = TVD::minmod(-dBy, dBx);
                corr.p = 0.5 * _1_4pi * dBx * W.B;
                corr.E = 0.5 * _1_4pi * (W.B.y * W.v.y * lim_y + W.B.z * W.v.z * lim_z);
                corr.B.y = -0.25 * dt_dz * (E0[i,j,k+1].x - E0[i,j,k].x + E0[i,j+1,k+1].x - E0[i,j+1,k].x) + 0.5 * W.v.y * lim_y;
                corr.B.z = 0.25 * dt_dy * (E0[i,j+1,k].x - E0[i,j,k].x + E0[i,j+1,k+1].x - E0[i,j,k+1].x) + 0.5 * W.v.z * lim_z;
                uL = _xL[i,j,k] + corr;
                uR = _xR[i,j,k] + corr;
                uL.B.x = Bhalf[i,j,k].x;
                uR.B.x = Bhalf[i+1,j,k].x;
                _xL[i,j,k] = uL; _xR[i,j,k] = uR;

                double lim_x = TVD::minmod(-dBz, dBy);
                lim_z = TVD::minmod(-dBx, dBy);
                corr = ConservativeState();
                corr.p = 0.5* (_1_4pi * dBy) * W.B;
                corr.E = 0.5 * _1_4pi * (W.B.x * W.v.x * lim_x + W.B.z * W.v.z * lim_z);
                corr.B.x = 0.25 * dt_dz * (E0[i,j,k+1].y - E0[i,j,k].y + E0[i+1,j,k+1].y - E0[i+1,j,k].y) + 0.5 * W.v.x * lim_x;
                corr.B.z = -0.25 * dt_dx * (E0[i+1,j,k].y - E0[i,j,k].y + E0[i+1,j,k+1].y - E0[i,j,k+1].y) + 0.5 * W.v.z * lim_z;
                uL = _yL[i,j,k] + corr;
                uR = _yR[i,j,k] + corr;
                uL.B.y = Bhalf[i,j,k].y;
                uR.B.y = Bhalf[i,j+1,k].y;
                _yL[i,j,k] = uL; _yR[i,j,k] = uR;
                
               
                lim_x = TVD::minmod(-dBy, dBz);
                lim_y = TVD::minmod(-dBx, dBz);
                corr = ConservativeState();
                corr.p = 0.5 * (_1_4pi * dBz) * W.B;
                corr.E = 0.5 * _1_4pi * (W.B.x * W.v.x * lim_x + W.B.y * W.v.y * lim_y);
                corr.B.x =  -0.25 * dt_dy * (E0[i,j+1,k].z - E0[i,j,k].z + E0[i+1,j+1,k].z - E0[i+1,j,k].z) + 0.5 * W.v.x * lim_x;
                corr.B.y =  0.25 * dt_dx * (E0[i+1,j,k].z - E0[i,j,k].z + E0[i+1,j+1,k].z - E0[i,j+1,k].z) + 0.5 * W.v.y * lim_y;
                uL = _zL[i,j,k] + corr;
                uR = _zR[i,j,k] + corr;
                uL.B.z = Bhalf[i,j,k].z;
                uR.B.z = Bhalf[i,j,k+1].z;
                _zL[i,j,k] = uL; _zR[i,j,k] = uR;
            }
        }
    }
}
#endif

//MARK: CTU MHD 4-Solve (2D)
#ifdef MHD
void ctu_sweep_MHD(FluidArray2D& _xL, FluidArray2D& _xR, FluidArray2D& _yL, FluidArray2D& _yR, const MagneticArray2D& A, const MagneticArray2D &B, const FluidArray2D& w, MagneticArray2D& E, double dt, double dx, double dy){
    const int nx = _xL.getSizeX(), ny = _xL.getSizeY(), g = _xL.getGhosts();
    const double dt_dx = dt/dx, dt_dy = dt/dy;

    //Preliminary Fluxes
        auto __fluxes = DRAGONWING::requestFluxArrays(2, nx, ny, g);
    FluxArray2D& F_X = *__fluxes[0];
    FluxArray2D& F_Y = *__fluxes[1];
    computeFlux_X(_xL, _xR, F_X, -2, nx+2, -2, ny+2, dt_dx);
    computeFlux_Y(_yL, _yR, F_Y, -2, nx+2, -2, ny+2, dt_dy);
    
    //Preliminary CT Update
        auto __mags = DRAGONWING::requestVec3Arrays(3, nx+1, ny+1, g);
    MagneticArray2D &_A = *__mags[0], &E0 = *__mags[1], &Bhalf = *__mags[2];
    //Compute E
    CT::computeElectric(E0, F_X, F_Y, 1);
    CT::bodyElectric(w, E,1); //Use Ehalf for Eref since we don't need it yet
    CT::upwindElectric(E0, F_X, F_Y, E,1);
    //Compute A and B
    _A.clone(A);
    CT::updatePotential(_A, E0, dt/2,1);
    CT::computeFaceFields(_A, Bhalf, dx, dy);
    //Construct Ehalf
        auto __whalf = DRAGONWING::requestPrimitiveArrays(1, nx, ny,g);
    auto& whalf = *__whalf[0];
    applyFluxes(w, whalf, F_X, F_Y, 0.5*dt_dx, 0.5*dt_dy, 1);
    CT::computeBodyFields(Bhalf, whalf);
    CT::bodyElectric(whalf, E);
        __whalf.release();
    
    
    //CTU corrections (Fluid components)
    correctState(_xL, _xR, F_Y, 0.5*dt_dy, 1);
    correctState(_yL, _yR, F_X, 0.5*dt_dx, 0);
    //CT::copyFaceFields_X(_xL, Bhalf, _xR);
    //CT::copyFaceFields_Y(_yL, Bhalf, _yR);

    //Correct Magnetic fields
    for(int i=-1; i<nx+1; i++){
        for(int j=-1; j<ny+1; j++){
            double dBx = (B[i+1,j].x - B[i,j].x) * dt_dx;
            double dBy = (B[i,j+1].y - B[i,j].y) * dt_dy;

            ConservativeState uL, uR;
            
            auto& W = w[i,j];
            
            
            auto corr = ConservativeState();
            double lim_z = TVD::minmod(-dBy, dBx);
            corr.p = 0.5 * _1_4pi * dBx * W.B;
            corr.E = 0.5 * _1_4pi * (W.B.z * W.v.z * lim_z);
            corr.B.z = 0.5 * dt_dy * (E0[i,j+1].x - E0[i,j].x) + 0.5 * W.v.z * lim_z;
            uL = _xL[i,j] + corr;
            uR = _xR[i,j] + corr;
            uL.B.x = Bhalf[i,j].x;
            uR.B.x = Bhalf[i+1,j].x;
            _xL[i,j] = uL; _xR[i,j] = uR;

            lim_z = TVD::minmod(-dBx, dBy);
            corr = ConservativeState();
            corr.p = 0.5* (_1_4pi * dBy) * W.B;
            corr.E = 0.5 * _1_4pi * (W.B.z * W.v.z * lim_z);
            corr.B.z = -0.5 * dt_dx * (E0[i+1,j].y - E0[i,j].y) + 0.5 * W.v.z * lim_z;
            uL = _yL[i,j] + corr;
            uR = _yR[i,j] + corr;
            uL.B.y = Bhalf[i,j].y;
            uR.B.y = Bhalf[i,j+1].y;
            _yL[i,j] = uL; _yR[i,j] = uR;
        }
    }
}
#endif



//MARK: CTU Corrections


void correctState(FluidArray2D& _L, FluidArray2D& _R, const FluxArray2D& F, double dt_dL, int dim){
    const int xL = -1, xR = F.getSizeX()+1, yL = -1, yR = F.getSizeY()+1;
    //Extract direction encoding
    int isX = dim%2 == 0 ? 1 : 0;
    int isY = dim%2 == 1 ? 1 : 0;
    //Cycle
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<yR; j++){
            auto trans = (F[i+isX,j+isY] - F[i,j]) * dt_dL;
            #ifdef MHD
            trans.B = {0,0,0};
            #endif
            _L[i,j] = _L[i,j] -  trans;
            _R[i,j] = _R[i,j] -  trans;
        }
    }
}
void correctState(FluidArray3D& _L, FluidArray3D& _R, const FluxArray3D& F, double dt_dL, int dim){
    const int xL = -1, xR = F.getSizeX()+1, yL = -1, yR = F.getSizeY()+1, zL = -1, zR = F.getSizeZ()+1;
    //Extract direction encoding
    int isX = dim%3 == 0 ? 1 : 0;
    int isY = dim%3 == 1 ? 1 : 0;
    int isZ = dim%3 == 2 ? 1 : 0;
    //cycle
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<yR; j++){
            for(int k=zL; k<zR; k++){
                auto trans = (F[i+isX, j+isY, k+isZ] - F[i,j,k]) * dt_dL;
                #ifdef MHD
                trans.B = {0,0,0};
                #endif
                _L[i,j,k] = _L[i,j,k] -  trans;
                _R[i,j,k] = _R[i,j,k] -  trans;
            }
        }
    }
}

//MARK: CTU Correction -> Flux
//Compute X fluxes between Right(_R) and Left (_L) half-states, as corrected by transverse fluxes FYZ
//Equivalent to correctState -> computeFlux_X but without needing intermediate arrays
void computeCTUFlux_X(const FluidArray3D& _L, const FluidArray3D& _R, const FluxArray3D& FYZ, FluxArray3D& F,  double dt_dx, double dt_dy, int dim){
    const int xL = -1, xR = F.getSizeX()+1, yL = -1, yR = F.getSizeY()+1, zL = -1, zR = F.getSizeZ()+1;
    //Extract direction encoding
    int isY = dim%3 == 1 ? 1 : 0;
    int isZ = dim%3 == 2 ? 1 : 0;
    //cycle
    for(int j=yL; j<yR; j++){
        for(int k=zL; k<zR; k++){
            //Compute Left half-state on the first cell
            PrimitiveState  __R, __L =  _R[xL-1, j,k] - (FYZ[xL-1, j+isY, k+isZ] - FYZ[xL-1,j,k]) * dt_dy;
            for(int i=xL; i<=xR; i++){
                auto trans = (FYZ[i, j+isY, k+isZ] - FYZ[i,j,k]) * dt_dy;
                __R = _L[i,j,k] -  trans;//Compute Right half-state for this flux
                F[i,j,k] = Riemann(__L, __R).flux_X(dt_dx);
                __L = _R[i,j,k] -  trans;//Compute Left half-state for the next flux
            }
        }
    }
}
//Compute Y fluxes between Right(_R) and Left (_L) half-states, as corrected by transverse fluxes FXZ
//Equivalent to correctState -> computeFlux_Y but without needing intermediate arrays
void computeCTUFlux_Y(const FluidArray3D& _L, const FluidArray3D& _R, const FluxArray3D& FXZ, FluxArray3D& F, double dt_dy, double dt_dz, int dim){
    const int xL = -1, xR = F.getSizeX()+1, yL = -1, yR = F.getSizeY()+1, zL = -1, zR = F.getSizeZ()+1;

    //Extract direction encoding
    int isX = dim%3 == 0 ? 1 : 0;
    int isZ = dim%3 == 2 ? 1 : 0;
    //cycle
    for(int i=xL; i<xR; i++){
        for(int k=zL; k<zR; k++){
            //Compute Left half-state on the first cell
            PrimitiveState  __R, __L =  _R[i,yL-1,k] - (FXZ[i+isX,yL-1, k+isZ] - FXZ[i,yL-1,k]) * dt_dz;
            for(int j=yL; j<=yR; j++){
                auto trans = (FXZ[i+isX, j, k+isZ] - FXZ[i,j,k]) * dt_dz;
                __R = _L[i,j,k] -  trans;//Compute Right half-state for this flux
                F[i,j,k] = Riemann(__L, __R).flux_Y(dt_dy);
                __L = _R[i,j,k] -  trans;//Compute Left half-state for the next flux
            }
        }
    }
}
//Compute Z fluxes between Right(_R) and Left (_L) half-states, as corrected by transverse fluxes FXY
//Equivalent to correctState -> computeFlux_Z but without needing intermediate arrays
void computeCTUFlux_Z(const FluidArray3D& _L, const FluidArray3D& _R, const FluxArray3D& FXY, FluxArray3D& F, double dt_dz, double dt_dy, int dim){
    const int xL = -1, xR = F.getSizeX()+1, yL = -1, yR = F.getSizeY()+1, zL = -1, zR = F.getSizeZ()+1;
    //Extract direction encoding
    int isX = dim%3 == 0 ? 1 : 0;
    int isY = dim%3 == 1 ? 1 : 0;
    //cycle
    for(int i=xL; i<xR; i++){
        for(int j=yL; j<yR; j++){
            //Compute Left half-state on the first cell
            PrimitiveState  __R, __L =  _R[i,j,zL-1] - (FXY[i+isX,j+isY,zL-1] - FXY[i,j,zL-1]) * dt_dy;
            for(int k=zL; k<=zR; k++){
                auto trans = (FXY[i+isX, j+isY, k] - FXY[i,j,k]) * dt_dy;
                __R = _L[i,j,k] -  trans;//Compute Right half-state for this flux
                F[i,j,k] = Riemann(__L, __R).flux_Z(dt_dz);
                __L = _R[i,j,k] -  trans;//Compute Left half-state for the next flux
            }
        }
    }
}
