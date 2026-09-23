//
//  CTU.cpp
//  DRAGON/Hydro/Godunov
//
//  Created by Bobbie Markwick on 10/06/2026.
//  Implementation based on  Gardiner and Stone (2008) https://arxiv.org/abs/0712.2634
//


#include "Hydro/Grid.hpp"
#include "Godunov.hpp"

#include "MHD/CT.hpp"   //For MHD
#include "Config.h"     //For MHD vs Hydro

#include "Hydro/Riemann/Riemann.hpp"    //For Riemann Solvers
#include "DragonWing.hpp"               //For Memory Management
#include "Hydro/Reconstruction/TVD.hpp" //For MINMOD
#include "Constants.h"                  //For _1_8pi
using namespace DRAGON;

#ifdef CTU
#ifndef MUSCL_Hancock
#error CTU requires MUSCL reconstruction. Please enable MUSCL_Hancock in Config.h
#endif
#endif

#ifdef MHD
//MARK: CTU MHD 6-Solve (3D)
void Godunov::ctu_sweep_MHD(FluidArray3D& _xL, FluidArray3D& _xR, FluidArray3D& _yL, FluidArray3D& _yR, FluidArray3D& _zL, FluidArray3D& _zR, const FluxArray3D& F_X, const FluxArray3D& F_Y, const FluxArray3D& F_Z, const MagneticArray3D& Ehalf, const MagneticArray3D &Bhalf,   const FluidArray3D& w0, const MagneticArray3D &B0, double dt_dx, double dt_dy, double dt_dz){
    const int nx = _xL.getSizeX(), ny = _xL.getSizeY(), nz = _xL.getSizeZ();
    
    //CTU corrections
    for (int i=-1; i<nx+1; i++) {
        for (int j=-1; j<ny+1; j++) {
            for (int k=-1; k<nz+1; k++) {
                double dBx = (B0[i+1,j,k].x - B0[i,j,k].x) * dt_dx;
                double dBy = (B0[i,j+1,k].y - B0[i,j,k].y) * dt_dy;
                double dBz = (B0[i,j,k+1].z - B0[i,j,k].z) * dt_dz;
                

                ConservativeState uL, uR;
                
                auto& W = w0[i,j,k];
                
                //Fluxes for the hydro components
                auto corr =  0.5 * ((F_Y[i,j,k] - F_Y[i,j+1,k]) * dt_dy + (F_Z[i,j,k] - F_Z[i,j,k+1]) * dt_dz);
                corr.B = {0,0,0};
                //MHD Source Terms (see GS08 eqs 51-59)
                double lim_y = TVD::minmod(-dBz, dBx);
                double lim_z = TVD::minmod(-dBy, dBx);
                //corr.mom += _1_8pi * dBx * W.B; //Our MUSCL implemenation makes the momentum correction (GS08-51) unnecessary (& wrong)
                corr.E += _1_8pi * (W.B.y * W.v.y * lim_y + W.B.z * W.v.z * lim_z); //But we do need the energy correction (GS08-52)
                corr.B.y = -0.25 * dt_dz * (Ehalf[i,j,k+1].x - Ehalf[i,j,k].x + Ehalf[i,j+1,k+1].x - Ehalf[i,j+1,k].x) + 0.5 * W.v.y * lim_y;
                corr.B.z = 0.25 * dt_dy * (Ehalf[i,j+1,k].x - Ehalf[i,j,k].x + Ehalf[i,j+1,k+1].x - Ehalf[i,j,k+1].x) + 0.5 * W.v.z * lim_z;
                //Apply fluxes, copy normal B's at Faces
                uL = _xL[i,j,k] + corr;
                uR = _xR[i,j,k] + corr;
                uL.B.x = Bhalf[i,j,k].x;
                uR.B.x = Bhalf[i+1,j,k].x;
                _xL[i,j,k] = uL; _xR[i,j,k] = uR;

                //Fluxes for the hydro components
                corr =  0.5 * ((F_X[i,j,k] - F_X[i+1,j,k]) * dt_dx + (F_Z[i,j,k] - F_Z[i,j,k+1]) * dt_dz);
                corr.B = {0,0,0};
                //MHD Source Terms (see GS08 eqs 51-59)
                double lim_x = TVD::minmod(-dBz, dBy);
                lim_z = TVD::minmod(-dBx, dBy);
                //corr.mom += _1_8pi * dBy * W.B; //Our MUSCL implemenation makes the momentum correction (GS08-51) unnecessary (& wrong)
                corr.E += _1_8pi * (W.B.x * W.v.x * lim_x + W.B.z * W.v.z * lim_z); //But we do need the energy correction (GS08-52)
                corr.B.x = 0.25 * dt_dz * (Ehalf[i,j,k+1].y - Ehalf[i,j,k].y + Ehalf[i+1,j,k+1].y - Ehalf[i+1,j,k].y) + 0.5 * W.v.x * lim_x;
                corr.B.z = -0.25 * dt_dx * (Ehalf[i+1,j,k].y - Ehalf[i,j,k].y + Ehalf[i+1,j,k+1].y - Ehalf[i,j,k+1].y) + 0.5 * W.v.z * lim_z;
                //Apply fluxes, copy normal B's at Faces
                uL = _yL[i,j,k] + corr;
                uR = _yR[i,j,k] + corr;
                uL.B.y = Bhalf[i,j,k].y;
                uR.B.y = Bhalf[i,j+1,k].y;
                _yL[i,j,k] = uL; _yR[i,j,k] = uR;
                
                //Fluxes
                corr =  0.5 * ((F_X[i,j,k] - F_X[i+1,j,k]) * dt_dx + (F_Y[i,j,k] - F_Y[i,j+1,k]) * dt_dy);
                corr.B = {0,0,0};
                //MHD Source Terms (see GS08 eqs 51-59)
                lim_x = TVD::minmod(-dBy, dBz);
                lim_y = TVD::minmod(-dBx, dBz);
                //corr.mom += _1_8pi * dBz * W.B; //Our MUSCL implemenation makes the momentum correction (GS08-51) unnecessary (& wrong)
                corr.E += _1_8pi * (W.B.x * W.v.x * lim_x + W.B.y * W.v.y * lim_y); //But we do need the energy correction (GS08-52)
                corr.B.x =  -0.25 * dt_dy * (Ehalf[i,j+1,k].z - Ehalf[i,j,k].z + Ehalf[i+1,j+1,k].z - Ehalf[i+1,j,k].z) + 0.5 * W.v.x * lim_x;
                corr.B.y =  0.25 * dt_dx * (Ehalf[i+1,j,k].z - Ehalf[i,j,k].z + Ehalf[i+1,j+1,k].z - Ehalf[i,j+1,k].z) + 0.5 * W.v.y * lim_y;
                //Apply fluxes, copy normal B's at Faces
                uL = _zL[i,j,k] + corr;
                uR = _zR[i,j,k] + corr;
                uL.B.z = Bhalf[i,j,k].z;
                uR.B.z = Bhalf[i,j,k+1].z;
                _zL[i,j,k] = uL; _zR[i,j,k] = uR;
            }
        }
    }
}

//MARK: CTU MHD 4-Solve (2D)
void Godunov::ctu_sweep_MHD(FluidArray2D& _xL, FluidArray2D& _xR, FluidArray2D& _yL, FluidArray2D& _yR, const FluxArray2D& F_X, const FluxArray2D& F_Y, const MagneticArray2D& Ehalf, const MagneticArray2D &Bhalf,   const FluidArray2D& w0, const MagneticArray2D &B0, double dt_dx, double dt_dy){
    const int nx = _xL.getSizeX(), ny = _xL.getSizeY();
    
    //CTU Corrections
    for (int i=-1; i<nx+1; i++) {
        for (int j=-1; j<ny+1; j++) {
            double dBx = (B0[i+1,j].x - B0[i,j].x) * dt_dx;
            double dBy = (B0[i,j+1].y - B0[i,j].y) * dt_dy;

            ConservativeState uL, uR;
            
            auto& W = w0[i,j];
            
            //Fluxes
            auto corr =  0.5 * (F_Y[i,j] - F_Y[i,j+1]) * dt_dy;
            corr.B = {0,0,0};
            //MHD Source Terms
            double lim_z = TVD::minmod(-dBy, dBx);
            //corr.mom +=  _1_8pi * dBx * W.B;
            corr.E += _1_8pi * (W.B.z * W.v.z * lim_z);
            corr.B.z = 0.5 * dt_dy * (Ehalf[i,j+1].x - Ehalf[i,j].x) + 0.5 * W.v.z * lim_z;
            //Apply Fluxes & Face Fields
            uL = _xL[i,j] + corr;
            uR = _xR[i,j] + corr;
            uL.B.x = Bhalf[i,j].x;
            uR.B.x = Bhalf[i+1,j].x;
            _xL[i,j] = uL; _xR[i,j] = uR;

            //Fluxes
            corr =  0.5 * (F_X[i,j] - F_X[i+1,j]) * dt_dx;
            corr.B = {0,0,0};
            //MHD Source Terms
            lim_z = TVD::minmod(-dBx, dBy);
            //corr.mom +=  _1_8pi * dBy * W.B;
            corr.E +=  _1_8pi * (W.B.z * W.v.z * lim_z);
            corr.B.z = -0.5 * dt_dx * (Ehalf[i+1,j].y - Ehalf[i,j].y) + 0.5 * W.v.z * lim_z;
            //Apply Fluxes & Face Fields
            uL = _yL[i,j] + corr;
            uR = _yR[i,j] + corr;
            uL.B.y = Bhalf[i,j].y;
            uR.B.y = Bhalf[i,j+1].y;
            _yL[i,j] = uL; _yR[i,j] = uR;
        }
    }
}

#else
//MARK: CTU 2D Hydro
void Godunov::ctu_sweep_hydro(FluidArray2D& _xL, FluidArray2D& _xR, FluidArray2D& _yL, FluidArray2D& _yR,  FluxArray2D& F_X, FluxArray2D& F_Y,  double dt_dx, double dt_dy){
    correctState(_xL, _xR, F_Y, (0.5*dt_dy), 1);
    correctState(_yL, _yR, F_X, (0.5*dt_dx), 0);
}


//MARK: CTU Hydro 12-Solve (3D)
//See Gardiner and Stone (2008) section 5.1
void Godunov::ctu_sweep_hydro(FluidArray3D& _xL, FluidArray3D& _xR, FluidArray3D& _yL, FluidArray3D& _yR, FluidArray3D& _zL, FluidArray3D& _zR, FluxArray3D& F_X, FluxArray3D& F_Y, FluxArray3D& F_Z, double dt_dx, double dt_dy, double dt_dz){
    const int nx = _xL.getSizeX(), ny = _xL.getSizeY(), nz = _xL.getSizeZ(), ghosts = _xL.getGhosts();
        
    //Compute Edge-correct the fluxes
        auto __fluxes = DRAGONWING::requestFluxArrays(2, nx, ny, nz, ghosts);
    FluxArray3D& F_Xz = *__fluxes[0];
    computeCTUFlux_X(_xL, _xR, F_Z, F_Xz, dt_dx, (dt_dz/3.0), 2);
    FluxArray3D& F_Yz = *__fluxes[1];
    computeCTUFlux_Y(_yL, _yR, F_Z, F_Yz, dt_dy, (dt_dz/3.0), 2);

    FluxArray3D& F_Xy = F_Z; //F_Z isn't used again before it gets recomputed, so we can resue it here
    computeCTUFlux_X(_xL, _xR, F_Y, F_Xy, dt_dx, (dt_dy/3.0), 1);
    FluxArray3D& F_Zy = F_Y;  //In computeCTUFlux, the last [i,j,k] read is before the [i,j,k] write, so we can likewise reuse F_Y here
    computeCTUFlux_Z(_zL, _zR, F_Y, F_Zy, dt_dz, (dt_dy/3.0),1);

    //Update the X half states based on the YZ corner fluxes
    correctState(_xL, _xR, F_Yz, (0.5*dt_dy), 0,1);
    correctState(_xL, _xR, F_Zy, (0.5*dt_dz), 0,2);
    //Doing this early means we can reuse the F_Yz and F_Zy grids for F_Yx and F_Zx
    
    FluxArray3D& F_Yx = F_Yz; //Reuse existing grid that has already served its purpose
    computeCTUFlux_Y(_yL, _yR, F_X, F_Yx, dt_dy, (dt_dx/3.0), 0);
    FluxArray3D& F_Zx = F_Zy; //Reuse existing grid that has already served its purpose
    computeCTUFlux_Z(_zL, _zR, F_X, F_Zx, dt_dz, (dt_dx/3.0),0);

    //Update the Y half states based on the XZ corner fluxes
    correctState(_yL, _yR, F_Xz, (0.5*dt_dx), 1,0);
    correctState(_yL, _yR, F_Zx, (0.5*dt_dz), 1,2);
    //Update the Z half states based on the XY corner fluxes
    correctState(_zL, _zR, F_Xy, (0.5*dt_dx), 2,0);
    correctState(_zL, _zR, F_Yx, (0.5*dt_dy), 2,1);
}

//MARK: CTU Corrections
void Godunov::correctState(FluidArray2D& _L, FluidArray2D& _R, const FluxArray2D& F, double dt_dL, int dim){
    //Extract direction encoding: 0 -> x, 1 -> y
    const int isX = dim%2 == 0 ? 1 : 0;
    const int isY = dim%2 == 1 ? 1 : 0;
    const int xL = -isY, xR = F.getSizeX()+isY, yL = -isX, yR = F.getSizeY()+isX;
    
    //Cycle
    for (int i=xL; i<xR; i++) {
        for (int j=yL; j<yR; j++) {
            auto trans = (F[i+isX,j+isY] - F[i,j]) * dt_dL;
            _L[i,j] = _L[i,j] -  trans;
            _R[i,j] = _R[i,j] -  trans;
        }
    }
}

void Godunov::correctState(FluidArray3D& _L, FluidArray3D& _R, const FluxArray3D& F, double dt_dL, int state_dim, int flux_dim){
    //Extract direction encoding: 0 -> x, 1 -> y, 2 -> z
    const int isX = flux_dim%3 == 0 ? 1 : 0;
    const int isY = flux_dim%3 == 1 ? 1 : 0;
    const int isZ = flux_dim%3 == 2 ? 1 : 0;
    const int isX_ = state_dim%3 == 0 ? 1 : 0;
    const int isY_ = state_dim%3 == 1 ? 1 : 0;
    const int isZ_ = state_dim%3 == 2 ? 1 : 0;
    
    const int xL = -isX_, xR = F.getSizeX()+isX_, yL = -isY_, yR = F.getSizeY()+isY_, zL = -isZ_, zR = F.getSizeZ()+isZ_;
    
    //cycle
    for (int i=xL; i<xR; i++) {
        for (int j=yL; j<yR; j++) {
            for (int k=zL; k<zR; k++) {
                auto trans = (F[i+isX, j+isY, k+isZ] - F[i,j,k]) * dt_dL;
                _L[i,j,k] = _L[i,j,k] -  trans;
                _R[i,j,k] = _R[i,j,k] -  trans;
            }
        }
    }
}
    
//MARK: 12-Solve Correction -> Flux
//Compute X fluxes between _R[i-1] and _L[i] half-states, as corrected by transverse fluxes FYZ
//Equivalent to correctState (with wider bounds) -> computeFlux_X but without needing intermediate arrays
void Godunov::computeCTUFlux_X(const FluidArray3D& _L, const FluidArray3D& _R, const FluxArray3D& FYZ, FluxArray3D& F,  double dt_dx, double dt_dy, int dim){
    //Extract direction encoding
    int isY = dim%3 == 1 ? 1 : 0;
    int isZ = dim%3 == 2 ? 1 : 0;
    
    const int xL = 0, xR = F.getSizeX(), yL = -1+isY, yR = F.getSizeY()+1-isY, zL = -1+isZ, zR = F.getSizeZ()+1-isZ;
    
    //cycle
    for (int j=yL; j<yR; j++) {
        for (int k=zL; k<zR; k++) {
            //Compute Left half-state on the first cell
            PrimitiveState  __R, __L =  _R[xL-1, j,k] - (FYZ[xL-1, j+isY, k+isZ] - FYZ[xL-1,j,k]) * dt_dy;
            for (int i=xL; i<=xR; i++) {
                auto trans = (FYZ[i, j+isY, k+isZ] - FYZ[i,j,k]) * dt_dy;
                __R = _L[i,j,k] -  trans;//Compute Right half-state for this flux
                F[i,j,k] = Riemann(__L, __R).flux_X(dt_dx);
                __L = _R[i,j,k] -  trans;//Compute Left half-state for the next flux
            }
        }
    }
}
//Compute Y fluxes between _R[j-1] and L[j] half-states, as corrected by transverse fluxes FXZ
//Equivalent to correctState (with wider bounds) -> computeFlux_Y but without needing intermediate arrays
void Godunov::computeCTUFlux_Y(const FluidArray3D& _L, const FluidArray3D& _R, const FluxArray3D& FXZ, FluxArray3D& F, double dt_dy, double dt_dz, int dim){
    //Extract direction encoding
    const int isX = dim%3 == 0 ? 1 : 0;
    const int isZ = dim%3 == 2 ? 1 : 0;
    const int xL = -1+isX, xR = F.getSizeX()+1-isX, yL = 0, yR = F.getSizeY(), zL = -1+isZ, zR = F.getSizeZ()+1-isZ;
    
    
    //cycle
    for (int i=xL; i<xR; i++) {
        for (int k=zL; k<zR; k++) {
            //Compute Left half-state on the first cell
            PrimitiveState  __R, __L =  _R[i,yL-1,k] - (FXZ[i+isX,yL-1, k+isZ] - FXZ[i,yL-1,k]) * dt_dz;
            for (int j=yL; j<=yR; j++) {
                auto trans = (FXZ[i+isX, j, k+isZ] - FXZ[i,j,k]) * dt_dz;
                __R = _L[i,j,k] -  trans;//Compute Right half-state for this flux
                F[i,j,k] = Riemann(__L, __R).flux_Y(dt_dy);
                __L = _R[i,j,k] -  trans;//Compute Left half-state for the next flux
            }
        }
    }
}
//Compute Z fluxes between _R[k-1] and _L[k] half-states, as corrected by transverse fluxes FXY
//Equivalent to correctState (with wider bounds) -> computeFlux_Z but without needing intermediate arrays
void Godunov::computeCTUFlux_Z(const FluidArray3D& _L, const FluidArray3D& _R, const FluxArray3D& FXY, FluxArray3D& F, double dt_dz, double dt_dy, int dim){
    //Extract direction encoding
    int isX = dim%3 == 0 ? 1 : 0;
    int isY = dim%3 == 1 ? 1 : 0;
    const int xL = -1+isX, xR = F.getSizeX()+1-isX, yL = -1+isY, yR = F.getSizeY()+1-isY, zL = 0, zR = F.getSizeZ();
    //cycle
    for (int i=xL; i<xR; i++) {
        for (int j=yL; j<yR; j++) {
            //Compute Left half-state on the first cell
            PrimitiveState  __R, __L =  _R[i,j,zL-1] - (FXY[i+isX,j+isY,zL-1] - FXY[i,j,zL-1]) * dt_dy;
            for (int k=zL; k<=zR; k++) {
                auto trans = (FXY[i+isX, j+isY, k] - FXY[i,j,k]) * dt_dy;
                __R = _L[i,j,k] -  trans;//Compute Right half-state for this flux
                F[i,j,k] = Riemann(__L, __R).flux_Z(dt_dz);
                __L = _R[i,j,k] -  trans;//Compute Left half-state for the next flux
            }
        }
    }
}
#endif
