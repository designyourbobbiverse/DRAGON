//
//  Unsplit.hpp
//  DRAGON/Hydro/Godunov
//
//  Created by Bobbie Markwick on 13/07/2026.
//

#ifndef Unsplit_h
#define Unsplit_h

#include "Grid.hpp"
#include "CT.hpp"


//2D
void computeFlux_X(const FluidArray2D& _L, const FluidArray2D& _R, FluxArray2D& F, int xL, int xR, int yL, int yR, double dt_dx);
void computeFlux_Y(const FluidArray2D& _L, const FluidArray2D& _R, FluxArray2D& F, int xL, int xR, int yL, int yR, double dt_dy);
void applyFluxes(const FluidArray2D& w, FluidArray2D& _w, const FluxArray2D& F_X, const FluxArray2D& F_Y,  double dt_dx, double dt_dy,  int g = 0);
void computeHalfStates_X(FluidArray2D& _L, const Grid2D& _W,  FluidArray2D& _R, double dt);
void computeHalfStates_Y( FluidArray2D& _L,const Grid2D& _W, FluidArray2D& _R, double dt);
//2D CTU
#ifdef MHD
void ctu_sweep_MHD(FluidArray2D& _xL, FluidArray2D& _xR, FluidArray2D& _yL, FluidArray2D& _yR, const MagneticArray2D& A, const MagneticArray2D &B, const FluidArray2D& w, MagneticArray2D& Ehalf, double dt, double dx, double dy);
#else
void ctu_sweep_hydro(FluidArray2D& _xL, FluidArray2D& _xR, FluidArray2D& _yL, FluidArray2D& _yR,  double dt_dx, double dt_dy);
#endif
void correctState(FluidArray2D& _L, FluidArray2D& _R, const FluxArray2D& F, double dt_dL, int dim);



//3D
void computeFlux_X(const FluidArray3D& _L, const FluidArray3D& _R, FluxArray3D& F, int xL, int xR, int yL, int yR, int zL, int zR, double dt_dx);
void computeFlux_Y(const FluidArray3D& _L, const FluidArray3D& _R, FluxArray3D& F, int xL, int xR, int yL, int yR, int zL, int zR, double dt_dy);
void computeFlux_Z(const FluidArray3D& _L, const FluidArray3D& _R, FluxArray3D& F, int xL, int xR, int yL, int yR, int zL, int zR, double dt_dz);
void applyFluxes(const FluidArray3D& w, FluidArray3D& _w, const FluxArray3D& F_X, const FluxArray3D& F_Y, const FluxArray3D& F_Z, double dt_dx, double dt_dy, double dt_dz, int g=0);
void computeHalfStates_X(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, double dt);
void computeHalfStates_Y(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, double dt);
void computeHalfStates_Z(FluidArray3D& _L, const Grid3D& _W, FluidArray3D& _R, double dt);
//3DCTU
#ifdef MHD
void ctu_sweep_MHD(FluidArray3D& _xL, FluidArray3D& _xR, FluidArray3D& _yL, FluidArray3D& _yR, FluidArray3D& _zL, FluidArray3D& _zR, const MagneticArray3D& A, const MagneticArray3D &B, const FluidArray3D& w, MagneticArray3D& Ehalf, double dt, double dx, double dy, double dz);
void correctState(FluidArray3D& _L, FluidArray3D& _R, const FluxArray3D& F, double dt_dL, int dim);
#else
void ctu_sweep_hydro(FluidArray3D& _xL, FluidArray3D& _xR, FluidArray3D& _yL, FluidArray3D& _yR, FluidArray3D& _zL, FluidArray3D& _zR, double dt_dx, double dt_dy, double dt_dz);
void correctState(FluidArray3D& _L, FluidArray3D& _R, const FluxArray3D& F, double dt_dL, int state_dim, int flux_dim);
void computeCTUFlux_X(const FluidArray3D& _L, const FluidArray3D& _R, const FluxArray3D& Ftrans, FluxArray3D& F, double dt_dx, double dt_dy, int dim);
void computeCTUFlux_Y(const FluidArray3D& _L, const FluidArray3D& _R, const FluxArray3D& Ftrans, FluxArray3D& F,  double dt_dy, double dt_dz, int dim);
void computeCTUFlux_Z(const FluidArray3D& _L, const FluidArray3D& _R, const FluxArray3D& Ftrans, FluxArray3D& F, double dt_dz, double dt_dy, int dim);
#endif

#endif /* Unsplit_h */
