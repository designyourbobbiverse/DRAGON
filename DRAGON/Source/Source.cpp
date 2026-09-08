//
//  Source.cpp
//  DRAGON/Source
//
//  Created by Bobbie Markwick on 04/09/2026.
//

#include "Source.hpp"
#include "Grid.hpp" //Grid::source_step


#include "Config.h"     //For integration scheme
#include "Constants.h"  //For gamma

using namespace DRAGON;
using namespace Source;


//MARK: Grid Sweep
void Grid::source_step(double dt){
    if(sources.count() > 0) sources.advance(*this, dt);
}
void SourceTerm::advance(Grid& grid, double dt){
    Grid3D* grid3D = dynamic_cast<Grid3D*>(&grid);
    if (grid3D) {
        advance(*grid3D, dt);
        return;
    }
    Grid2D* grid2D = dynamic_cast<Grid2D*>(&grid);
    if (grid2D) {
        advance(*grid2D, dt);
        return;
    }
    Grid1D* grid1D = dynamic_cast<Grid1D*>(&grid);
    if (grid1D) {
        advance(*grid1D, dt);
        return;
    }
}

void SourceTerm::advance(Grid1D& w, double dt){
    const int nx = w.getSize();
    for(int i=0; i<nx; i++){
        w[i] += integrate(dt, w[i]);
    }
}
void SourceTerm::advance(Grid2D& w, double dt){
     const int nx = w.getSizeX(), ny = w.getSizeY();
     for(int i=0; i<nx; i++){
         for(int j=0; j<ny; j++){
             w[i,j] +=  integrate(dt, w[i,j]);
         }
     }
 }
void SourceTerm::advance(Grid3D& w, double dt){
     const int nx = w.getSizeX(), ny = w.getSizeY(), nz = w.getSizeZ();
     for(int i=0; i<nx; i++){
         for(int j=0; j<ny; j++){
             for(int k=0; k<nz; k++){
                 w[i,j,k] +=  integrate(dt, w[i,j,k]);
             }
         }
     }
 }


//MARK: Integration
ConservativeState SourceTerm::integrate(double dt, const PrimitiveState& w0, double t0){
    #if SRC_SPLIT_INTEGRATION == CHOOSE_RUNTIME
    switch (Config::src_integration_choice) {
    case RK2: return rk2(dt, w0, t0);
    case RK4: return rk4(dt, w0, t0);
    default: return rk2(dt, w0, t0);
    }
    #elif SRC_SPLIT_INTEGRATION == RK2
    return rk2(dt, w0, t0);
    #elif SRC_SPLIT_INTEGRATION == RK4
    return rk4(dt, w0, t0);
    #endif
}
ConservativeState SourceTerm::rk2(double dt, const PrimitiveState& w0, double t0){
    ConservativeState k1 = source_density(w0, t0);
    ConservativeState k2 = source_density(w0 + k1*dt, t0 + dt);
    return (k1 + k2) * dt / 2.0;
}
ConservativeState SourceTerm::rk4(double dt, const PrimitiveState& w0, double t0){
    ConservativeState k1 = source_density(w0, t0);
    ConservativeState k2 = source_density(w0 + 0.5*k1*dt, t0 + 0.5*dt);
    ConservativeState k3 = source_density(w0 + 0.5*k2*dt, t0 + 0.5*dt);
    ConservativeState k4 = source_density(w0 + k3*dt, t0 + dt);
    return (k1 + 2*k2 + 2*k3 + k4) * dt / 6.0;
}



//MARK: Energy/Force/Mass Types
ConservativeState EnergySource::source_density(const PrimitiveState &w, double t){
    ConservativeState S{};
    S.E = energy(w,t);
    return S;
}
ConservativeState MomentumSource::source_density(const PrimitiveState &w, double t){
    vec3 f = force(w,t);

    ConservativeState S{};
    S.mom = f;
    S.E = w.v * f;
    return S;
}
ConservativeState MassSource::source_density(const PrimitiveState &w, double t){
    double S_rho = density(w,t);
    vec3 v = velocity(w, t);
    
    ConservativeState S{};
    S.rho = S_rho;
    S.mom = S_rho * v;
    S.E = S_rho * (thermal_energy(w, t) +  0.5*v*v);
    return S;
}
vec3 MassSource::velocity(const PrimitiveState& w, double t){
    return w.v;
}
double MassSource::thermal_energy(const PrimitiveState &w, double t){
    return w.p/((_gamma-1) * w.rho);
}
