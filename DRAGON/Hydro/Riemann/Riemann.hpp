//
//  Riemann.hpp
//  DRAGON/Hydro/Riemann
//
//  Created by Bobbie Markwick on 28/04/2026.
//  Based mostly on Toro (2009). https://doi.org/10.1007/b79761
//

#ifndef Riemann_hpp
#define Riemann_hpp

#include "FluidElement.hpp"
#include <cmath>
#include "Config.h"


/*
 Example API usage.
 To calculate the exact flux between state wL on the left and state wR on the right:
 Riemann(wL,wR).exact().flux()
 
 Some solvers skip the full solution and directly return the flux. For example
 Riemann(wL,wR).HLLC()
 Riemann(wL,wR).HLLC(SL,SR)
 Riemann(wL,wR).Roe()
 
 To use the pre-selected solver, use
 Riemann(wL,wR).flux()
 or
 Riemann(wL,wR).flux(dt/dx) //does a safety check, pivot to exact on failure
 */

struct Riemann;
struct RiemannSolution;

struct Riemann{
    PrimitiveState L, R;
    Riemann(PrimitiveState L, PrimitiveState R);
    
    //Computes the flux using whatever method was selected in Config.h
    //dt_dx is used in verifying physicality of solution (see RIEMANN_VERIFY_FALLBACK in Config.h)
    ConservativeState flux(double dt_dx = 0);
    ConservativeState flux_X(double dt_dx = 0);
    ConservativeState flux_Y(double dt_dy = 0);
    ConservativeState flux_Z(double dt_dz = 0);

    
//MARK: Exact Riemann Solver
    RiemannSolution exact();
    RiemannSolution exact(double pGuess);
    
    static double f(double p, PrimitiveState w);
    double exact_StarP();
    double exact_StarP(double pGuess);
    double exact_StarV(double pStar);
    double exact_StarRho(PrimitiveState L, double pStar);

//MARK: Approximate Solvers
    ConservativeState HLL(); // Harten, Lax, & van Leer (1983). https://doi.org/10.1137/1025002
    ConservativeState HLL(double SL, double SR);
#ifdef HYDRO_AVAILABLE
    ConservativeState HLLC(); // Toro, Spruce, & Speares (1994). https://doi.org/10.1007/BF01414629
    ConservativeState HLLC(double SL, double SR);
#endif
#ifdef MHD
    ConservativeState HLLD(); // Miyoshi and Kusano (2005). https://doi.org/10.1016/j.jcp.2005.02.017
#endif
    ConservativeState HLLE(); // Einfeldt (1988). https://doi.org/10.1137/0725021
#ifdef HYDRO_AVAILABLE
    ConservativeState Roe(); // Roe (1981). https://doi.org/10.1016/0021-9991(81)90128-5
#endif
    
private:
    RiemannSolution TRRS();//Two-Rarefaction Riemann Solver
    
    void verify_and_fallback(ConservativeState& flux, double dt_dx);
};
    
//MARK: Riemann Solution State
struct RiemannSolution{
    PrimitiveState wL, wR; //Left and Right Initial States
    PrimitiveState sL, sR; //Left and Right Star Regions
    
    //Creates a Riemann Solution from a pair of initial values
    //Copies the y and z velocities, but does not solve for rho/v.x/p
    RiemannSolution(Riemann problem);
    
//MARK: Solution Sampling
    //Returns the state along the line x/t
    PrimitiveState sample(double x_t);
    //Returns the flux through x/t
    ConservativeState flux(double x_t);
    //Flux through x=0
    ConservativeState flux();
 
private:
    //Swap left and right states, negate all v.x.
    void mirror();
};




    
    
#endif /* Solvers_hpp */
