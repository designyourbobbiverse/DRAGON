//
//  Riemann.hpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#ifndef Riemann_hpp
#define Riemann_hpp

#include <stdio.h>
#include "FluidElement.hpp"
#include "Config.h"


/*
 Example API usage.
 To calculate the exact flux between state wL on the left and state wR on the right:
 Riemann(wL,wR).exact().flux()
 
 Some solvers skip the full solution and directly return the flux. For example
 Riemann(wL,wR).HLLC()
 Riemann(wL,wR).HLLC(SL,SR)
 Riemann(wL,wR).Roe()
  
 */

struct Riemann;
struct RiemannSolution;

struct Riemann{
    PrimitiveState L, R;
    Riemann(PrimitiveState L, PrimitiveState R);

    
    //Computes the flux using whatever method was selected in Constants.h
    ConservativeState flux();
    ConservativeState flux_X();
    ConservativeState flux_Y();
    ConservativeState flux_Z();

    
//MARK: Exact Riemann Solvers
    RiemannSolution exact();
    RiemannSolution exact(double pGuess);
    
    static double f(double p, PrimitiveState w);
    double exact_StarP();
    double exact_StarP(double pGuess);
    double exact_StarV(double pStar);
    double exact_StarRho(PrimitiveState L, double pStar);

//MARK: HLL/HLLC
    ConservativeState HLL();
    ConservativeState HLL(double SL, double SR);
    ConservativeState HLLC();
    ConservativeState HLLC(double SL, double SR);

//MARK: Roe
    ConservativeState Roe();

    
private:
    //Internal Implementations
    RiemannSolution TRRS();//Two-Rarefaction Riemann Solver
    RiemannSolution TRRS(double aL, double aR);
};
    
//MARK: Riemann Solution State
struct RiemannSolution{
    PrimitiveState wL, wR; //Left and Right Initial States
    PrimitiveState sL, sR; //Left and Right Star Regions
    
    //Creates a Riemann Solution from a pair of initial values
    //Copies the y and z velocities, but does not solve for rho/vx/p
    RiemannSolution(Riemann problem);
    
//MARK: Solution Sampling
    //Returns the state along the line x/t
    PrimitiveState sample(double x_t);
    //Returns the flux through x/t
    ConservativeState flux(double x_t);
    //Flux through x=0
    ConservativeState flux();
 
private:
    //Swap left and right states, negate all vx.
    void mirror();
};




    
    
#endif /* Solvers_hpp */
