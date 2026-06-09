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
    
//MARK: Other Approximate Solvers
    ///Star Region Approximation Solvers (Toro Ch 9)
    RiemannSolution PVRS();//Primative Variable Riemann Solver
    RiemannSolution TRRS();//Two-Rarefaction Riemann Solver
    RiemannSolution TSRS(); //Two-Shock Riemann Solver
    
//MARK: Adaptive Solvers
    //Adaptive Iterative Riemann Solver: PVRS if within limits, otherwise pivots to iterative
    RiemannSolution PVRS_Iter();
    //Adaptive Noniterative Riemann Solver: PVRS if within limits, otherwise pivots to TRRS/TSRS
    RiemannSolution PVRS_TXRS();
    
private:
    //Internal Implementations
    RiemannSolution PVRS(double aL, double aR, double p_pvrs);
    RiemannSolution TRRS(double aL, double aR);
    RiemannSolution TSRS(double aL, double aR, double pGuess);
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
