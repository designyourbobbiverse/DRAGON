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
    PrimativeState L, R;

    ///Exact Riemann Solvers
    RiemannSolution exact();
    RiemannSolution exact(double pGuess);

    double exact_StarP();
    double exact_StarP(double pGuess);
    double exact_StarV(double pStar);
    double exact_StarRho(PrimativeState L, double pStar);

    ///Flux Approximation Solvers

    //HLL Solvers
    ConservativeState HLL();
    ConservativeState HLL(double SL, double SR);
    ConservativeState HLLC();
    ConservativeState HLLC(double SL, double SR);

    //Roe Solver
    ConservativeState Roe();

    ///Star Region Approximation Solvers (Toro Ch 9)
    //Primative Variable Riemann Solver
    RiemannSolution PVRS();
    //Two-Rarefaction Riemann Solver
    RiemannSolution TRRS();
    //Two-Shock Riemann Solver
    RiemannSolution TSRS();
    
    ///Adaptive Solvers
    //Two-Rarefaction if conditions met, otherwise pivots to iterative
    RiemannSolution TRRS_Iter();
    RiemannSolution TRRS_Iter(double pGuess);
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
    

struct RiemannSolution{
    PrimativeState wL, wR; //Left and Right Initial States
    PrimativeState sL, sR; //Left and Right Star Regions
    
    //Returns the state along the line x/t
    PrimativeState sample(double x_t);
    //Returns the flux through x/t
    ConservativeState flux(double x_t);
    //Flux through x=0
    ConservativeState flux();


    //Creates a Riemann Solution from a pair of initial values
    //Copies the y and z velocities, but does not solve for rho/vx/p
    RiemannSolution(Riemann problem);
 
private:
    //Swap left and right states, negate all vx.
    void mirror();
};

    
    
#endif /* Solvers_hpp */
