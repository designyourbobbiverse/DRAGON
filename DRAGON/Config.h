//
//  Config.h
//  DRAGON
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#ifndef Config_h
#define Config_h

#define TESTMODE

//MARK: Algorithm Choices
#define CHOOSE_RUNTIME -1

//Riemann Solver
    #define RIEMANN_EXACT 0
        #define Exact_Rarefactions_Check //Checks for the 2-rarefaction case before attempting an iterative procedure
    #define RIEMANN_HLL 1
    #define RIEMANN_HLLC 2
    #define RIEMANN_ROE 3
        // Comment the following line to skip the Harten-Hyman Entropy Fix when using the Roe solver
        #define Harten_Hyman // If not using Roe, this line does nothing
#define RIEMANN_DEFAULT HLLC

//CFL Limit Calculation
    #define CFL_ADD 1 //Adds the speeds for dimension together
    #define CFL_MAX 0 //Uses max[ (|vx|+a)/dx, (|vy|+a)/dy, ... ]
    //Can set to any p>0 to use [ ((|vx|+a)/dx)^p + ((|vy|+a)/dy)^p +... ]^(1/p)
#define CFL_CALCULATION CFL_ADD

//Comment the following line to only use first order Godunov
#define MUSCL_Hancock
//Default to Unsplit advancement for multidimensional flows
#define DimensionUnsplit
//CTU
#define CTU








//MARK: Choose at Runtime
namespace CONFIG{
#if RIEMANN_DEFAULT == CHOOSE_RUNTIME
extern int riemann_choice;
#endif
#if CFL_CALCULATION == CHOOSE_RUNTIME || defined(TESTMODE)
//If set to p>0,  CFL uses  speed = [ ((|vx|+a)/dx)^p + ((|vy|+a)/dy)^p +... ]^(1/p)
//If set to 0, CFL uses max[ (|vx|+a)/dx, (|vy|+a)/dy, ... ]
extern int cfl_choice;
#endif
}

//MARK: Numerical Parameters

constexpr int bin_size = 125;

constexpr double CFL_coeff = 0.3;
constexpr double Timestep_Tolerance = 1e-14;

constexpr double ExactRiemann_Tolerance = 1E-12;
constexpr double ExactRiemann_MaxIters = 0; //Use a nonpostiive value for unlimited iterations

constexpr double Adaptive_PVRS_Ratio = 2.0;



#endif /* Constants_h */
