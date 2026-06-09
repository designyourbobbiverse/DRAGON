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
        // Comment the following line to skip the Harten_Hyman Entropy Fix when using the Roe solver
        #define Harten_Hyman // If not using Roe, this line does nothing

#define RIEMANN_DEFAULT RIEMANN_HLLC





//Comment the following line to only use first order Godunov
#define MUSCL_Hancock


//MARK: Numerical Parameters
constexpr double CFL = 0.5;
constexpr double Timestep_Tolerance = 1e-14;


constexpr double ExactSolver_Tolerance = 1E-12;
constexpr double ExactSolver_MaxIters = 0; //Use a nonpostiive value for unlimited iterations

constexpr double Adaptive_PVRS_Ratio = 2.0;



#endif /* Constants_h */
