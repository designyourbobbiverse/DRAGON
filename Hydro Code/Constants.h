//
//  Constants.h
//  Hydro Code
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#ifndef Constants_h
#define Constants_h


//MARK: Physical Constants

constexpr double _gamma = 5.0/3.0;

//MARK: Grid Parameters
const double stepSize = 0.01;
const double timeStep = 0.01;

//MARK: Algorithm Choices
#define CHOOSE_RUNTIME -1

//Riemann Solver algorithm
    #define RIEMANN_EXACT 0
    #define RIEMANN_HLL 1
    #define RIEMANN_HLLC 2
    #define RIEMANN_ROE 3
#define RIEMANN_DEFAULT RIEMANN_HLLC

#if RIEMANN_DEFAULT == CHOOSE_RUNTIME
int dynamicSolverChoice = RIEMANN_EXACT;
#endif

// Comment the following line to skip the Harten_Hyman Entropy Fix when using the Roe solver
// If not using Roe, this line does nothing
#define Harten_Hyman




//Comment the following line to only use first order Godunov
#define MUSCL_Hancock


//MARK: Numerical Parameters

const double ExactSolver_Tolerance = 1E-13;
const double ExactSolver_MaxIters = 0; //Use a nonpostiive value for unlimited iterations

const double Adaptive_PVRS_Ratio = 2.0;











//MARK: Gamma-dependent constants
constexpr double  _Gm1_2G = (_gamma - 1.0)/(2.0 * _gamma);
constexpr double _Gp1_2G = (_gamma + 1.0)/(2.0*_gamma);
constexpr double _2G_Gm1 = 2.0*_gamma/(_gamma - 1.0);
constexpr double _G_Gm1 = _gamma / (_gamma - 1.0);
constexpr double _Gm1_2 = (_gamma - 1.0)/2.0;
constexpr double _2_Gm1 = 2.0/(_gamma - 1.0);
constexpr double _2_Gp1 = 2.0/(_gamma + 1.0);
constexpr double _Gm1_Gp1 = (_gamma - 1.0)/(_gamma + 1.0);
constexpr double _Ginv = 1.0/_gamma;

#endif /* Constants_h */
