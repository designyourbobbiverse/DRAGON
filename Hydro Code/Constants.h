//
//  Constants.h
//  Hydro Code
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#ifndef Constants_h
#define Constants_h


//Physical Constants

const double gamma = 5.0/3.0;


//Numerical Parameters

const double ExactSolver_Tolerance = 1E-13;
const double ExactSolver_MaxIters = 0; //Use a nonpostiive value for unlimited iterations

const double Adaptive_PVRS_Ratio = 2.0;


//Comment to use 
#define Roe_Pike
#define Harten_Hyman




#endif /* Constants_h */
