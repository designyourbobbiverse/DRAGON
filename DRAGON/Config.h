//
//  Config.h
//  DRAGON
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#ifndef Config_h
#define Config_h


#define MHD


#define CHOOSE_RUNTIME -1

//MARK: Riemann Solver
//DRAGON offers several different choices of Riemann Solver in Hydrodynamic mode, and choice of HLL/D/E in MHD
    #define RIEMANN_EXACT 0 //Produces an exact solution to the Euler Equations using an iterative procedure
        #define Exact_Rarefactions_Check //Checks for the 2-rarefaction case before attempting an iterative procedure
        constexpr double ExactRiemann_Tolerance = 1E-12; //Defines the convergence threshold for the iterative procedure
        constexpr double ExactRiemann_MaxIters = 0; //Use a nonpostiive value for unlimited iterations
    #define RIEMANN_HLL 1 // Harten, Lax, and van Leer (1983). https://doi.org/10.1137/1025002
    #define RIEMANN_HLLC 2 // Toro, Spruce, and Speares (1994). https://doi.org/10.1007/BF01414629
    #define RIEMANN_HLLD 3 // Miyoshi and Kusano (2005). https://doi.org/10.1016/j.jcp.2005.02.017
    #define RIEMANN_HLLE 4 //Einfeldt (1988). https://doi.org/10.1137/0725021
    #define RIEMANN_HLLX 5 //HLLC for Pure Hydro, HLLD for MHD
    #define RIEMANN_ROE 6 // Roe (1981). https://doi.org/10.1016/0021-9991(81)90128-5
        #define Harten_Hyman //Harten and Hyman (1983). https://doi.org/10.1016/0021-9991(83)90066-9
#define RIEMANN_DEFAULT RIEMANN_HLLX

//When using an approximate solver, verify that L-F*dt/dx and R+F*dt/dx are physical + recalculate using HLLE or Exact if not
#define RIEMANN_VERIFY_FALLBACK
    constexpr double Riemann_ExactFallback_Parameter = 1.0; //Scales F*dt/dx for the purpose of physicality verificaiton


//MARK: CFL Calculation
//Courant, Friedrichs, and Lewy (1928). https://doi.org/10.1007/BF01448839
//The following variants are available for computing the CFL time limit for each cell in a multidimensional grid
    #define CFL_ADD 1 //Adds the speeds for each dimension together: (|v.x|+a)/dx + (|v.y|+a)/dy + ...
    #define CFL_MAX 0 //Uses max[ (|v.x|+a)/dx, (|v.y|+a)/dy, ...]
    //Can set to any p>0 to use [ ((|v.x|+a)/dx)^p + ((|v.y|+a)/dy)^p + ... ]^(1/p)
#define CFL_CALCULATION CFL_ADD

constexpr double CFL_coeff = 0.3; //The Coefficient used together with the above to determine the maximum timestep size
constexpr double Timestep_Tolerance = 1e-14; //Timesteps smaller than this are treated as zero

//MARK: MUSCL Reconstruction
//Comment the following line to only use a first order Godunov scheme
#define MUSCL_Hancock // van Leer (1979). https://doi.org/10.1016/0021-9991(79)90145-1
        #define LIMITER_MINMOD 0
        #define LIMITER_MC 1
        #define LIMITER_VANLEER 2
        #define LIMITER_SUPERBEE 3
        #define LIMITER_VANALBADA 4
    #define MUSCL_DEFAULT_LIMITER LIMITER_MINMOD

//MARK: Multi-dimension
#define DimensionUnsplit //Use an Unsplit advancement scheme for multidimensional flows
    #define CTU //Corner Transport Upwind.  Colella (1990). https://doi.org/10.1016/0021-9991(90)90233-Q


//MARK: Grid Operation

constexpr int bin_size = 250; //AMRGrid larger than this will split into child grids of this size or smaller





//MARK: Choose at Runtime
namespace CONFIG{
#if RIEMANN_DEFAULT == CHOOSE_RUNTIME || defined(TESTMODE)
extern int riemann_choice;
#endif
#if MUSCL_DEFAULT_LIMITER == CHOOSE_RUNTIME || defined(TESTMODE)
extern int limiter_choice;
#endif
#if CFL_CALCULATION == CHOOSE_RUNTIME || defined(TESTMODE)
//If set to p>0,  CFL uses  speed = [ ((|v.x|+a)/dx)^p + ((|v.y|+a)/dy)^p +... ]^(1/p)
//If set to 0, CFL uses max[ (|v.x|+a)/dx, (|v.y|+a)/dy, ... ]
extern int cfl_choice;
#endif
}

//Macro which tells hydro only code to not compile if using MHD, except for unit testing where both modes are needed
#if !defined(MHD) || defined(TESTMODE)
#define HYDRO_AVAILABLE
#endif

#endif /* Constants_h */
