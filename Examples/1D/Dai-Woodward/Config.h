//
//
//  Config.h
//  DRAGON/Examples
//
//  Created by Bobbie Markwick on 7/07/2026.
//

#ifndef Config_h
#define Config_h

#include <string>

//MARK: How to use this file
//If something is indented once and sits above something else, it is one of the options for the something else
    //For example, RIEMANN default has several choices listed above it, such as RIEMANN_EXACT, RIEMANN_HLLC, RIEMANN_ROE, etc
    //For most of these choice items, you also have this option available, which allows you to choose the option at runtime (see the bottom of this header file)
        #define CHOOSE_RUNTIME -1 //Not available for choices in the File I/O section
//If something is indented and sits below something else, it is an option specific to that setting.
    //For example, RIEMANN_ROE has Harten_Hyman  below it.
        //If Harten_Hyman is enabled, the roe solver includes the entropy fix. If you turn it off, it skips that step
        //If you aren't using the Roe solver, Harten_Hyman doesn't do antying
    //As another example, RIEMANN_VERIFY_FALLBACK has the constant Riemann_ExactFallback_Parameter
        //This constant is used only in the execution of the RIEMANN_VERIFY_FALLBACK procedure
    //Indented and below generally priority over Indented and above (e.g. ExactRiemann_MaxIters is a parameter for RIEMANN_EXACT, not one of the options for RIEMANN_HLL [which doesn't even have options])


//MARK: Top level Options

#define MHD //Determines whether the simulation is run using MHD or Pure Hydrodynamics

#define DIMENSION_UNSPLIT //Use an Unsplit advancement scheme for multidimensional flows
    #define CTU //Corner Transport Upwind.  Colella (1990). https://doi.org/10.1016/0021-9991(90)90233-Q

namespace CONFIG{
//MARK: Riemann Solver
//DRAGON offers several different choices of Riemann Solver in Hydrodynamic mode, and choice of HLL/D/E in MHD
    #define RIEMANN_EXACT 0 //Produces an exact solution to the Hydrodynamic Euler Equations using an iterative procedure
        #define Exact_Rarefactions_Check //Checks for the 2-rarefaction case before attempting an iterative procedure
        constexpr double ExactRiemann_Tolerance = 1E-12; //Defines the convergence threshold for the iterative procedure
        constexpr double ExactRiemann_MaxIters = 6; //Use a nonpostiive value for unlimited iterations
    #define RIEMANN_HLL 1 // Harten, Lax, and van Leer (1983). https://doi.org/10.1137/1025002
    #define RIEMANN_HLLC 2 // Toro, Spruce, and Speares (1994). https://doi.org/10.1007/BF01414629
    #define RIEMANN_HLLD 3 // Miyoshi and Kusano (2005). https://doi.org/10.1016/j.jcp.2005.02.017
        #define HLLD_PHYSICAL_SAFETY //Falls back to HLLE if star state is nonphysical
    #define RIEMANN_HLLE 4 //Einfeldt (1988). https://doi.org/10.1137/0725021
    #define RIEMANN_HLLX 5 //HLLC for Pure Hydro, HLLD for MHD
    #define RIEMANN_ROE 6 // Roe (1981). https://doi.org/10.1016/0021-9991(81)90128-5
        #define Harten_Hyman //Entropy Fix: Harten and Hyman (1983). https://doi.org/10.1016/0021-9991(83)90066-9
#define RIEMANN_DEFAULT RIEMANN_HLLX

//When using an approximate solver, verify that L-F*dt/dx and R+F*dt/dx are physical
//Check if physical. If not, recalculate tries Exact (Hydro only). Failing that, forces a restart
#ifndef MHD //I found this option to be somewhat detrimental if HLLD_PHYSICAL_SAFETY is on
#define RIEMANN_VERIFY_FALLBACK
    constexpr double Riemann_ExactFallback_Parameter = 1.0; //Scales F*dt/dx for the purpose of physicality verificaiton
//    #define RIEMANN_FALLBACK_TRY_HLLE //try HLLE before Exact (Hydro) or restart (MHD)
#endif

//MARK: Time Control

constexpr double final_time = 0.2;
constexpr double dt = 0.1;

constexpr double CFL_coeff = 0.3; //The Coefficient used together with the above to determine the maximum timestep size
constexpr double Timestep_Tolerance = 1e-14; //Timesteps smaller than this are treated as zero

//Courant, Friedrichs, and Lewy (1928). https://doi.org/10.1007/BF01448839
//The following variants are available for computing the CFL heuristic for each cell in a multidimensional grid
    #define CFL_ADD 1 //Adds the speeds for each dimension together: (|v.x|+a)/dx + (|v.y|+a)/dy + ...
    #define CFL_MAX 0 //Uses max[ (|v.x|+a)/dx, (|v.y|+a)/dy, ...]
    //Can set to any p>0 to use [ ((|v.x|+a)/dx)^p + ((|v.y|+a)/dy)^p + ... ]^(1/p)
#define CFL_CALCULATION CFL_ADD


//MARK: MUSCL Reconstruction
//Comment the following line to only use a first order Godunov scheme
#define MUSCL_Hancock // van Leer (1979). https://doi.org/10.1016/0021-9991(79)90145-1
        #define LIMITER_MINMOD 0 //Most diffusive, prioritizes robustness
        #define LIMITER_MC 1 //Moderately compressive, sharper than minmod while still fairly well-behaved.
        #define LIMITER_VANLEER 2 //Smooth, less diffusive than minmod with gentle behavior across smooth gradients
        #define LIMITER_SUPERBEE 3 //Very compressive, sharply preserves discontinuities but can be aggressive
        #define LIMITER_VANALBADA 4 //Smooth, reduces clipping near smooth extrema while remaining shock-safe
    #define MUSCL_DEFAULT_LIMITER LIMITER_MINMOD


//MARK: File I/O
//For configuration choices in this section, CHOOSE_RUNTIME is not available
inline std::string output_base_name = "frame"; //File names will be <output_base_name>_#####.hdf
inline std::string output_dir = "/Users/bobbiemarkwick/DRAGON_OUT"; //Set this to your output directory

//#define RESTART_FROM_FILE //Attempt to restart from the output of a previous run
    #define RESTART_FRAME -1 //Use a number <0 to automatically find the latest file and restart from that

//Output parameters
#define HDF5_COMPRESSION_LEVEL 4 //HDF5 compression level. Set to a nonpositive number for no compression

    #define HDF5_WRITE_PRIMITIVE 1 //Only write the primitive values to the output files
    #define HDF5_WRITE_CONSERVATIVE 6 //Only write the conservative values to the output files
    #define HDF5_WRITE_PRIMITIVE_AND_ENERGY 3 //Write the primitive values, plus the energy density
    #define HDF5_WRITE_PRIMITIVE_AND_CONSERVATIVE 7 //Write both Primitive & Conservative values (produces larger files)
#define HDF5_WRITE_OPTION HDF5_WRITE_PRIMITIVE_AND_ENERGY //Determines which fluid representation to use for output files

//Writes ghost cells to file. Useful for debugging, but typically you want to disable this
//#define WRITE_GHOSTS_TO_FILE

    #define HDF5_WRITE_OMIT 0 //Omits B values if they can be calculated from A
    #define HDF5_WRITE_FLOAT 1 //Writes B values as Floats if they can be calculated from A
    #define HDF5_WRITE_DOUBLE 2 //Writes all B values as Doubles
#define HDF5_REDUNDANT_VALS_OPTION HDF5_WRITE_FLOAT //Output precision of certain calculable values
    //Specifically, applies to Bx (2D+3D), By (2D+3D), Bz (3D only), and E (HDF5_WRITE_PRIMITIVE_AND_ENERGY only)

//MARK: Performance

constexpr int core_count = 16; //Helps the Root level grid decide how many children to split into


//******************************************************************//
//*DO NOT TOUCH anything below this line (but feel free to read it)*//
//******************************************************************//

//MARK: Choose at Runtime
#if RIEMANN_DEFAULT == CHOOSE_RUNTIME || defined(TESTMODE)
inline int riemann_choice = RIEMANN_EXACT;
#endif
#if MUSCL_DEFAULT_LIMITER == CHOOSE_RUNTIME || defined(TESTMODE)
inline int limiter_choice = LIMITER_MINMOD;
#endif
#if CFL_CALCULATION == CHOOSE_RUNTIME || defined(TESTMODE)
//If set to p>0,  CFL uses  speed = [ ((|v.x|+a)/dx)^p + ((|v.y|+a)/dy)^p +... ]^(1/p)
//If set to 0, CFL uses max[ (|v.x|+a)/dx, (|v.y|+a)/dy, ... ]
inline int cfl_choice = CFL_ADD;
#endif
}

#endif
