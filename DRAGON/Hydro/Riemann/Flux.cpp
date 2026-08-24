//
//  Flux.cpp
//  DRAGON/Hydro/Riemann
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#include "Hydro/Riemann/Riemann.hpp"

#include "Config.h"

#include <stdexcept> //For error messages
using namespace DRAGON;


//MARK: Selected Flux algorithm

//Make sure to set RIEMANN_DEFAULT in Config.h
ConservativeState Riemann::flux(double dt_dx){
    #if defined(TESTMODE) //Runnign unit tests, need all solvers to be available
    ConservativeState flux;
    switch (Config::riemann_choice) {
    case RIEMANN_EXACT: flux = exact().flux(); break;
    case RIEMANN_HLL: flux =  HLL(); break;
    case RIEMANN_HLLE: flux =  HLLE(); break;
    case RIEMANN_ROE: flux =  Roe(); break;
    #ifdef MHD
    case RIEMANN_HLLX:
    case RIEMANN_HLLD: flux = HLLD(); break;
    #else
    case RIEMANN_HLLX:
    #endif
    case RIEMANN_HLLC: flux =  HLLC(); break;
    default: flux = exact().flux();
    }
    #elif !defined(MHD) && RIEMANN_DEFAULT == CHOOSE_RUNTIME //Hydro, user wants to choose the Riemann sovler at runtime
    ConservativeState flux;
    switch (Config::riemann_choice) {
    case RIEMANN_HLLX:
    case RIEMANN_HLLC: flux =  HLLC(); break;
    case RIEMANN_EXACT: flux = exact().flux(); break;
    case RIEMANN_HLL: flux =  HLL(); break;
    case RIEMANN_HLLE: flux =  HLLE(); break;
    case RIEMANN_ROE: flux =  Roe(); break;
    default: flux = exact().flux();
    }
    #elif defined(MHD) && RIEMANN_DEFAULT == CHOOSE_RUNTIME //MHD, user wants to choose the Riemann solver at runtime
    ConservativeState flux;
    switch (Config::riemann_choice) {
    case RIEMANN_HLLX:
    case RIEMANN_HLLD: flux =  HLLD(); break;
    case RIEMANN_HLL: flux =  HLL(); break;
    case RIEMANN_HLLE: flux =  HLLE(); break;
    default: flux =  HLLD(); break;
    }
    //User chose to hardwire a particular Riemann solver
    #elif !defined(MHD) && RIEMANN_DEFAULT == RIEMANN_EXACT //Exact Solver (Hydro Only)
    auto flux = exact().flux();
    #elif RIEMANN_DEFAULT == RIEMANN_HLL //HLL Solver
    auto flux =  HLL();
    #elif !defined(MHD) && (RIEMANN_DEFAULT == RIEMANN_HLLC || RIEMANN_DEFAULT == RIEMANN_HLLX) //HLLC Solver (Hydro Only)
    auto flux =  HLLC();
    #elif defined(MHD) && (RIEMANN_DEFAULT == RIEMANN_HLLD || RIEMANN_DEFAULT == RIEMANN_HLLX) //HLLD Solver (MHD Only)
    auto flux =  HLLD();
    #elif RIEMANN_DEFAULT == RIEMANN_HLLE //HLLE Solver
    auto flux =  HLLE();
    #elif !defined(MHD) && RIEMANN_DEFAULT == RIEMANN_ROE //Roe Solver (Hydro Only)
    auto flux =  Roe();
    #endif

    #ifdef RIEMANN_VERIFY_FALLBACK
    if (dt_dx > 0) verify_and_fallback(flux, dt_dx);
    #endif
    return flux;
}
//MARK: Fallback to HLLE/Exact
//Verify that the solution produces a physical result, fallback if not
void Riemann::verify_and_fallback(ConservativeState& flux, double dt_dx){
#ifdef RIEMANN_VERIFY_FALLBACK
    dt_dx *= Config::riemann_fallback_param;//Scale time by the desired amount
    //Check whether both states would still be physical after update
    if ((L - flux*dt_dx).isPhysical() &&  (R+flux*dt_dx).isPhysical()) return;
    #ifdef RIEMANN_FALLBACK_TRY_HLLE
    //If not, try HLLE
    flux = HLLE();
    //Check whether both states would still be physical after update
    if ((L - flux*dt_dx).isPhysical() &&  (R+flux*dt_dx).isPhysical()) return;
    #endif
    #if defined(HYDRO_AVAILABLE) && RIEMANN_DEFAULT != RIEMANN_EXACT
    //If not, try Exact
    flux = exact().flux();
    //Check whether both states would still be physical after update
    if ((L - flux*dt_dx).isPhysical() &&  (R+flux*dt_dx).isPhysical()) return;
    #endif
    if (!L.isPhysical() || !R.isPhysical()) throw std::runtime_error("Riemann Solution is unphysical because inputs were not physical");
    else throw std::runtime_error("Riemann Solution to physical inputs is unphysical");
#endif
}


//MARK: Dimension Convenience
ConservativeState Riemann::flux_X(double dt_dx){ return flux(dt_dx); }
ConservativeState Riemann::flux_Y(double dt_dy){
    L.swapXY(); R.swapXY();//Swaps XY components
    auto f = flux(dt_dy); //Solve the problem as if it were X
    f.swapXY();//Swap the output back
    L.swapXY(); R.swapXY();//Swap Back Inputs to be a good citizen, even though 99% of the time we don't actually care
    return f;
}
ConservativeState Riemann::flux_Z(double dt_dz){
    L.swapXZ(); R.swapXZ();//Swaps XZ components
    auto f = flux(dt_dz); //Solve the problem as if it were X
    f.swapXZ();//Swap the output back
    L.swapXZ(); R.swapXZ();//Swap Back Inputs to be a good citizen, even though 99% of the time we don't actually care
    return f;
}


