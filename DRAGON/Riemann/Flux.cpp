//
//  Riemann/Flux.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#include "Riemann.hpp"
#include <math.h>
#include "Constants.h"
#include "Config.h"


//MARK: Selected Flux algorithm
#if RIEMANN_DEFAULT == CHOOSE_RUNTIME || defined(TESTMODE)
namespace CONFIG {
    int riemann_choice = RIEMANN_EXACT;
}
#endif

//Make sure to set RIEMANN_DEFAULT in Config.h
ConservativeState Riemann::flux(double dt_dx){
#if defined(TESTMODE)
    ConservativeState flux;
    switch (CONFIG::riemann_choice){
        case RIEMANN_EXACT: flux = exact().flux();
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
#elif !defined(MHD) && RIEMANN_DEFAULT == CHOOSE_RUNTIME
    ConservativeState flux;
    switch (CONFIG::riemann_choice){
        case RIEMANN_HLLX:
        case RIEMANN_HLLC: flux =  HLLC(); break;
        case RIEMANN_EXACT: flux = exact().flux();
        case RIEMANN_HLL: flux =  HLL(); break;
        case RIEMANN_HLLE: flux =  HLLE(); break;
        case RIEMANN_ROE: flux =  Roe(); break;
        default: flux = exact().flux();
    }
#elif defined(MHD) && RIEMANN_DEFAULT == CHOOSE_RUNTIME
    ConservativeState flux;
    switch (CONFIG::riemann_choice){
        case RIEMANN_HLLX:
        case RIEMANN_HLLD: flux =  HLLD(); break;
        case RIEMANN_HLL: flux =  HLL(); break;
        case RIEMANN_HLLE: flux =  HLLE(); break;
        default: flux =  HLLD(); break;
    }
#elif !defined(MHD) && RIEMANN_DEFAULT == RIEMANN_EXACT
    auto flux = exact().flux();
#elif RIEMANN_DEFAULT == RIEMANN_HLL
    auto flux =  HLL();
#elif !defined(MHD) && (RIEMANN_DEFAULT == RIEMANN_HLLC || RIEMANN_DEFAULT == RIEMANN_HLLX)
    auto flux =  HLLC();
#elif defined(MHD) && (RIEMANN_DEFAULT == RIEMANN_HLLD || RIEMANN_DEFAULT == RIEMANN_HLLX)
    auto flux =  HLLD();
#elif RIEMANN_DEFAULT == RIEMANN_HLLE
    auto flux =  HLLE();
#elif !defined(MHD) && RIEMANN_DEFAULT == RIEMANN_ROE
    auto flux =  Roe();
#endif
#ifdef RIEMANN_VERIFY_FALLBACK
    if(dt_dx > 0) verify_and_fallback(flux, dt_dx);
#endif
    return flux;
}
//MARK: Fallback to Exact
void Riemann::verify_and_fallback(ConservativeState& flux, double dt_dx){
    dt_dx *= Riemann_ExactFallback_Parameter;
    if((L - flux*dt_dx).isPhysical() &&  (R+flux*dt_dx).isPhysical()) return;
    flux = HLLE();
    if((L - flux*dt_dx).isPhysical() &&  (R+flux*dt_dx).isPhysical()) return;
    #if defined(TESTMODE) || (!defined(MHD) && RIEMANN_DEFAULT != RIEMANN_EXACT)
    flux = exact().flux();
    if((L - flux*dt_dx).isPhysical() &&  (R+flux*dt_dx).isPhysical()) return;
    #endif
    //TODO: Throw Exception
}


//MARK: Dimension Convenience
ConservativeState Riemann::flux_X(double dt_dx){ return flux(dt_dx); }
ConservativeState Riemann::flux_Y(double dt_dy){
    L.swapXY(); R.swapXY();
    auto f = flux(dt_dy);
    L.swapXY(); R.swapXY();//Swap Back to be a good citizen, even though 99% of the time we don't actually care
    f.swapXY();
    return f;
}
ConservativeState Riemann::flux_Z(double dt_dz){
    L.swapXZ(); R.swapXZ();
    auto f = flux(dt_dz);
    L.swapXZ(); R.swapXZ();//Swap Back to be a good citizen, even though 99% of the time we don't actually care
    f.swapXZ();
    return f;
}


//MARK: Solution Sampling
//  Based mostly on Toro (2009). https://doi.org/10.1007/b79761
ConservativeState RiemannSolution::flux(){ return flux(0); }
ConservativeState RiemannSolution::flux(double x_t){
    PrimitiveState w = sample(x_t);
    return ConservativeState(w).flux(w.v);
}


PrimitiveState RiemannSolution::sample(double x_t){
    PrimitiveState state;
    //Handle left vs Right side
    bool isLeft = x_t < sR.v.x;
    if(isLeft){ mirror(); x_t=-x_t; }
    
    double a = sqrt(_gamma * wR.p / wR.rho);
    
    //Determine Zone
    int zone = 0; // 1 = outside, 2 = fan, 3 = star
    if (sR.p > wR.p){ //shock
        double scale = sqrt(_Gp1_2G*sR.p/wR.p + _Gm1_2G );
        zone = ( (x_t-wR.v.x)  > scale * a ) ? 1 : 3;
    } else{ //Rarefraction
        if ( (x_t-wR.v.x) > a ) zone = 1;
        else if ( (x_t-sR.v.x) > a ) zone = 2;
        else if ( (x_t-sR.v.x) > a * pow(sR.p/wR.p, _Gm1_2G) ) zone = 2;
        else zone = 3;
    }
    //Calculate the State at x/t
    switch(zone){
        case 1: state=wR; break; //Outer Region
        case 3: state=sR;break; //Star Region
        case 2://Fan
            state = wR;
            double scale = _2_Gp1 - _Gm1_Gp1 * (wR.v.x-x_t)/a;
            state.rho *= pow(scale, _2_Gm1);
            state.p *= pow(scale,_2G_Gm1);
            state.v.x = _2_Gp1 * (x_t-a + wR.v.x * _Gm1_2);
            break;
    }
    //Be a good citizen, restore original state if we mirrored
    if(isLeft){ mirror(); state.v.x *= -1;}

    return state;
}

void RiemannSolution::mirror(){
    PrimitiveState temp = wL;
    wL = wR;
    wR = temp;
    
    temp = sL;
    sL = sR;
    sR = temp;
    
    wL.v.x *= -1;
    wR.v.x *= -1;
    sL.v.x *= -1;
    sR.v.x *= -1;
#ifdef MHD
    wL.B.x *= -1;
    wR.B.x *= -1;
    sL.B.x *= -1;
    sR.B.x *= -1;
#endif
    
}
