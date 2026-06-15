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
#if RIEMANN_DEFAULT_HYDRO == CHOOSE_RUNTIME || defined(TESTMODE)
namespace CONFIG {
    int riemann_choice = RIEMANN_EXACT;
}
#endif

//Make sure to set RIEMANN_DEFAULT in Config.h
ConservativeState Riemann::flux(double dt_dx){
#if RIEMANN_DEFAULT_HYDRO == CHOOSE_RUNTIME || defined(TESTMODE)
    ConservativeState flux;
    switch (CONFIG::riemann_choice){
        case RIEMANN_HLL: flux =  HLL(); break;
        case RIEMANN_HLLC: flux =  HLLC(); break;
        case RIEMANN_ROE: flux =  Roe(); break;
        default: flux = exact().flux();
    }
#elif RIEMANN_DEFAULT_HYDRO == RIEMANN_EXACT
    auto flux = exact().flux();
#elif RIEMANN_DEFAULT_HYDRO == RIEMANN_HLL
    auto flux =  HLL();
#elif RIEMANN_DEFAULT_HYDRO == RIEMANN_HLLC
    auto flux =  HLLC();
#elif RIEMANN_DEFAULT_HYDRO == RIEMANN_ROE
    auto flux =  Roe();
#endif
//MARK: Fallback to Exact
#if RIEMANN_DEFAULT_HYDRO != RIEMANN_EXACT //Verify physicality, fallback to exact if needed
    if(dt_dx == 0) return flux;
    dt_dx *= Riemann_ExactFallback_Parameter;
    if(!(L - flux*dt_dx).isPhysical() ||  !(R+flux*dt_dx).isPhysical()){
        flux = exact().flux();
        if(!(L - flux*dt_dx).isPhysical() ||  !(R+flux*dt_dx).isPhysical()){
            //TODO: Throw Exception
        }
    }
#endif
    return flux;
}



//MARK: Dimension Convenience
ConservativeState Riemann::flux_X(double dt_dx){ return flux(dt_dx); }
ConservativeState Riemann::flux_Y(double dt_dy){  return Riemann(L.swapXY(),R.swapXY()).flux(dt_dy).swapXY(); }
ConservativeState Riemann::flux_Z(double dt_dz){  return Riemann(L.swapXZ(),R.swapXZ()).flux(dt_dz).swapXZ(); }


//MARK: Solution Sampling
//  Based mostly on Toro (2009). https://doi.org/10.1007/b79761
ConservativeState RiemannSolution::flux(){ return flux(0); }
ConservativeState RiemannSolution::flux(double x_t){
    PrimitiveState w = sample(x_t);
    return ConservativeState(w).flux(w.v.x);
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
}
