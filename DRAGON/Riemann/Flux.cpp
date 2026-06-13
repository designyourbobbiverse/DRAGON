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
//Make sure to set RIEMANN_DEFAULT in Constants.h
ConservativeState Riemann::flux(){
#if RIEMANN_DEFAULT == RIEMANN_EXACT || defined(TEST_MODE)
    return exact().flux();
#elif RIEMANN_DEFAULT == RIEMANN_HLL
    return HLL();
#elif RIEMANN_DEFAULT == RIEMANN_HLLC
    return HLLC();
#elif RIEMANN_DEFAULT == RIEMANN_ROE
    return Roe();
#elif RIEMANN_DEFAULT == CHOOSE_RUNTIME
    switch (CONFIG::riemann_choice){
        case RIEMANN_HLL: return HLL();
        case RIEMANN_HLLC: return HLLC();
        case RIEMANN_ROE: return Roe();
        default: return exact().flux();
    }
#endif
}

#if RIEMANN_DEFAULT == CHOOSE_RUNTIME
namespace CONFIG {
    int riemann_choice = RIEMANN_EXACT;
}
#endif

//MARK: Dimension Convenience + Error Checking
ConservativeState Riemann::flux_X(double dt_dx){
    auto f = flux();
    #if RIEMANN_DEFAULT != RIEMANN_EXACT
    if(dt_dx > 0 && (!(L - f*dt_dx).isPhysical() ||  !(R+f*dt_dx).isPhysical())){
        f = exact().flux();
    }
    #endif
    return flux();
}
ConservativeState Riemann::flux_Y(double dt_dy){  return Riemann(L.swapXY(),R.swapXY()).flux_X(dt_dy).swapXY(); }
ConservativeState Riemann::flux_Z(double dt_dz){  return Riemann(L.swapXZ(),R.swapXZ()).flux_X(dt_dz).swapXZ(); }


//MARK: Solution Sampling
ConservativeState RiemannSolution::flux(){ return flux(0); }
ConservativeState RiemannSolution::flux(double x_t){
    PrimitiveState w = sample(x_t);
    return ConservativeState(w).flux(w.vx);
}


PrimitiveState RiemannSolution::sample(double x_t){
    PrimitiveState state;
    //Handle left vs Right side
    bool isLeft = x_t < sR.vx;
    if(isLeft){ mirror(); x_t=-x_t; }
    
    double a = sqrt(_gamma * wR.p / wR.rho);
    
    //Determine Zone
    int zone = 0; // 1 = outside, 2 = fan, 3 = star
    if (sR.p > wR.p){ //shock
        double scale = sqrt(_Gp1_2G*sR.p/wR.p + _Gm1_2G );
        zone = ( (x_t-wR.vx)  > scale * a ) ? 1 : 3;
    } else{ //Rarefraction
        if ( (x_t-wR.vx) > a ) zone = 1;
        else if ( (x_t-sR.vx) > a ) zone = 2;
        else if ( (x_t-sR.vx) > a * pow(sR.p/wR.p, _Gm1_2G) ) zone = 2;
        else zone = 3;
    }
    //Calculate the State at x/t
    switch(zone){
        case 1: state=wR; break; //Outer Region
        case 3: state=sR;break; //Star Region
        case 2://Fan
            state = wR;
            double scale = _2_Gp1 - _Gm1_Gp1 * (wR.vx-x_t)/a;
            state.rho *= pow(scale, _2_Gm1);
            state.p *= pow(scale,_2G_Gm1);
            state.vx = _2_Gp1 * (x_t-a + wR.vx * _Gm1_2);
            break;
    }
    //Be a good citizen, restore original state if we mirrored
    if(isLeft){ mirror(); state.vx *= -1;}

    return state;
}

void RiemannSolution::mirror(){
    PrimitiveState temp = wL;
    wL = wR;
    wR = temp;
    
    temp = sL;
    sL = sR;
    sR = temp;
    
    wL.vx *= -1;
    wR.vx *= -1;
    sL.vx *= -1;
    sR.vx *= -1;
}
