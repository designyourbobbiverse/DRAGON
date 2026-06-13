//
//  TVD.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 10/06/2026.
//  Implementation based mostly on Toro (2009). https://doi.org/10.1007/b79761
//      van Leer (1979). https://doi.org/10.1016/0021-9991(79)90145-1
//

#include "TVD.hpp"
#include <math.h>
#include "Config.h"


//MARK: MUSCL-Hancock

void TVD::MUSCL(const PrimitiveState& wL, PrimitiveState& _L, const PrimitiveState& wC, PrimitiveState& _R, const PrimitiveState& wR, double dt_dL){
#ifdef MUSCL_Hancock
    //Compute State
    PrimitiveState dW = TVD::limit(wC - wL, wR - wC);
    _L = wC - dW/2;
    _R = wC + dW/2;
    ConservativeState UL = ConservativeState(_L), UR = ConservativeState(_R);
    ConservativeState correction = (UR.flux() - UL.flux()) * (0.5 * dt_dL);
    _L = UL - correction;
    _R = UR - correction;
    
    if(_L.isPhysical() && _R.isPhysical() ) return; //Check to make sure this is physical, fallback to First order if not
#endif
    //First order version if MUSCL is disabled or didn't work.
    _L = wC;
    _R = wC;
}


//MARK: Limiter Selection

#if MUSCL_DEFAULT_LIMITER == CHOOSE_RUNTIME || defined(TESTMODE)
namespace CONFIG {
    int limiter_choice = LIMITER_MINMOD;
}
#endif

PrimitiveState TVD::limit(const PrimitiveState& a, const PrimitiveState& b) {
#if MUSCL_DEFAULT_LIMITER == CHOOSE_RUNTIME || defined(TESTMODE)
    switch(CONFIG::limiter_choice){
        case LIMITER_MINMOD: return TVD::minmod(a, b);
        case LIMITER_MC: return TVD::MC(a,b);
        case LIMITER_VANLEER: return TVD::vanLeer(a, b);
        case LIMITER_SUPERBEE: return superbee(a, b);
        case LIMITER_VANALBADA: return TVD::vanAlbada(a, b);
        default: return TVD::minmod(a, b);
    }
#elif MUSCL_DEFAULT_LIMITER == LIMITER_MINMOD
    return TVD::minmod(a, b);
#elif MUSCL_DEFAULT_LIMITER == LIMITER_MC
    return TVD::MC(a, b);
#elif MUSCL_DEFAULT_LIMITER == LIMITER_VANLEER
    return TVD::vanLeer(a, b);
#elif MUSCL_DEFAULT_LIMITER == LIMITER_SUPERBEE
    return TVD::superbee(a, b);
#endif
}



//MARK: MINMOD
double TVD::minmod(double a, double b) {
    if(a*b <= 0) return 0;
    if(fabs(a) < fabs(b)) return a;
    return a < 0 ? -fabs(b) : fabs(b);
}
PrimitiveState TVD::minmod(const PrimitiveState& a, const PrimitiveState& b) {
    PrimitiveState W;
    W.rho = minmod(a.rho,b.rho);
    W.vx = minmod(a.vx,b.vx);
    W.vy = minmod(a.vy,b.vy);
    W.vz = minmod(a.vz,b.vz);
    W.p = minmod(a.p,b.p);
    return W;
}

//MARK: Monotonised Central
double TVD::MC(double a, double b) {
    if(a*b <= 0) return 0;
    double c = fabs(a+b);
    if(c < 4*fabs(a) && c < 4*fabs(b)) return c/2;
    else return 2*minmod(a, b);
}
PrimitiveState TVD::MC(const PrimitiveState& a, const PrimitiveState& b) {
    PrimitiveState W;
    W.rho = MC(a.rho,b.rho);
    W.vx = MC(a.vx,b.vx);
    W.vy = MC(a.vy,b.vy);
    W.vz = MC(a.vz,b.vz);
    W.p = MC(a.p,b.p);
    return W;
}

//MARK: van Leer
double TVD::vanLeer(double a, double b) {
    if(a*b <= 0) return 0;
    return  2*a*b / (a + b);
}
PrimitiveState TVD::vanLeer(const PrimitiveState& a, const PrimitiveState& b) {
    PrimitiveState W;
    W.rho = vanLeer(a.rho,b.rho);
    W.vx = vanLeer(a.vx,b.vx);
    W.vy = vanLeer(a.vy,b.vy);
    W.vz = vanLeer(a.vz,b.vz);
    W.p = vanLeer(a.p,b.p);
    return W;
}

//MARK: SUPERBEE
double TVD::superbee(double a, double b) {
    if(a*b <= 0) return 0;
    double _a = fabs(a), _b = fabs(b);
    double s1 = fmin(2*_a,_b), s2 = fmin(_a,2*_b);
    double mag = fmax(s1, s2);
    return a < 0 ? -mag : mag;
}
PrimitiveState TVD::superbee(const PrimitiveState& a, const PrimitiveState& b) {
    PrimitiveState W;
    W.rho = superbee(a.rho,b.rho);
    W.vx = superbee(a.vx,b.vx);
    W.vy = superbee(a.vy,b.vy);
    W.vz = superbee(a.vz,b.vz);
    W.p = superbee(a.p,b.p);
    return W;
}
//MARK: van Albada
double TVD::vanAlbada(double a, double b) {
    if(a*b <= 0) return 0;
    return  a*b * (a + b) / (a*a + b*b);
}
PrimitiveState TVD::vanAlbada(const PrimitiveState& a, const PrimitiveState& b) {
    PrimitiveState W;
    W.rho = vanAlbada(a.rho,b.rho);
    W.vx = vanAlbada(a.vx,b.vx);
    W.vy = vanAlbada(a.vy,b.vy);
    W.vz = vanAlbada(a.vz,b.vz);
    W.p = vanAlbada(a.p,b.p);
    return W;
}




