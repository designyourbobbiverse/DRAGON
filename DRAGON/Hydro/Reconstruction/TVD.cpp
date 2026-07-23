//
//  TVD.cpp
//  DRAGON/Hydro/TVD
//
//  Created by Bobbie Markwick on 10/06/2026.
//  Implementation based mostly on Toro (2009). https://doi.org/10.1007/b79761
//      van Leer (1979). https://doi.org/10.1016/0021-9991(79)90145-1
//

#include "TVD.hpp"
#include <math.h>
#include "Config.h"
#include "Constants.h"


//MARK: MUSCL-Hancock
void TVD::MUSCL(const PrimitiveState& wL, PrimitiveState& _L, const PrimitiveState& wC, PrimitiveState& _R, const PrimitiveState& wR, double dt_dL, const vec3& dB){
#ifdef MUSCL_Hancock
    // Reconstruct a limited primitive-variable slope
    PrimitiveState leftDiff, rightDiff;
    leftDiff.rho = wC.rho - wL.rho; rightDiff.rho = wR.rho - wC.rho;
    leftDiff.v = wC.v - wL.v; rightDiff.v = wR.v - wC.v;
    leftDiff.p = wC.p - wL.p; rightDiff.p = wR.p - wC.p;
#ifdef MHD
    leftDiff.B = wC.B - wL.B; rightDiff.B = wR.B - wC.B;
#endif
    PrimitiveState dW = TVD::limit(leftDiff, rightDiff);
    //Spatial half step using the limited slope
    _L.rho = wC.rho - 0.5*dW.rho;
    _R.rho = wC.rho + 0.5*dW.rho;
    _L.v = wC.v - 0.5*dW.v;
    _R.v = wC.v + 0.5*dW.v;
    _L.p = wC.p - 0.5*dW.p;
    _R.p = wC.p + 0.5*dW.p;
#ifdef MHD
    _L.B = wC.B - 0.5*dW.B;
    _R.B = wC.B + 0.5*dW.B;
#endif
    
    //Time half step (MUSCL-Hancock Predictor)
    ConservativeState UL = ConservativeState(_L), UR = ConservativeState(_R);
    ConservativeState correction = (UR.flux(_R.v) - UL.flux(_L.v)) * (0.5 * dt_dL);
    #ifdef MHD//MHD Source Terms
    const double a_s = (_R.B.x - _L.B.x) * dt_dL;
    const double vy = 0.5*(_L.v.y + _R.v.y), vz = 0.5*(_L.v.z + _R.v.z);
    const double dBy_ = 0.5 * wC.v.y * TVD::minmod(dB.x, -dB.y) - 0.5*vy*a_s;
    const double dBz_ = 0.5 * wC.v.z * TVD::minmod(dB.x, -dB.z) - 0.5*vz*a_s;

    correction.B.y -= dBy_;
    correction.B.z -= dBz_;
    correction.E   -= _1_4pi * (wC.B.y*dBy_ + wC.B.z*dBz_);   // hold thermal p fixed
    #endif
    // Correction is applied in conservative, then automatically converted back to primitive
    _L = UL - correction;
    _R = UR - correction;
    
    if(_L.isPhysical() && _R.isPhysical() ) return; //Check to make sure this is physical, fallback to First order if not
#endif
    //First order version if MUSCL is disabled or didn't work.
    _L = wC;
    _R = wC;
}


//MARK: Limiter Selection

PrimitiveState TVD::limit(const PrimitiveState& a, const PrimitiveState& b) {
#if MUSCL_DEFAULT_LIMITER == CHOOSE_RUNTIME || defined(TESTMODE)
    //User wants to choose their limiter at runtime (or is unit testing and needs all of them)
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
#elif MUSCL_DEFAULT_LIMITER == LIMITER_VANALBADA
    return TVD::vanAlbada(a, b);
#endif
}



//MARK: MINMOD
double TVD::minmod(double a, double b) {
    if(a*b <= 0) return 0;
    return (fabs(a) < fabs(b)) ? a : b;
}
//Apply Limiter to each component of a vector
vec3 TVD::minmod(const vec3& a, const vec3& b) {
    vec3 v;
    v.x = minmod(a.x,b.x);
    v.y = minmod(a.y,b.y);
    v.z = minmod(a.z,b.z);
    return v;
}
PrimitiveState TVD::minmod(const PrimitiveState& a, const PrimitiveState& b) {
    PrimitiveState W;
    W.rho = minmod(a.rho,b.rho);
    W.v = minmod(a.v,b.v);
    W.p = minmod(a.p,b.p);
#ifdef MHD
    W.B = minmod(a.B,b.B);
#endif
    return W;
}

//MARK: Monotone-Centred
double TVD::MC(double a, double b) {
    if(a*b <= 0) return 0;
    double c = fabs(a+b);
    if(c < 4*fabs(a) && c < 4*fabs(b)) return a < 0 ? -c/2 : c/2;
    else return 2*minmod(a, b);
}
//Apply Limiter to each component of a vector
vec3 TVD::MC(const vec3& a, const vec3& b) {
    vec3 v;
    v.x = MC(a.x,b.x);
    v.y = MC(a.y,b.y);
    v.z = MC(a.z,b.z);
    return v;
}
PrimitiveState TVD::MC(const PrimitiveState& a, const PrimitiveState& b) {
    PrimitiveState W;
    W.rho = MC(a.rho,b.rho);
    W.v = MC(a.v,b.v);
    W.p = MC(a.p,b.p);
#ifdef MHD
    W.B = MC(a.B,b.B);
#endif
    return W;
}

//MARK: van Leer
double TVD::vanLeer(double a, double b) {
    if(a*b <= 0) return 0;
    return  2*a*b / (a + b);
}
//Apply Limiter to each component of a vector
vec3 TVD::vanLeer(const vec3& a, const vec3& b) {
    vec3 v;
    v.x = vanLeer(a.x,b.x);
    v.y = vanLeer(a.y,b.y);
    v.z = vanLeer(a.z,b.z);
    return v;
}
PrimitiveState TVD::vanLeer(const PrimitiveState& a, const PrimitiveState& b) {
    PrimitiveState W;
    W.rho = vanLeer(a.rho,b.rho);
    W.v = vanLeer(a.v,b.v);
    W.p = vanLeer(a.p,b.p);
#ifdef MHD
    W.B = vanLeer(a.B,b.B);
#endif
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
//Apply Limiter to each component of a vector
vec3 TVD::superbee(const vec3& a, const vec3& b) {
    vec3 v;
    v.x = superbee(a.x,b.x);
    v.y = superbee(a.y,b.y);
    v.z = superbee(a.z,b.z);
    return v;
}
PrimitiveState TVD::superbee(const PrimitiveState& a, const PrimitiveState& b) {
    PrimitiveState W;
    W.rho = superbee(a.rho,b.rho);
    W.v = superbee(a.v,b.v);
    W.p = superbee(a.p,b.p);
#ifdef MHD
    W.B = superbee(a.B,b.B);
#endif
    return W;
}
//MARK: van Albada
double TVD::vanAlbada(double a, double b) {
    if(a*b <= 0) return 0;
    return  a*b * (a + b) / (a*a + b*b);
}
//Apply Limiter to each component of a vector
vec3 TVD::vanAlbada(const vec3& a, const vec3& b) {
    vec3 v;
    v.x = vanAlbada(a.x,b.x);
    v.y = vanAlbada(a.y,b.y);
    v.z = vanAlbada(a.z,b.z);
    return v;
}
PrimitiveState TVD::vanAlbada(const PrimitiveState& a, const PrimitiveState& b) {
    PrimitiveState W;
    W.rho = vanAlbada(a.rho,b.rho);
    W.v = vanAlbada(a.v,b.v);
    W.p = vanAlbada(a.p,b.p);
#ifdef MHD
    W.B = vanAlbada(a.B,b.B);
#endif
    return W;
}




