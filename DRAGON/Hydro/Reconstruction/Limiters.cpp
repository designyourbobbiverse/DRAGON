//
//  Limiters.cpp
//  DRAGON/Hydro/Reconstruction
//
//  Created by Bobbie Markwick on 31/07/2026.
//

#include "TVD.hpp"

#include "Config.h"
#include <cmath> //For std:abs
#include <algorithm> //For std::min/max
using namespace DRAGON;



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
    return (std::abs(a) < std::abs(b)) ? a : b;
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
    double c = std::abs(a+b);
    if(c < 4*std::abs(a) && c < 4*std::abs(b)) return a < 0 ? -c/2 : c/2;
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
    double _a = std::abs(a), _b = std::abs(b);
    double s1 = std::min(2*_a,_b), s2 = std::min(_a,2*_b);
    double mag = std::max(s1, s2);
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




