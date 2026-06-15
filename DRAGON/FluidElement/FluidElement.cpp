//
//  FluidElement.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#include "FluidElement.hpp"
#include "Constants.h"
#include <math.h>

PrimitiveState::PrimitiveState(){
    rho = 0;
    v.x = 0;
    v.y = 0;
    v.z = 0;
    p = 0;
}
ConservativeState::ConservativeState(){
    rho = 0;
    p.x = 0;
    p.y = 0;
    p.z = 0;
    E = 0;
}

//MARK: Type Conversion


ConservativeState::ConservativeState(PrimitiveState state){
    rho = state.rho;
    p.x = rho * state.v.x;
    p.y = rho * state.v.y;
    p.z = rho * state.v.z;
    E = state.energy();
}

PrimitiveState::PrimitiveState(ConservativeState state){
    rho = state.rho;
    v.x = state.p.x / rho;
    v.y = state.p.y / rho;
    v.z = state.p.z / rho;
    p = state.pressure();
}

double PrimitiveState::energy() const {
    return p/(_gamma - 1.0) + (rho/2)*(v.x*v.x + v.y*v.y + v.z*v.z);
}
double ConservativeState::pressure() const {
    return (_gamma - 1.0) * (E - (p.x*p.x + p.y*p.y + p.z*p.z) / (2*rho));
}

//enthalpy Density
double PrimitiveState::enthalpy() const {
    return (_G_Gm1*p/rho) + (v.x*v.x + v.y*v.y + v.z*v.z)/2;
}

/*
//Speed plus sound speed
double PrimativeState::S(){
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z) + sqrt(gamma * p/rho);
}*/

//MARK: Physical validity
bool PrimitiveState::isPhysical() const {
    return isfinite(rho) && isfinite(v.x)  && isfinite(v.y)  && isfinite(v.z)  && isfinite(p)   && rho > 0.0 && p > 0.0;
}
bool ConservativeState::isPhysical()  const {
    return isfinite(rho) && isfinite(p.x)  && isfinite(p.y)  && isfinite(p.z)  && isfinite(E)   && rho > 0.0 && pressure() > 0.0;
}
bool ConservativeState::isFinite()  const {
    return isfinite(rho) && isfinite(p.x)  && isfinite(p.y)  && isfinite(p.z)  && isfinite(E);
}


//MARK: Flux
ConservativeState PrimitiveState::flux() const { return ConservativeState(*this).flux(v.x); }
ConservativeState ConservativeState::flux() const { return flux(p.x / rho); }
ConservativeState ConservativeState::flux(double v) const {
    ConservativeState F = ConservativeState();
    double _p = pressure();
    F.rho = p.x;
    F.p.x = p.x * v + _p;
    F.p.y = p.y * v;
    F.p.z = p.z * v;
    F.E = v * (E + _p);
    return F;
}

PrimitiveState& operator+=(PrimitiveState &W, const ConservativeState &dU){
    ConservativeState U(W);
    U += dU;
    W = PrimitiveState(U);
    return W;
}




//MARK: 2D + 3D
vec3 vec3::swapXY() const {
    vec3 wT = *this;
    wT.x = y; wT.y = x;
    return wT;
}
vec3 vec3::swapXZ() const {
    vec3 wT = *this;
    wT.x = z; wT.z = x;
    return wT;
}
vec3 vec3::swapYZ() const {
    vec3 wT = *this;
    wT.z = y; wT.y = z;
    return wT;
}
PrimitiveState PrimitiveState::swapXY() const {
    PrimitiveState wT = *this;
    wT.v = v.swapXY();
    return wT;
}
PrimitiveState PrimitiveState::swapXZ() const {
    PrimitiveState wT = *this;
    wT.v = v.swapXZ();
    return wT;
}
PrimitiveState PrimitiveState::swapYZ() const {
    PrimitiveState wT = *this;
    wT.v = v.swapYZ();
    return wT;
}
ConservativeState ConservativeState::swapXY() const {
    ConservativeState wT = *this;
    wT.p = p.swapXY();
    return wT;
}
ConservativeState ConservativeState::swapXZ() const {
    ConservativeState wT = *this;
    wT.p = p.swapXZ();
    return wT;
}
ConservativeState ConservativeState::swapYZ() const {
    ConservativeState wT = *this;
    wT.p = p.swapYZ();
    return wT;
}

