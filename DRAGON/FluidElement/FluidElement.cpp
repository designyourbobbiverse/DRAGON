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
    v = {0,0,0};
    p = 0;
#ifdef MHD
    B = {0,0,0};
#endif
}
ConservativeState::ConservativeState(){
    rho = 0;
    p = {0,0,0};
    E = 0;
#ifdef MHD
    B = {0,0,0};
#endif
}

//MARK: Type Conversion


ConservativeState::ConservativeState(PrimitiveState state){
    rho = state.rho;
    p = rho * state.v;
    E = state.energy();
#ifdef MHD
    B = state.B;
#endif
}

PrimitiveState::PrimitiveState(ConservativeState state){
    rho = state.rho;
    v = state.p / rho;
    p = state.pressure();
#ifdef MHD
    B = state.B;
#endif
}

double PrimitiveState::energy() const {
    double E =  p/(_gamma - 1.0) + (rho/2)*(v*v);
#ifdef MHD
    E += (B*B)/(8*M_PI);
#endif
    return E;
}
double ConservativeState::pressure() const {
    double E_thermal = this->E - (p*p)/(2*rho);
#ifdef MHD
     E_thermal -= (B*B)/(8*M_PI);
#endif
    return  (_gamma - 1.0) * E_thermal;
}

//Enthalpy Density
double PrimitiveState::enthalpy() const {
    double h = energy() + p;
#ifdef MHD
    h += (B*B)/(8*M_PI);
#endif
    return h / rho;
}

//MARK: Physical validity
bool PrimitiveState::isPhysical() const {
#ifdef MHD
    if( !isfinite(B.x) || !isfinite(B.y) || !isfinite(B.z)) return false;
#endif
    return isfinite(rho) && isfinite(v.x)  && isfinite(v.y)  && isfinite(v.z)  && isfinite(p)   && rho > 0.0 && p > 0.0;
}
bool ConservativeState::isFinite()  const {
#ifdef MHD
    if( !isfinite(B.x) || !isfinite(B.y) || !isfinite(B.z)) return false;
#endif
    return isfinite(rho) && isfinite(p.x)  && isfinite(p.y)  && isfinite(p.z)  && isfinite(E);
}
bool ConservativeState::isPhysical()  const {
    return isFinite() && rho > 0.0 && pressure() > 0.0;
}

//MARK: Wave Speeds
double PrimitiveState::cs(){
    return sqrt(_gamma * p / rho);
}
#ifdef MHD
double PrimitiveState::c_alfven(){
    return sqrt( (B*B) / (4*M_PI*rho) );
}
double PrimitiveState::c_fast(){
    double c = cs(), a = c_alfven(), c2 = pow(c,2), a2 = pow(a,2), ax2 = B.x*B.x/(4*M_PI*rho);
    double disc = fmax(0, pow(c2+a2,2)*0.25 - ax2*c2);
    return sqrt( (c2 + a2)*0.5 + sqrt(disc) );
}
#endif

//MARK: Flux
ConservativeState PrimitiveState::flux() const { return ConservativeState(*this).flux(v); }
ConservativeState ConservativeState::flux() const {
#ifdef MHD
    return flux(p / rho);
#else //vy and vz are never used in hydro, so don't bother calculating them
    return flux({p.x / rho, 0 , 0});
#endif
}
ConservativeState ConservativeState::flux(vec3 v) const {
    ConservativeState F = ConservativeState();
    double _p = pressure();
    F.rho = p.x;
    F.p = p * v.x;
    F.p.x += _p;
    F.E =  (E + _p) * v.x;
#ifdef MHD
    //Magnetic Pressure
    double p_mag = B*B / (8*M_PI);
    F.p.x += p_mag;
    F.E += p_mag * v.x;
    //Magnetic Tension
    vec3 mag_ten = (B.x * B) / (4*M_PI);
    F.p -= mag_ten;
    F.E -= v * mag_ten; //Poynting Flux
    //Magnetic Induction
    F.B = v.x * B - B.x * v;
#endif
    
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
#ifdef MHD
    wT.B = B.swapXY();
#endif
    return wT;
}
PrimitiveState PrimitiveState::swapXZ() const {
    PrimitiveState wT = *this;
    wT.v = v.swapXZ();
#ifdef MHD
    wT.B = B.swapXZ();
#endif
    return wT;
}
PrimitiveState PrimitiveState::swapYZ() const {
    PrimitiveState wT = *this;
    wT.v = v.swapYZ();
#ifdef MHD
    wT.B = B.swapYZ();
#endif
    return wT;
}
ConservativeState ConservativeState::swapXY() const {
    ConservativeState wT = *this;
    wT.p = p.swapXY();
#ifdef MHD
    wT.B = B.swapXY();
#endif
    return wT;
}
ConservativeState ConservativeState::swapXZ() const {
    ConservativeState wT = *this;
    wT.p = p.swapXZ();
#ifdef MHD
    wT.B = B.swapXZ();
#endif
    return wT;
}
ConservativeState ConservativeState::swapYZ() const {
    ConservativeState wT = *this;
    wT.p = p.swapYZ();
#ifdef MHD
    wT.B = B.swapYZ();
#endif
    return wT;
}

