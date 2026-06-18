//
//  FluidElement.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#include "FluidElement.hpp"
#include "Constants.h"
#include <math.h>
#include <utility>

//MARK: Empty Constructors
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
    p = rho * state.v; //Momentum Density = Mass Denisty * Velocity
    E = state.energy();
#ifdef MHD
    B = state.B;
#endif
}

PrimitiveState::PrimitiveState(ConservativeState state){
    rho = state.rho;
    v = state.p / rho; //Bulk Velocity = Momentum Density / Mass Denisty
    p = state.pressure();
#ifdef MHD
    B = state.B;
#endif
}

double PrimitiveState::energy() const {
    double E =  p/(_gamma - 1.0) + (rho/2)*(v*v); //Thermal + Bulk Kinetic Energy Densities
#ifdef MHD
    E += (B*B)/(8*M_PI); //Magnetic Energy Density (if applicable)
#endif
    return E;
}
double ConservativeState::pressure() const {
    double E_thermal = this->E - (p*p)/(2*rho); //Remove Bulk Kinetic Energy Density
#ifdef MHD
     E_thermal -= (B*B)/(8*M_PI); //Remove Magnetic Energy Density (if applicable)
#endif
    return  (_gamma - 1.0) * E_thermal; //Convert thermal energy density to pressure
}

//Enthalpy per mass
double PrimitiveState::enthalpy() const {
    double h = energy() + p; //Energy + (Thermal) Pressure
#ifdef MHD
    h += (B*B)/(8*M_PI); //Don't forget magnetic pressure (p is thermal only)
#endif
    return h / rho; //Divide by density to get per particle
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
double PrimitiveState::cs() const {
    return sqrt(_gamma * p / rho); //Speed of sound waves in a pure hydro fluid
}
#ifdef MHD
double PrimitiveState::c_alfven() const {
    return sqrt( (B*B) / (4*M_PI*rho) ); //Alfven speed, given in gaussian units
}
double PrimitiveState::c_fast() const {
    //Precompute the sound and alvfen speeds
    double c = cs(), a = c_alfven(), c2 = pow(c,2), a2 = pow(a,2), ax2 = B.x*B.x/(4*M_PI*rho);
    // fmax(0,__) to protect against numerical issues
    double disc = fmax(0, pow(c2+a2,2)*0.25 - ax2*c2);
    //Compute fast speed
    return sqrt( (c2 + a2)*0.5 + sqrt(disc) );
}
#endif

//MARK: Flux
ConservativeState PrimitiveState::flux() const { return ConservativeState(*this).flux(v); }
ConservativeState ConservativeState::flux() const {
#ifdef MHD
    return flux(p / rho); //Compute the velocity then call the flux method
#else
    //vy and vz are never used in hydro flux, so don't waste a division on them
    return flux({p.x / rho, 0 , 0});
#endif
}
ConservativeState ConservativeState::flux(vec3 v) const {
    ConservativeState F = ConservativeState();
    double _p = pressure();
    F.rho = p.x; //Mass Flux = rho*v = Momentum density
    F.p = p * v.x; //Momentum flux from advection = \rho (v\otimes v)_x = (Momentum density) * x-velocity
    F.p.x += _p; //Thermal pressure contributes to the normal component of momentum flux
    F.E =  (E + _p) * v.x; //Energy Flux = (Enthalpy density) * normal velocity [1st law of thermodynamics]
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




//MARK: Dimension Swapping
//swap__() swaps the components in place
//swapped__() returns copies with the components swapped
void vec3::swapXY(){ std::swap(x,y); }
vec3 vec3::swappedXY() const {
    vec3 wT = *this; //Copy
    wT.swapXY(); //Swap
    return wT; //Return
}
void vec3::swapXZ(){ std::swap(x,z); }
vec3 vec3::swappedXZ() const {
    vec3 wT = *this; //Copy
    wT.swapXZ(); //Swap
    return wT; //Return
}
void vec3::swapYZ(){ std::swap(y,z); }
vec3 vec3::swappedYZ() const {
    vec3 wT = *this;
    wT.swapYZ();
    return wT;
}


void PrimitiveState::swapXY() {
    v.swapXY();
#ifdef MHD
    B.swapXY();
#endif
}
PrimitiveState PrimitiveState::swappedXY() const {
    PrimitiveState wT = *this; //Copy
    wT.swapXY();//Swap
    return wT;//Return
}

void PrimitiveState::swapXZ() {
    v.swapXZ();
#ifdef MHD
    B.swapXZ();
#endif
}
PrimitiveState PrimitiveState::swappedXZ() const {
    PrimitiveState wT = *this;//Copy
    wT.swapXZ(); //Swap
    return wT; //Return
}

void PrimitiveState::swapYZ() {
    v.swapYZ();
#ifdef MHD
    B.swapYZ();
#endif
}
PrimitiveState PrimitiveState::swappedYZ() const {
    PrimitiveState wT = *this; //Copy
    wT.swapYZ(); //Swap
    return wT; //Return
}





void ConservativeState::swapXY() {
    p.swapXY();
#ifdef MHD
    B.swapXY();
#endif
}
ConservativeState ConservativeState::swappedXY() const {
    ConservativeState uT = *this; //Copy
    uT.swapXY();//Swap
    return uT;//Return
}

void ConservativeState::swapXZ() {
    p.swapXZ();
#ifdef MHD
    B.swapXZ();
#endif
}
ConservativeState ConservativeState::swappedXZ() const {
    ConservativeState uT = *this;//Copy
    uT.swapXZ(); //Swap
    return uT; //Return
}

void ConservativeState::swapYZ() {
    p.swapYZ();
#ifdef MHD
    B.swapYZ();
#endif
}
ConservativeState ConservativeState::swappedYZ() const {
    ConservativeState uT = *this; //Copy
    uT.swapYZ(); //Swap
    return uT; //Return
}
