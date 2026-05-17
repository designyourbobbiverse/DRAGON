//
//  FluidElement.cpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#include "FluidElement.hpp"
#include "Constants.h"
#include <math.h>

PrimitiveState::PrimitiveState(){
    rho = 0;
    vx = 0;
    vy = 0;
    vz = 0;
    p = 0;
}
ConservativeState::ConservativeState(){
    rho = 0;
    px = 0;
    py = 0;
    pz = 0;
    E = 0;
}

//MARK: Type Conversion


ConservativeState::ConservativeState(PrimitiveState state){
    rho = state.rho;
    px = rho * state.vx;
    py = rho * state.vy;
    pz = rho * state.vz;
    E = state.energy();
}

PrimitiveState::PrimitiveState(ConservativeState state){
    rho = state.rho;
    vx = state.px / rho;
    vy = state.py / rho;
    vz = state.pz / rho;
    p = state.pressure();
}

double PrimitiveState::energy(){
    return p/(_gamma - 1.0) + (rho/2)*(vx*vx + vy*vy + vz*vz);
}
double ConservativeState::pressure(){
    return (_gamma - 1.0) * (E - (px*px + py*py + pz*pz) / (2*rho));
}

//Enthalpy Density
double PrimitiveState::enthalpy(){
    return (_G_Gm1*p/rho) + (vx*vx + vy*vy + vz*vz)/2;
}

/*
//Speed plus sound speed
double PrimativeState::S(){
    return sqrt(vx*vx + vy*vy + vz*vz) + sqrt(gamma * p/rho);
}*/

//MARK: Physical validity
bool PrimitiveState::isPhysical() {
    return isfinite(rho) && isfinite(vx)  && isfinite(vy)  && isfinite(vz)  && isfinite(p)   && rho > 0.0 && p > 0.0;
}
bool ConservativeState::isPhysical() {
    return isfinite(rho) && isfinite(px)  && isfinite(py)  && isfinite(pz)  && isfinite(E)   && rho > 0.0 && E > 0.0;
}



//MARK: Flux
ConservativeState PrimitiveState::flux(){ return ConservativeState(*this).flux(vx); }
ConservativeState ConservativeState::flux(){ return flux(px / rho); }
ConservativeState ConservativeState::flux(double v){
    ConservativeState F = ConservativeState();
    double p = pressure();
    F.rho = px;
    F.px = px * v + p;
    F.py = py * v;
    F.pz = pz * v;
    F.E = v * (E + p);
    return F;
}

PrimitiveState& operator+=(PrimitiveState &W, const ConservativeState &dU){
    ConservativeState U(W);
    U += dU;
    W = PrimitiveState(U);
    return W;
}




//MARK: 2D + 3D
PrimitiveState PrimitiveState::swapXY(){
    PrimitiveState wT = *this;
    wT.vx = vy; wT.vy = vx;
    return wT;
}
PrimitiveState PrimitiveState::swapXZ(){
    PrimitiveState wT = *this;
    wT.vx = vz; wT.vz = vx;
    return wT;
}
PrimitiveState PrimitiveState::swapYZ(){
    PrimitiveState wT = *this;
    wT.vz = vy; wT.vy = vz;
    return wT;
}
ConservativeState ConservativeState::swapXY(){
    ConservativeState wT = *this;
    wT.px = py; wT.py = px;
    return wT;
}
ConservativeState ConservativeState::swapXZ(){
    ConservativeState wT = *this;
    wT.px = pz; wT.pz = px;
    return wT;
}
ConservativeState ConservativeState::swapYZ(){
    ConservativeState wT = *this;
    wT.pz = py; wT.py = pz;
    return wT;
}



//MARK: Arithmetic

//Arithmetic (+): Add two states together
ConservativeState operator+(ConservativeState X, const ConservativeState &Y){
    X += Y;
    return X;
}
ConservativeState& operator+=(ConservativeState &X, const ConservativeState &Y){
    X.rho += Y.rho;
    X.px += Y.px;
    X.py += Y.py;
    X.pz += Y.pz;
    X.E += Y.E;
    return X;
}

//Arithmetic (-): Subtract two states
ConservativeState operator-(ConservativeState X, const ConservativeState &Y){
    X -= Y;
    return X;
}
ConservativeState& operator-=(ConservativeState &X, const ConservativeState &Y){
    X.rho -= Y.rho;
    X.px -= Y.px;
    X.py -= Y.py;
    X.pz -= Y.pz;
    X.E -= Y.E;
    return X;
}
//Arithmetic (*): Multiply State by some scalar
ConservativeState operator*(const double &a, ConservativeState X){ return X*a; }
ConservativeState operator*(ConservativeState X, double a){
    X *= a;
    return X;
}
ConservativeState& operator*=(ConservativeState &X, double a){
    X.rho *= a;
    X.px *= a;
    X.py *= a;
    X.pz *= a;
    X.E *= a;
    return X;
}
//Arithmetic (/): Divide state by some scalar
ConservativeState operator/(ConservativeState X, double a){
    X /= a;
    return X;
}
ConservativeState& operator/=(ConservativeState &X, double a){
    X.rho /= a;
    X.px /= a;
    X.py /= a;
    X.pz /= a;
    X.E /= a;
    return X;
}



//MARK: Primitive Arithmetic (Don't)
#ifdef PRIMITIVE_ARITHMETIC_ALLOWED
PrimitiveState operator+(const PrimitiveState &X, const PrimitiveState &Y){
    PrimitiveState z = X;
    z.rho += Y.rho;
    z.vx += Y.vx;
    z.vy += Y.vy;
    z.vz += Y.vz;
    z.p += Y.p;
    return z;
}
PrimitiveState operator-(const PrimitiveState &X, const PrimitiveState &Y){
    PrimitiveState z = X;
    z.rho -= Y.rho;
    z.vx -= Y.vx;
    z.vy -= Y.vy;
    z.vz -= Y.vz;
    z.p -= Y.p;
    return z;
}
PrimitiveState operator*(const double &a, const PrimitiveState &X){ return X * a; }
PrimitiveState operator*(const PrimitiveState &X, double a){
    PrimitiveState z = X;
    z.rho *= a;
    z.vx *= a;
    z.vy *= a;
    z.vz *= a;
    z.p *= a;
    return z;
}
PrimitiveState operator/(const PrimitiveState &X, double a){
    PrimitiveState z = X;
    z.rho /= a;
    z.vx /= a;
    z.vy /= a;
    z.vz /= a;
    z.p /= a;
    return z;
}
#endif
