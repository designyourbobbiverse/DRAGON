//
//  FluidElement.cpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#include "FluidElement.hpp"
#include "Constants.h"
#include <math.h>
#include <simd/simd.h>

PrimativeState::PrimativeState(){
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


ConservativeState::ConservativeState(PrimativeState state){
    rho = state.rho;
    px = rho * state.vx;
    py = rho * state.vy;
    pz = rho * state.vz;
    E = state.energy();
}

PrimativeState::PrimativeState(ConservativeState state){
    rho = state.rho;
    vx = state.px / rho;
    vy = state.py / rho;
    vz = state.pz / rho;
    p = state.pressure();
}

double PrimativeState::energy(){
    return p/(gamma-1.0) + (rho/2)*(vx*vx + vy*vy + vz*vz);
}
double ConservativeState::pressure(){
    return (gamma-1.0) * (E - (px*px + py*py + pz*pz) / (2*rho));
}

//Enthalpy Density
const double _G_Gm1 = gamma / (gamma - 1.0);
double PrimativeState::enthalpy(){
    return (_G_Gm1*p/rho) + (vx*vx + vy*vy + vz*vz)/2;
}

/*
//Speed plus sound speed
double PrimativeState::S(){
    return sqrt(vx*vx + vy*vy + vz*vz) + sqrt(gamma * p/rho);
}*/

//MARK: Flux
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





//MARK: 2D + 3D
PrimativeState PrimativeState::swapXY(){
    PrimativeState wT = *this;
    wT.vx = vy; wT.vy = vx;
    return wT;
}
PrimativeState PrimativeState::swapXZ(){
    PrimativeState wT = *this;
    wT.vx = vz; wT.vz = vx;
    return wT;
}
PrimativeState PrimativeState::swapYZ(){
    PrimativeState wT = *this;
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
ConservativeState operator*(ConservativeState X, const double &a){
    X *= a;
    return X;
}
ConservativeState& operator*=(ConservativeState &X, const double &a){
    X.rho *= a;
    X.px *= a;
    X.py *= a;
    X.pz *= a;
    X.E *= a;
    return X;
}
//Arithmetic (/): Divide state by some scalar
ConservativeState operator/(ConservativeState X, const double &a){
    X /= a;
    return X;
}
ConservativeState& operator/=(ConservativeState &X, const double &a){
    X.rho /= a;
    X.px /= a;
    X.py /= a;
    X.pz /= a;
    X.E /= a;
    return X;
}



//MARK: Primative Arithmetic (Don't)
#ifdef PRIMATIVE_ARITHMETIC_ALLOWED
PrimativeState operator+(const PrimativeState &X, const PrimativeState &Y){
    PrimativeState z = X;
    z.rho += Y.rho;
    z.vx += Y.vx;
    z.vy += Y.vy;
    z.vz += Y.vz;
    z.p += Y.p;
    return z;
}
PrimativeState operator-(const PrimativeState &X, const PrimativeState &Y){
    PrimativeState z = X;
    z.rho -= Y.rho;
    z.vx -= Y.vx;
    z.vy -= Y.vy;
    z.vz -= Y.vz;
    z.p -= Y.p;
    return z;
}
PrimativeState operator*(const double &a, const PrimativeState &X){ return X * a; }
PrimativeState operator*(const PrimativeState &X, const double &a){
    PrimativeState z = X;
    z.rho *= a;
    z.vx *= a;
    z.vy *= a;
    z.vz *= a;
    z.p *= a;
    return z;
}
PrimativeState operator/(const PrimativeState &X, const double &a){
    PrimativeState z = X;
    z.rho /= a;
    z.vx /= a;
    z.vy /= a;
    z.vz /= a;
    z.p /= a;
    return z;
}
#endif
