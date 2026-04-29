//
//  FluidElement.hpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#ifndef FluidElement_hpp
#define FluidElement_hpp

#include <stdio.h>

struct PrimativeState;
struct ConservativeState;


struct PrimativeState {
    double rho,  vx, vy, vz,  p;
    PrimativeState();
    PrimativeState(ConservativeState state);
    double energy();
    double enthalpy();
    //double S();
    PrimativeState swapXY();
    PrimativeState swapXZ();
    PrimativeState swapYZ();
};
struct ConservativeState {
    double rho,  px, py, pz, E;
    ConservativeState();
    ConservativeState(PrimativeState prim);
    double pressure();
    ConservativeState flux(double v);
    ConservativeState swapXY();
    ConservativeState swapXZ();
    ConservativeState swapYZ();
};



//Arithmetic (+): Add two states together
ConservativeState operator+(ConservativeState X, const ConservativeState &Y);
ConservativeState& operator+=(ConservativeState &X, const ConservativeState &Y);
//Arithmetic (-): Subtract two states
ConservativeState operator-(ConservativeState X, const ConservativeState &Y);
ConservativeState& operator-=(ConservativeState &X, const ConservativeState &Y);
//Arithmetic (*): Multiply State by some scalar
ConservativeState operator*(ConservativeState X, const double &a);
ConservativeState operator*(const double &a, ConservativeState X);
ConservativeState& operator*=(ConservativeState &X, const double &a);
//Arithmetic (/): Divide state by some scalar
ConservativeState operator/(ConservativeState X, const double &a);
ConservativeState& operator/=(ConservativeState &X, const double &a);






//Arithmetic operators for primatives
//Don't do this unless you know what you are doing.
//If you do need this, Uncomment the next line
//#define PRIMATIVE_ARITHMETIC_ALLOWED

#ifdef PRIMATIVE_ARITHMETIC_ALLOWED
PrimativeState operator+(const PrimativeState &X, const PrimativeState &Y);
PrimativeState operator-(const PrimativeState &X, const PrimativeState &Y);
PrimativeState operator*(const PrimativeState &X, const double &a);
PrimativeState operator*(const double &a, const PrimativeState &X);
PrimativeState operator/(const PrimativeState &X, const double &a);
#endif


#endif /* FluidElement_hpp */
