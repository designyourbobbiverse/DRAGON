//
//  FluidElement.hpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#ifndef FluidElement_hpp
#define FluidElement_hpp

#include <stdio.h>

struct PrimitiveState;
struct ConservativeState;


struct PrimitiveState {
    double rho,  vx, vy, vz,  p;
    PrimitiveState();
    PrimitiveState(ConservativeState state);
    double energy();
    double enthalpy();
    ConservativeState flux();

    PrimitiveState swapXY();
    PrimitiveState swapXZ();
    PrimitiveState swapYZ();
    
    bool isPhysical();
};
struct ConservativeState {
    double rho,  px, py, pz, E;
    ConservativeState();
    ConservativeState(PrimitiveState prim);
    double pressure();
    
    ConservativeState flux();
    ConservativeState flux(double v);

    ConservativeState swapXY();
    ConservativeState swapXZ();
    ConservativeState swapYZ();
    
    bool isPhysical();

};



//Arithmetic (+): Add two states together
ConservativeState operator+(ConservativeState X, const ConservativeState &Y);
ConservativeState& operator+=(ConservativeState &X, const ConservativeState &Y);
//Arithmetic (-): Subtract two states
ConservativeState operator-(ConservativeState X, const ConservativeState &Y);
ConservativeState& operator-=(ConservativeState &X, const ConservativeState &Y);
//Arithmetic (*): Multiply State by some scalar
ConservativeState operator*(ConservativeState X, double a);
ConservativeState operator*(const double &a, ConservativeState X);
ConservativeState& operator*=(ConservativeState &X, double a);
//Arithmetic (/): Divide state by some scalar
ConservativeState operator/(ConservativeState X, double a);
ConservativeState& operator/=(ConservativeState &X, double a);

//Arithmetic (+): Adjust primative state by flux
PrimitiveState& operator+=(PrimitiveState &X, const ConservativeState &dU);





//Arithmetic operators for primatives
//Don't do this unless you know what you are doing.
//If you do need this, Uncomment the next line
//#define PRIMITIVE_ARITHMETIC_ALLOWED

#ifdef PRIMITIVE_ARITHMETIC_ALLOWED
PrimitiveState operator+(const PrimitiveState &X, const PrimitiveState &Y);
PrimitiveState operator-(const PrimitiveState &X, const PrimitiveState &Y);
PrimitiveState operator*(const PrimitiveState &X, double a);
PrimitiveState operator*(const double &a, const PrimitiveState &X);
PrimitiveState operator/(const PrimitiveState &X, double a);
#endif


#endif /* FluidElement_hpp */
