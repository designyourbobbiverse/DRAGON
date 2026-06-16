//
//  FluidElement.hpp
//  DRAGON
//  User-Facing Header file
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#ifndef FluidElement_hpp
#define FluidElement_hpp

#include "Config.h"

struct PrimitiveState;
struct ConservativeState;

struct vec3{
    double x, y, z;
    vec3 swapXY() const;
    vec3 swapXZ() const;
    vec3 swapYZ() const;
};


struct PrimitiveState {
    double rho, p;
    vec3 v;
#ifdef MHD
    vec3 B;
#endif
    
    PrimitiveState();
    PrimitiveState(ConservativeState state);
    double energy() const;
    double enthalpy() const;
    
    double cs();
#ifdef MHD
    double c_alfven();
    double c_fast();
#endif
    
    ConservativeState flux() const;

    PrimitiveState swapXY() const;
    PrimitiveState swapXZ() const;
    PrimitiveState swapYZ() const;
    
    bool isPhysical() const;
};
struct ConservativeState {
    double rho, E;
    vec3 p;
#ifdef MHD
    vec3 B;
#endif
    
    ConservativeState();
    ConservativeState(PrimitiveState prim);
    double pressure() const;
    ConservativeState flux() const;
    ConservativeState flux(vec3 v) const; //More effiecient version if you already know vx

    ConservativeState swapXY() const;
    ConservativeState swapXZ() const;
    ConservativeState swapYZ() const;
    
    bool isPhysical() const;
    bool isFinite() const;
};



bool operator==(const PrimitiveState &X, const PrimitiveState &Y);
bool operator==(const ConservativeState &X, const ConservativeState &Y);
bool operator==(const vec3 &v, const vec3 &w);


//Arithmetic (+): Add two states together
ConservativeState operator+(ConservativeState X, const ConservativeState &Y);
ConservativeState& operator+=(ConservativeState &X, const ConservativeState &Y);
vec3 operator+(vec3 v, const vec3 &w);
vec3& operator+=(vec3 &v, const vec3 &w);
//Arithmetic (-): Subtract two states
ConservativeState operator-(ConservativeState X, const ConservativeState &Y);
ConservativeState& operator-=(ConservativeState &X, const ConservativeState &Y);
vec3 operator-(vec3 v, const vec3 &w);
vec3& operator-=(vec3 &v, const vec3 &w);
//Arithmetic (*): Multiply State by some scalar
ConservativeState operator*(ConservativeState X, double a);
ConservativeState operator*(const double &a, ConservativeState X);
ConservativeState& operator*=(ConservativeState &X, double a);
vec3 operator*(vec3 v, double a);
vec3 operator*(const double &a, vec3 v);
vec3& operator*=(vec3 &v, double a);
//Arithmetic (*): Dot product
double operator*(vec3 v, vec3 w);

//Arithmetic (/): Divide state by some scalar
ConservativeState operator/(ConservativeState X, double a);
ConservativeState& operator/=(ConservativeState &X, double a);
vec3 operator/(vec3 v, double a);
vec3& operator/=(vec3 &v, double a);

//Arithmetic (+): Adjust primative state by flux
PrimitiveState& operator+=(PrimitiveState &X, const ConservativeState &dU);






#endif /* FluidElement_hpp */
