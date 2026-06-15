//
//  FluidElement/Arithmetic.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 29/05/2026.
//

#include "FluidElement.hpp"
#include <math.h>

//MARK: (+): Add two states together
ConservativeState operator+(ConservativeState X, const ConservativeState &Y){
    X += Y;
    return X;
}
ConservativeState& operator+=(ConservativeState &X, const ConservativeState &Y){
    X.rho += Y.rho;
    X.p += Y.p;
    X.E += Y.E;
#ifdef MHD
    X.B += Y.B;
#endif
    return X;
}
vec3 operator+(vec3 v, const vec3 &w){
    v += w;
    return v;
}
vec3& operator+=(vec3 &v, const vec3 &w){
    v.x += w.x;
    v.y += w.y;
    v.z += w.z;
    return v;
}


//MARK: (-): Subtract two states
ConservativeState operator-(ConservativeState X, const ConservativeState &Y){
    X -= Y;
    return X;
}
ConservativeState& operator-=(ConservativeState &X, const ConservativeState &Y){
    X.rho -= Y.rho;
    X.p -= Y.p;
    X.E -= Y.E;
#ifdef MHD
    X.B -= Y.B;
#endif
    return X;
}
vec3 operator-(vec3 v, const vec3 &w){
    v -= w;
    return v;
}
vec3& operator-=(vec3 &v, const vec3 &w){
    v.x -= w.x;
    v.y -= w.y;
    v.z -= w.z;
    return v;
}
//MARK: (*): Multiply State by scalar
ConservativeState operator*(const double &a, ConservativeState X){ return X*a; }
ConservativeState operator*(ConservativeState X, double a){
    X *= a;
    return X;
}
ConservativeState& operator*=(ConservativeState &X, double a){
    X.rho *= a;
    X.p *= a;
    X.E *= a;
    return X;
}
vec3 operator*(const double &a, vec3 v){ return v*a; }
vec3 operator*(vec3 v, double a){
    v *= a;
    return v;
}
vec3& operator*=(vec3 &v, double a){
    v.x *= a;
    v.y *= a;
    v.z *= a;
    return v;
}
//MARK: (*): Dot product
double operator*(vec3 v, vec3 w){
    return v.x*w.x + v.y*w.y + v.z*w.z;
}


//MARK: (/): Divide state by scalar
ConservativeState operator/(ConservativeState X, double a){
    X /= a;
    return X;
}
ConservativeState& operator/=(ConservativeState &X, double a){
    X.rho /= a;
    X.p /= a;
    X.E /= a;
    return X;
}
vec3 operator/(vec3 v, double a){
    v /= a;
    return v;
}
vec3& operator/=(vec3 &v, double a){
    v.x /= a;
    v.y /= a;
    v.z /= a;
    return v;
}


//MARK: (==) Equality
bool operator==(const PrimitiveState &X, const PrimitiveState &Y){
    if(fabs(X.rho - Y.rho) >= 1e-12) return false;
    if(X.v != Y.v) return false;
    if(fabs(X.p - Y.p) >= 1e-12) return false;
#ifdef MHD
    if(X.B != Y.B) return false;
#endif
    return true;
}
bool operator==(const ConservativeState &X, const ConservativeState &Y){
    if(fabs(X.rho - Y.rho) >= 1e-12) return false;
    if(X.p != Y.p) return false;
    if(fabs(X.E - Y.E) >= 1e-12) return false;
#ifdef MHD
    if(X.B != Y.B) return false;
#endif
    return true;
}
bool operator==(const vec3 &X, const vec3 &Y){
    if(fabs(X.x - Y.x) >= 1e-12) return false;
    if(fabs(X.y - Y.y) >= 1e-12) return false;
    if(fabs(X.z - Y.z) >= 1e-12) return false;
    return true;
}
