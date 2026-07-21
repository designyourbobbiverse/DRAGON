//
//  Arithmetic.cpp
//  DRAGON/FluidElement
//
//  Created by Bobbie Markwick on 29/05/2026.
//

#include "FluidElement.hpp"
#include <math.h>



//MARK: (+): Add two states together
ConservativeState operator+(ConservativeState X, const ConservativeState &Y){
    //Compiler already made a copy of X, so we can just add to that and return it
    //This approach ensures consistency with +=
    X += Y;
    return X;
}
ConservativeState& operator+=(ConservativeState &X, const ConservativeState &Y){
    X.rho += Y.rho;
    X.mom += Y.mom;
    X.E += Y.E;
#ifdef MHD
    X.B += Y.B;
#endif
    return X;
}
vec3 operator+(vec3 v, const vec3 &w){
    v += w;//Compiler already made a copy of V, so we can just add to that and return it
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
    X -= Y; //Compiler already made a copy of X, so we can just subtract from that and return it
    return X;
}
ConservativeState& operator-=(ConservativeState &X, const ConservativeState &Y){
    X.rho -= Y.rho;
    X.mom -= Y.mom;
    X.E -= Y.E;
#ifdef MHD
    X.B -= Y.B;
#endif
    return X;
}
vec3 operator-(vec3 v, const vec3 &w){
    v -= w;//Compiler already made a copy of v, so we can just subtract from that and return it
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
    X *= a; //Compiler already made a copy of X, so we can just multiply that and return it
    return X;
}
ConservativeState& operator*=(ConservativeState &X, double a){
    X.rho *= a;
    X.mom *= a;
    X.E *= a;
#ifdef MHD
    X.B *= a;
#endif
    return X;
}
vec3 operator*(const double &a, vec3 v){ return v*a; }
vec3 operator*(vec3 v, double a){
    v *= a; //Compiler already made a copy of v, so we can just multiply that and return it
    return v;
}
vec3& operator*=(vec3 &v, double a){
    v.x *= a;
    v.y *= a;
    v.z *= a;
    return v;
}
//MARK: Dot + Cross product
double operator*(const vec3& v, const vec3& w){
    return v.x*w.x + v.y*w.y + v.z*w.z;
}
vec3 cross(const vec3& u, const vec3& v){
    vec3 w;
    w.x = u.y * v.z - u.z * v.y;
    w.y = u.z * v.x - u.x * v.z;
    w.z = u.x * v.y - u.y * v.x;
    return w;
}



//MARK: (/): Divide state by scalar
ConservativeState operator/(ConservativeState X, double a){
    X /= a; //Compiler already made a copy of X, so we can just divide that and return it
    return X;
}
ConservativeState& operator/=(ConservativeState &X, double a){
    X.rho /= a;
    X.mom /= a;
    X.E /= a;
#ifdef MHD
    X.B /= a;
#endif
    return X;
}
vec3 operator/(vec3 v, double a){
    v /= a; //Compiler already made a copy of v, so we can just divide that and return it
    return v;
}
vec3& operator/=(vec3 &v, double a){
    v.x /= a;
    v.y /= a;
    v.z /= a;
    return v;
}


//MARK: (==) Equality
//Each of these checks to within a 1e-12 tolerance to allow for floating point errors
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
    if(X.mom != Y.mom) return false;
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
