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
    X.px += Y.px;
    X.py += Y.py;
    X.pz += Y.pz;
    X.E += Y.E;
    return X;
}

//MARK: (-): Subtract two states
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
//MARK: (*): Multiply State by scalar
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
//MARK: (/): Divide state by scalar
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

//MARK: (==) Equality
bool operator==(const PrimitiveState &X, const PrimitiveState &Y){
    if(fabs(X.rho - Y.rho) >= 1e-12) return false;
    if(fabs(X.vx - Y.vx) >= 1e-12) return false;
    if(fabs(X.vy - Y.vy) >= 1e-12) return false;
    if(fabs(X.vz - Y.vz) >= 1e-12) return false;
    if(fabs(X.p - Y.p) >= 1e-12) return false;

    return true;
}
bool operator==(const ConservativeState &X, const ConservativeState &Y){
    if(fabs(X.rho - Y.rho) >= 1e-12) return false;
    if(fabs(X.px - Y.px) >= 1e-12) return false;
    if(fabs(X.py - Y.py) >= 1e-12) return false;
    if(fabs(X.pz - Y.pz) >= 1e-12) return false;
    if(fabs(X.E - Y.E) >= 1e-12) return false;

    return true;
}

