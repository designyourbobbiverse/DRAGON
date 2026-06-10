//
//  FluidElement/Arithmetic.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 29/05/2026.
//

#include "FluidElement.hpp"


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
