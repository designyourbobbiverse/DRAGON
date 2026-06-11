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

