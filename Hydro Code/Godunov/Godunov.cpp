//
//  Godunov.cpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#include "Godunov.hpp"
#include "Constants.h"
#include "Riemann.hpp"
#include <math.h>
#include <cassert>

//MARK: 1D array

FluidLine::FluidLine(int size, int ghosts){
    w = new PrimitiveState[size+2*ghosts];
    this->size = size; this->ghosts = ghosts;
}
PrimitiveState& FluidLine::operator[](int k) {
    assert(k + ghosts >= 0 && k < size+ghosts);
    return w[k+ghosts];
}
const PrimitiveState& FluidLine::operator[](int k) const {
    assert(k + ghosts >= 0 && k < size+ghosts);
    return w[k+ghosts];
}

FluidLine::~FluidLine(){ delete w; }



//MARK: TVD

double minmod(double a, double b) {
    if(a*b <= 0) return 0;
    if(fabs(a) < fabs(b)) return a;
    return a < 0 ? -fabs(b) : fabs(b);
}
PrimitiveState minmod(const PrimitiveState& a, const PrimitiveState& b) {
    PrimitiveState W;
    W.rho = minmod(a.rho,b.rho);
    W.vx = minmod(a.vx,b.vx);
    W.vy = minmod(a.vy,b.vy);
    W.vz = minmod(a.vz,b.vz);
    W.p = minmod(a.p,b.p);
    return W;
}


//MARK: 1D Godunov Advance
void FluidLine::advance(double dt, double dx){
    FluidLine _L(size,ghosts), _R(size,ghosts);
    advance(dt, dx,_L,_R);
}
void FluidLine::advance(double dt, double dx, FluidLine& _L, FluidLine& _R){
    for(int i=0; i<size; i++) {
#ifdef MUSCL_Hancock
        PrimitiveState dW = minmod((*this)[i] - (*this)[i-1], (*this)[i+1] - (*this)[i]);
        _L[i] = (*this)[i] - dW/2;
        _R[i] = (*this)[i] + dW/2;
        ConservativeState UL = ConservativeState(_L[i]), UR = ConservativeState(_R[i]);
        ConservativeState correction = (UR.flux() - UL.flux()) * 0.5 * (dt/stepSize);
        _L[i] = PrimitiveState(UL - correction);
        _R[i] = PrimitiveState(UR - correction);
        
        if(_L[i].isPhysical() && _R[i].isPhysical() ) //Check to make sure this is physical, fallback to First order if not
#endif
        {
            _L[i] = (*this)[i];
            _R[i] = (*this)[i];
        }
    }
    //Compute Fluxes
    ConservativeState fL, fR;//Initialize fL based on left boundary condition
    for(int i=0; i<size; i++){
        fR = Riemann(_R[i], _L[i+1]).flux();
        (*this)[i] += (fL - fR) * (dt/dx);
        fL = fR; //Right flux on this cell must equal Left flux on next cell
    }
}


