//
//  Godunov.cpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#include "Godunov.hpp"
#include "Constants.h"
#include "Config.h"
#include "Riemann.hpp"
#include "Boundary.hpp"
#include <math.h>
#include <cassert>

//MARK: 1D array

Grid1D::Grid1D(int size, double dx_, int ghosts): dx(dx_) {
    w = new PrimitiveState[size+2*ghosts];
    this->size = size; this->ghosts = ghosts;
    
}
PrimitiveState& Grid1D::operator[](int k) {
    assert(k + ghosts >= 0 && k < size+ghosts);
    return w[k+ghosts];
}
const PrimitiveState& Grid1D::operator[](int k) const {
    assert(k + ghosts >= 0 && k < size+ghosts);
    return w[k+ghosts];
}
int Grid1D::getSize(){ return size; }
int Grid1D::getGhosts(){ return ghosts; }
Grid1D::~Grid1D(){ delete[] w; }



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
void Grid1D::advance(double dt){
    Grid1D _L(size,ghosts), _R(size,ghosts);
    
    if(boundary) boundary->apply(*this);
    god_sweep(dt,_L,_R);
}
void Grid1D::god_sweep(double dt, Grid1D& _L, Grid1D& _R){
    for(int i=-ghosts; i<size+ghosts; i++) {
#ifdef MUSCL_Hancock
        if(i==-ghosts || i == size+ghosts-1){
            _L[i] = (*this)[i];
            _R[i] = (*this)[i];
            continue;
        }
        PrimitiveState dW = minmod((*this)[i] - (*this)[i-1], (*this)[i+1] - (*this)[i]);
        _L[i] = (*this)[i] - dW/2;
        _R[i] = (*this)[i] + dW/2;
        ConservativeState UL = ConservativeState(_L[i]), UR = ConservativeState(_R[i]);
        ConservativeState correction = (UR.flux() - UL.flux()) * (0.5 * dt/dx);
        _L[i] = PrimitiveState(UL - correction);
        _R[i] = PrimitiveState(UR - correction);
        
        if(!_L[i].isPhysical() || !_R[i].isPhysical() ) //Check to make sure this is physical, fallback to First order if not
#endif
        {
            _L[i] = (*this)[i];
            _R[i] = (*this)[i];
        }
    }
    //Compute Fluxes
    ConservativeState fL = (ghosts > 0 ? _R[-1] : _L[0]), fR;
    for(int i=0; i<size; i++) {
        fR = Riemann(_R[i], _L[i+1]).flux();
        (*this)[i] += (fL - fR) * (dt/dx);
        fL = fR; //Right flux on this cell must equal Left flux on next cell
    }
}


