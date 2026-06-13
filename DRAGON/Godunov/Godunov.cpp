//
//  Godunov.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 12/05/2026.
//  Implementation based mostly on Toro (2009). https://doi.org/10.1007/b79761
//      Godunov (1959). Mat. Sb. 47(89):3, 271–306. https://www.mathnet.ru/eng/sm4873
//

#include "Grid.hpp"
#include "Config.h"
#include "Riemann.hpp"
#include "Boundary.hpp"
#include "CFL.hpp"
#include "TVD.hpp"
#include <cassert>

//MARK: 1D array

Grid1D::Grid1D(int s_, double dx_, int g_): dx(dx_), size(s_), ghosts(g_) {
#ifdef MUSCL_Hancock
    if(ghosts < 2) ghosts = 2;
#else
    if(ghosts < 1) ghosts = 1;
#endif
    w = new PrimitiveState[size+2*ghosts];
}
PrimitiveState& Grid1D::operator[](int k) {
    assert(k + ghosts >= 0 && k < size+ghosts);
    return w[k+ghosts];
}
const PrimitiveState& Grid1D::operator[](int k) const {
    assert(k + ghosts >= 0 && k < size+ghosts);
    return w[k+ghosts];
}
int Grid1D::getSize() const { return size; }
int Grid1D::getGhosts() const { return ghosts; }
Grid1D::~Grid1D() { delete[] w; }




//MARK: 1D Godunov Advance

void Grid1D::advance(double dt, bool check_cfl){
    Grid1D _L(size,dx,ghosts), _R(size,dx,ghosts);//Buffer Grids
    
    while(dt > Timestep_Tolerance){
        //Apply Boundary Conditions
        boundary.apply(*this);
        //CFL Time Constraint
        double t1 = check_cfl ? std::min(dt,CFL::cfl_time(*this)) : dt;
        if (t1 >= dt) t1 = dt;
        //Execute the Advancemtn
        god_sweep(t1,_L,_R);
        dt -= t1;
    }
}


//MARK: Godunov Sweep
void Grid1D::god_sweep(double dt, Grid1D& _L, Grid1D& _R){
    for(int i=-ghosts+1; i<size+ghosts-1; i++) {
        TVD::MUSCL((*this)[i-1], _L[i], (*this)[i], _R[i], (*this)[i+1], dt/dx);
    }
    //Compute Fluxes
    ConservativeState fL, fR;
    int g = ghosts > 1 ? 1 : 0;
    fL = Riemann(_R[-1-g], _L[-g]).flux_X(dt/dx);
    for(int i=-g; i<size+g; i++) {
        fR = Riemann(_R[i], _L[i+1]).flux_X(dt/dx);
        (*this)[i] += (fL - fR) * (dt/dx);
        fL = fR; //Right flux on this cell must equal Left flux on next cell
    }
}


