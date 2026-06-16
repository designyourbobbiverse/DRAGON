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

static int validGhosts(int g){
#ifdef MUSCL_Hancock
    return std::max(g, 2);
#else
    return std::max(g, 1);
#endif
}
//MARK: Array Wrapper
Grid1D::Grid1D(int s_, double dx_, int g_): w(s_, validGhosts(g_)), dx(dx_) { }
PrimitiveState& Grid1D::operator[](int k) { return w[k]; }
const PrimitiveState& Grid1D::operator[](int k) const { return w[k]; }
int Grid1D::getSize() const { return w.getSize(); }
int Grid1D::getGhosts() const { return w.getGhosts(); }
Grid1D::~Grid1D() { w.~ExtendedArray1D(); }


//MARK: 1D Godunov Advance

void Grid1D::advance(double dt, bool check_cfl){
    int size = getSize(), ghosts = getGhosts();
    ExtendedArray1D<PrimitiveState> _L(size,ghosts), _R(size,ghosts);//Buffer Grids
    
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
void Grid1D::god_sweep(double dt, ExtendedArray1D<PrimitiveState>& _L, ExtendedArray1D<PrimitiveState>& _R){
    int size = getSize(), ghosts = getGhosts();
    for(int i=-ghosts+1; i<size+ghosts-1; i++) {
        TVD::MUSCL((*this)[i-1], _L[i], (*this)[i], _R[i], (*this)[i+1], dt/dx);
    }
    //Compute Fluxes
    ConservativeState fL, fR;
    fL = Riemann(_R[-ghosts], _L[-ghosts+1]).flux_X(dt/dx);
    for(int i=-ghosts+1; i<size+ghosts-1; i++) {
        fR = Riemann(_R[i], _L[i+1]).flux_X(dt/dx);
        (*this)[i] += (fL - fR) * (dt/dx);
        fL = fR; //Right flux on this cell must equal Left flux on next cell
    }
}


