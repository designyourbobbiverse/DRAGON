//
//  TVD.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 10/06/2026.
//  Implementation based mostly on Toro (2009). https://doi.org/10.1007/b79761
//      van Leer (1979). https://doi.org/10.1016/0021-9991(79)90145-1
//

#include "TVD.hpp"
#include <math.h>
#include "Config.h"

double TVD::minmod(double a, double b) {
    if(a*b <= 0) return 0;
    if(fabs(a) < fabs(b)) return a;
    return a < 0 ? -fabs(b) : fabs(b);
}
PrimitiveState TVD::minmod(const PrimitiveState& a, const PrimitiveState& b) {
    PrimitiveState W;
    W.rho = minmod(a.rho,b.rho);
    W.vx = minmod(a.vx,b.vx);
    W.vy = minmod(a.vy,b.vy);
    W.vz = minmod(a.vz,b.vz);
    W.p = minmod(a.p,b.p);
    return W;
}

void TVD::MUSCL(const PrimitiveState& wL, PrimitiveState& _L, const PrimitiveState& wC, PrimitiveState& _R, const PrimitiveState& wR, double dt_dL){
#ifdef MUSCL_Hancock
    //Compute State
    PrimitiveState dW = TVD::minmod(wC - wL, wR - wC);
    _L = wC - dW/2;
    _R = wC + dW/2;
    ConservativeState UL = ConservativeState(_L), UR = ConservativeState(_R);
    ConservativeState correction = (UR.flux() - UL.flux()) * (0.5 * dt_dL);
    _L = UL - correction;
    _R = UR - correction;
    
    if(_L.isPhysical() && _R.isPhysical() ) return; //Check to make sure this is physical, fallback to First order if not
#endif
    _L = wC;
    _R = wC;
}

