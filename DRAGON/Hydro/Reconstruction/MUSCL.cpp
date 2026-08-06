//
//  MUSCL.cpp
//  DRAGON/Hydro/Reconstruction
//
//  Created by Bobbie Markwick on 10/06/2026.
//  Implementation based mostly on Toro (2009). https://doi.org/10.1007/b79761
//      van Leer (1979). https://doi.org/10.1016/0021-9991(79)90145-1
//

#include "TVD.hpp"

#include "Config.h"
#include "Constants.h" //For _1_4pi
using namespace DRAGON;

//MARK: MUSCL-Hancock
void TVD::MUSCL(const PrimitiveState& wL, PrimitiveState& _L, const PrimitiveState& wC, PrimitiveState& _R, const PrimitiveState& wR, double dt_dL, const vec3& dB){
#ifdef MUSCL_Hancock
    // Reconstruct a limited primitive-variable slope
    PrimitiveState leftDiff, rightDiff;
    leftDiff.rho = wC.rho - wL.rho; rightDiff.rho = wR.rho - wC.rho;
    leftDiff.v = wC.v - wL.v; rightDiff.v = wR.v - wC.v;
    leftDiff.p = wC.p - wL.p; rightDiff.p = wR.p - wC.p;
    #ifdef MHD
    leftDiff.B = wC.B - wL.B; rightDiff.B = wR.B - wC.B;
    #endif
    PrimitiveState dW = TVD::limit(leftDiff, rightDiff);
    //Spatial half step using the limited slope
    _L.rho = wC.rho - 0.5*dW.rho;
    _R.rho = wC.rho + 0.5*dW.rho;
    _L.v = wC.v - 0.5*dW.v;
    _R.v = wC.v + 0.5*dW.v;
    _L.p = wC.p - 0.5*dW.p;
    _R.p = wC.p + 0.5*dW.p;
    #ifdef MHD
    _L.B = wC.B - 0.5*dW.B;
    _R.B = wC.B + 0.5*dW.B;
    #endif
    
    //Time half step (MUSCL-Hancock Predictor)
    ConservativeState UL = ConservativeState(_L), UR = ConservativeState(_R);
    ConservativeState correction = (UR.flux(_R.v) - UL.flux(_L.v)) * (0.5 * dt_dL);
    #ifdef MHD//MHD Source Terms
    //dB.x = (Bx[x+1/2 CT face] - Bx[x-1/2 CT face]) * dt_dL, etc
    const double dBy_ = 0.5 * wC.v.y * (TVD::minmod(dB.x, -dB.y) - dW.B.x * dt_dL);
    const double dBz_ = 0.5 * wC.v.z * (TVD::minmod(dB.x, -dB.z) - dW.B.x * dt_dL);

    correction.B.y -= dBy_;
    correction.B.z -= dBz_;
    correction.E   -= _1_4pi * (wC.B.y*dBy_ + wC.B.z*dBz_);   // hold thermal p fixed
    #endif
    // Correction is applied in conservative, then automatically converted back to primitive
    _L = UL - correction;
    _R = UR - correction;
    
    //Check to make sure this didn't overshoot past positivity limit, fallback to First order if needed
    if(_L.isPhysical() && _R.isPhysical() ) return;
#endif
    //First order version if MUSCL is disabled or failed physicality check.
    _L = wC;
    _R = wC;
}
