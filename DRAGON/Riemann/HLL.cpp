//
//  Riemann/HLL.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 29/04/2026.
//  Implementation based mostly on Toro (2009). https://doi.org/10.1007/b79761
//      Harten, Lax, & van Leer (1983). https://doi.org/10.1137/1025002
//      Toro, Spruce, & Speares (1994). https://doi.org/10.1007/BF01414629
//



#include "Riemann.hpp"
#include <math.h>
#include "Constants.h"


//MARK: HLL
ConservativeState Riemann::HLL(){
    //Estimate the pressure
    double aL = sqrt(_gamma * L.p/L.rho), aR = sqrt(_gamma * R.p/R.rho); // Sound Speeds
    double p_pvrs = (L.p + R.p)/2 + (L.rho + R.rho)*(L.v.x - R.v.x)*(aL + aR)/8;
    //Left Speed
    double SL = aL;
    if(p_pvrs > L.p) SL *= sqrt(1 + _Gp1_2G*(p_pvrs/L.p - 1));
    SL = L.v.x - SL;
    //Right Speed
    double SR = aR;
    if(p_pvrs > R.p) SR *= sqrt(1 + _Gp1_2G*(p_pvrs/R.p - 1));
    SR = R.v.x + SR;
    
    return HLL(SL, SR);
}

ConservativeState Riemann::HLL(double sl, double sr){
    //Outside region
    if(sl >= 0) return ConservativeState(L).flux(L.v.x);
    if(sr <= 0) return ConservativeState(R).flux(R.v.x);
    
    ConservativeState UL = ConservativeState(L);
    ConservativeState UR = ConservativeState(R);
    ConservativeState FL = UL.flux(L.v.x);
    ConservativeState FR = UR.flux(R.v.x);

    return (sr*FL - sl*FR + sl*sr*(UR-UL))/(sr-sl);
}


//MARK: HLLC
ConservativeState Riemann::HLLC(){
    //Estimate the pressure
    double aL = sqrt(_gamma * L.p/L.rho), aR = sqrt(_gamma * R.p/R.rho); // Sound Speeds
    double p_pvrs = (L.p + R.p)/2 + (L.rho + R.rho)*(L.v.x - R.v.x)*(aL + aR)/8;
    //Left Speed
    double SL = aL;
    if(p_pvrs > L.p) SL *= sqrt(1 + _Gp1_2G*(p_pvrs/L.p - 1));
    SL = L.v.x - SL;
    //Right Speed
    double SR = aR;
    if(p_pvrs > R.p) SR *= sqrt(1 + _Gp1_2G*(p_pvrs/R.p - 1));
    SR = R.v.x + SR;
    
    return HLLC(SL, SR);
}

ConservativeState Riemann::HLLC(double sl, double sr){
    //Outside region
    if(sl >= 0) return ConservativeState(L).flux(L.v.x);
    if(sr <= 0) return ConservativeState(R).flux(R.v.x);
    //Calculate the contact wave
    double _pr = R.p + R.rho*R.v.x*(R.v.x - sr), _pl = L.p + L.rho*L.v.x*(L.v.x - sl);
    double sc = (_pr - _pl) / (R.rho*(R.v.x-sr) - L.rho*(L.v.x-sl));
    //Left or Right
    PrimitiveState X = (sc > 0 ? L : R);
    double sx = sc > 0 ? sl : sr;
    ConservativeState UX = ConservativeState(X);
    //Compute Star Region
    ConservativeState U = UX * (sx - X.v.x) / (sx - sc);
    U.p.x = U.rho * sc;
    U.E +=  U.rho*(sc - X.v.x)*(sc + X.p/(X.rho*(sx-X.v.x)) );
    //Compute Flux
    return UX.flux(X.v.x) + (U - UX) * sx;
}


