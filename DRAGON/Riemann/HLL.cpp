//
//  Riemann/HLL.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 29/04/2026.
//  Implementation based mostly on Toro (2009). https://doi.org/10.1007/b79761
//      Harten, Lax, & van Leer (1983). https://doi.org/10.1137/1025002
//      Toro, Spruce, & Speares (1994). https://doi.org/10.1007/BF01414629
//      Einfeldt (1988). https://doi.org/10.1137/0725021
//



#include "Riemann.hpp"
#include <math.h>
#include "Constants.h"



//MARK: HLL/E
ConservativeState Riemann::HLL(){
    double sql = sqrt(L.rho), sqr = sqrt(R.rho);
    PrimitiveState M = (sql*L + sqr*R) / (sql + sqr);
#ifdef MHD
    double aL = L.c_fast(), aR = R.c_fast(), aM = M.c_fast();
#else
    double aL = L.cs(), aR = R.cs(), aM = M.cs();
#endif
    double SL = fmin(L.v.x - aL, M.v.x-aM);
    double SR = fmin(R.v.x + aR, M.v.x+aM);
    
    return HLL(SL, SR);
}
ConservativeState Riemann::HLLE(){
    double sql = sqrt(L.rho), sqr = sqrt(R.rho);
#ifdef MHD
    double aL = L.c_fast(), aR = R.c_fast();
#else
    double aL = L.cs(), aR = R.cs();
#endif
    double _v = (R.v.x-L.v.x) / (sql + sqr);
    double d = (sql*aL*aL + sqr*aR*aR) / (sql + sqr)  +  0.5 * (sql * sqr) * _v * _v;
    d = sqrt(d);
    double u = 0.5*(R.v.x+L.v.x);
    
    return HLL(u-d, u+d);
}
ConservativeState Riemann::HLL(double sl, double sr){
    //Outside region
    if(sl >= 0) return ConservativeState(L).flux(L.v);
    if(sr <= 0) return ConservativeState(R).flux(R.v);
    
    ConservativeState UL = ConservativeState(L);
    ConservativeState UR = ConservativeState(R);
    ConservativeState FL = UL.flux(L.v);
    ConservativeState FR = UR.flux(R.v);

    return (sr*FL - sl*FR + sl*sr*(UR-UL))/(sr-sl);
}



#if !defined(MHD) || defined(TESTMODE)

//MARK: HLLC
ConservativeState Riemann::HLLC(){
    //Estimate the pressure
    double aL = sqrt(_gamma * L.p/L.rho), aR = sqrt(_gamma * R.p/R.rho); // Sound Speeds
    double p_pvrs = fmax(0,(L.p + R.p)/2 + (L.rho + R.rho)*(L.v.x - R.v.x)*(aL + aR)/8);
    //Left Speed
    double SL = aL;
    if(p_pvrs > L.p) SL *= sqrt(1 + _Gp1_2G*(p_pvrs/L.p - 1));
    SL = L.v.x - SL;
    if(SL > 0) return L.flux();
    //Right Speed
    double SR = aR;
    if(p_pvrs > R.p) SR *= sqrt(1 + _Gp1_2G*(p_pvrs/R.p - 1));
    SR = R.v.x + SR;
    if(SR < 0) return R.flux();

    
    return HLLC(SL, SR);
}

ConservativeState Riemann::HLLC(double sl, double sr){
    //Outside region
    if(sl >= 0) return L.flux();
    if(sr <= 0) return R.flux();
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
    return UX.flux(X.v) + (U - UX) * sx;
}

#elif RIEMANN_DEFAULT_HYDRO == RIEMANN_HLLC
#error HLLC Solver incompatible with MHD
#endif



