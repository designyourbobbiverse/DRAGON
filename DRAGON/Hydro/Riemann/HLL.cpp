//
//  HLL.cpp
//  DRAGON/Hydro/Riemann
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

namespace HLL{
PrimitiveState roeAvg(PrimitiveState L, PrimitiveState R){
    double sql = sqrt(L.rho), sqr = sqrt(R.rho);
    double sum = sql + sqr;
    sql /= sum;
    sqr /= sum;
    
    PrimitiveState w;
    w.rho = sql*L.rho + sqr*R.rho;
    w.v = sql*L.v + sqr*R.v;
    w.p = sql*L.p + sqr*R.p;
#ifdef MHD
    w.B = sql*L.B + sqr*R.B;
#endif
    return w;
}
}


//MARK: HLL/E
ConservativeState Riemann::HLL(){
    PrimitiveState M = HLL::roeAvg(L,R);
    //Compute Sound/Fast Speed
#ifdef MHD
    //Set Normal Magnetic Fields
    double Bx = (L.B.x+R.B.x)/2;
    L.B.x = Bx; R.B.x = Bx;
    //Calculate Fast Mode
    double aL = L.c_fast(), aR = R.c_fast(), aM = M.c_fast();
#else
    //Sound Speed
    double aL = L.cs(), aR = R.cs(), aM = M.cs();
#endif
    //Compare to v +- a
    double SL = fmin(L.v.x - aL, M.v.x-aM);
    double SR = fmax(R.v.x + aR, M.v.x+aM);
    
    return HLL(SL, SR);
}
ConservativeState Riemann::HLLE(){
    double sql = sqrt(L.rho), sqr = sqrt(R.rho);
#ifdef MHD
    //Set Normal Magnetic Fields
    double Bx = (L.B.x+R.B.x)/2;
    L.B.x = Bx; R.B.x = Bx;
    //Calculate Fast Mode
    double aL = L.c_fast(), aR = R.c_fast();
#else
    //Sound Speed
    double aL = L.cs(), aR = R.cs();
#endif
    double _v = (R.v.x-L.v.x) / (sql + sqr);
    double d = (sql*aL*aL + sqr*aR*aR) / (sql + sqr)  +  0.5 * (sql * sqr) * _v * _v;
    d = sqrt(d);
    double u = 0.5 * (L.v.x + R.v.x);
    
    //Compare to v +- a
    double SL = fmin(L.v.x - aL, u-d);
    double SR = fmax(R.v.x + aR, u+d);
    
    return HLL(SL, SR);
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



#ifdef HYDRO_AVAILABLE

//MARK: HLLC
ConservativeState Riemann::HLLC(){
    //Roe averages
    PrimitiveState M = HLL::roeAvg(L,R);
    //Compute Sound Speed
    double aL = L.cs(), aR = R.cs(), aM = M.cs();
    //Compare to v +- a
    double SL = fmin(L.v.x - aL, M.v.x-aM);
    double SR = fmax(R.v.x + aR, M.v.x+aM);
    
    return HLLC(SL, SR);
}

ConservativeState Riemann::HLLC(double sl, double sr){
    //Exterior region
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

#elif RIEMANN_DEFAULT == RIEMANN_HLLC
#error HLLC Solver incompatible with MHD
#endif



