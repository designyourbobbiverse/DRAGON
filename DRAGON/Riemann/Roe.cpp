//
//  Riemann/Roe.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 29/04/2026.
//  Implementation based mostly on Toro (2009). https://doi.org/10.1007/b79761
//      Roe (1981). https://doi.org/10.1016/0021-9991(81)90128-5
//      Roe & Pike (1984). Computing Methods in Applied Sciences and Engineering VI, pp. 499–518.
//      Harten and Hyman (1983). https://doi.org/10.1016/0021-9991(83)90066-9
//

#include "Riemann.hpp"
#include <math.h>
#include "Constants.h"
#include "Config.h"

#ifdef HYDRO_AVAILABLE
//MARK: Eigenvectors
inline ConservativeState _K1(vec3 _v, double _H, double _a){
    ConservativeState K = ConservativeState();
    K.rho = 1;
    K.p = _v - vec3{_a,0,0};
    K.E = _H - _v.x*_a;
    return K;
}
inline ConservativeState _K2(vec3 _v, double _V2){
    ConservativeState K = ConservativeState();
    K.rho = 1;
    K.p = _v;
    K.E = _V2/2;
    return K;
}
inline ConservativeState _K3(double _vy){
    ConservativeState K = ConservativeState();
    K.p.y = 1;
    K.E = _vy;
    return K;
}
inline ConservativeState _K4(double _vz){
    ConservativeState K = ConservativeState();
    K.p.z = 1;
    K.E = _vz;
    return K;
}
inline ConservativeState _K5(vec3 _v, double _H, double _a){
    ConservativeState K = ConservativeState();
    K.rho = 1;
    K.p = _v + vec3{_a,0,0};
    K.E = _H + _v.x*_a;
    return K;
}

//MARK: Solver
ConservativeState Riemann::Roe(){
    double sqL = sqrt(L.rho), sqR = sqrt(R.rho);
    //weighted averages
    double _rho = sqL * sqR;
    vec3 _v = (sqL * L.v + sqR * R.v)/(sqL+sqR);
    double _H = (sqL*L.enthalpy() + sqR*R.enthalpy()) / (sqL + sqR);
    double _V2 = _v.x*_v.x + _v.y*_v.y + _v.z*_v.z;
    double _a = sqrt((_gamma-1) * (_H - _V2/2));
    //eigenvalues + eigenvectors
    double lambda[5] = { _v.x - _a, _v.x, _v.x, _v.x,  _v.x + _a };
    ConservativeState K[5] = {
        _K1(_v, _H, _a),
        _K2(_v, _V2), _K3(_v.y), _K4(_v.z),
        _K5(_v, _H, _a)
    };
    //wave strengths
    ConservativeState UL= ConservativeState(L), UR = ConservativeState(R);
    double alpha[5] = {
        ((R.p - L.p) - _rho*_a*(R.v.x - L.v.x)) / (2*_a*_a),
        (R.rho - L.rho) - (R.p - L.p)/(_a*_a),  _rho * (R.v.y - L.v.y),  _rho * (R.v.z - L.v.z),
        ((R.p - L.p) + _rho*_a*(R.v.x - L.v.x)) / (2*_a*_a)
    };
#ifdef Harten_Hyman
    ConservativeState SL = UL + alpha[0]*K[0];
    double aL = sqrt(_gamma * L.p/L.rho), aSL = sqrt(_gamma*SL.pressure()/SL.rho);
    double lambdaL = L.v.x-aL, lambdaLS = SL.p.x/SL.rho - aSL;
    if(lambdaL < 0 &&  lambdaLS > 0 ) {//Left Rarefaction
        double _lambda = lambdaL * (lambda[0] - lambdaLS)/(lambdaL - lambdaLS);
        return UL.flux(L.v) + _lambda*alpha[0]*K[0];
    }
    ConservativeState SR = UR - alpha[4]*K[4];
    double aR = sqrt(_gamma * R.p/R.rho), aSR = sqrt(_gamma*SR.pressure()/SR.rho);
    double lambdaR = R.v.x+aR, lambdaRS = SR.p.x/SR.rho + aSR;
    if(lambdaR > 0 &&  lambdaRS  < 0) {//Right Rarefaction
        double _lambda = lambdaR * (lambda[4] - lambdaRS)/(lambdaR - lambdaRS);
        return UR.flux(R.v) - _lambda*alpha[4]*K[4];
    }
#endif
    //Combine the waves
    ConservativeState F = (UL.flux(L.v) + UR.flux(R.v));
    for(int i = 0; i < 5; i++) F -= alpha[i]*fabs(lambda[i])*K[i];
    return F/2;
}

#elif RIEMANN_DEFAULT == RIEMANN_ROE
#error Roe Solver incompatible with MHD
#endif
