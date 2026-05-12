//
//  Roe.cpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 29/04/2026.
//
#include "Riemann.hpp"
#include <math.h>
#include "Constants.h"

//MARK: Eigenvectors
inline ConservativeState _K1(double _vx, double _vy, double _vz, double _H, double _a){
    ConservativeState K = ConservativeState();
    K.rho = 1;
    K.px = _vx - _a;
    K.py = _vy;
    K.pz = _vz;
    K.E = _H - _vx*_a;
    return K;
}
inline ConservativeState _K2(double _vx, double _vy, double _vz, double _V2){
    ConservativeState K = ConservativeState();
    K.rho = 1;
    K.px = _vx;
    K.py = _vy;
    K.pz = _vz;
    K.E = _V2/2;
    return K;
}
inline ConservativeState _K3(double _vy){
    ConservativeState K = ConservativeState();
    K.py = 1;
    K.E = _vy;
    return K;
}
inline ConservativeState _K4(double _vz){
    ConservativeState K = ConservativeState();
    K.pz = 1;
    K.E = _vz;
    return K;
}
inline ConservativeState _K5(double _vx, double _vy, double _vz, double _H, double _a){
    ConservativeState K = ConservativeState();
    K.rho = 1;
    K.px = _vx + _a;
    K.py = _vy;
    K.pz = _vz;
    K.E = _H + _vx*_a;
    return K;
}

//MARK: Solver
ConservativeState Riemann::Roe(){
    double sqL = sqrt(L.rho), sqR = sqrt(R.rho);
    //weighted averages
    double _rho = sqL * sqR;
    double _vx = (sqL * L.vx + sqR * R.vx)/(sqL+sqR);
    double _vy = (sqL * L.vy + sqR * R.vy)/(sqL+sqR);
    double _vz = (sqL * L.vz + sqR * R.vz)/(sqL+sqR);
    double _H = (sqL*L.enthalpy() + sqR*R.enthalpy()) / (sqL + sqR);
    double _V2 = _vx*_vx + _vy*_vy + _vz*_vz;
    double _a = sqrt((_gamma-1) * (_H - _V2/2));
    //eigenvalues + eigenvectors
    double lambda[5] = { _vx - _a, _vx, _vx, _vx,  _vx + _a };
    ConservativeState K[5] = {
        _K1(_vx, _vy, _vz, _H, _a),
        _K2(_vx, _vy, _vz, _V2), _K3(_vy), _K4(_vz),
        _K5(_vx, _vy, _vz, _H, _a)
    };
    //wave strengths
    ConservativeState UL= ConservativeState(L), UR = ConservativeState(R);
    double alpha[5] = {
        ((R.p - L.p) - _rho*_a*(R.vx - L.vx)) / (2*_a*_a),
        (R.rho - L.rho) - (R.p - L.p)/(_a*_a),  _rho * (R.vy - L.vy),  _rho * (R.vz - L.vz),
        ((R.p - L.p) + _rho*_a*(R.vx - L.vx)) / (2*_a*_a)
    };
#ifdef Harten_Hyman
    ConservativeState SL = UL + alpha[0]*K[0];
    double aL = sqrt(_gamma * L.p/L.rho), aSL = sqrt(_gamma*SL.pressure()/SL.rho);
    double lambdaL = L.vx-aL, lambdaLS = SL.px/SL.rho - aSL;
    if(lambdaL < 0 &&  lambdaLS > 0 ) {//Left Rarefaction
        double _lambda = lambdaL * (lambda[0] - lambdaLS)/(lambdaL - lambdaLS);
        return UL.flux(L.vx) + _lambda*alpha[0]*K[0];
    }
    ConservativeState SR = UR - alpha[4]*K[4];
    double aR = sqrt(_gamma * R.p/R.rho), aSR = sqrt(_gamma*SR.pressure()/SR.rho);
    double lambdaR = R.vx+aR, lambdaRS = SR.px/SR.rho + aSR;
    if(lambdaR > 0 &&  lambdaRS  < 0) {//Right Rarefaction
        double _lambda = lambdaR * (lambda[4] - lambdaRS)/(lambdaR - lambdaRS);
        return UR.flux(R.vx) - _lambda*alpha[4]*K[4];
    }
#endif
    //Combine the waves
    ConservativeState F = (UL.flux(L.vx) + UR.flux(R.vx));
    for(int i = 0; i < 5; i++) F -= alpha[i]*fabs(lambda[i])*K[i];
    return F/2;
}
