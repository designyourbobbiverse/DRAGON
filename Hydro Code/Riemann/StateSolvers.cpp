//
//  OtherSolvers.cpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 29/04/2026.
//
#include "Riemann.hpp"
#include <math.h>
#include "Constants.h"




//MARK: Primative Variable
//Primative Variable Solver
RiemannSolution Riemann::PVRS(double aL, double aR, double p_pvrs){
    RiemannSolution s = RiemannSolution(*this);

    if (p_pvrs==0) p_pvrs = (L.p + R.p)/2 + (L.rho + R.rho)*(L.vx - R.vx)*(aL + aR)/8;
    //Pressure
    s.sL.p = p_pvrs;
        s.sR.p = p_pvrs;
    //Veloctiy
    s.sL.vx = (L.vx + R.vx)/2 + 2*(L.p - R.p) / ((L.rho + R.rho)*(aL + aR));
        s.sR.vx = s.sL.vx;
    //Density
    s.sL.rho = L.rho + (L.rho+R.rho) * ((L.vx - s.sL.vx) / (aL + aR));
    s.sR.rho = R.rho - (L.rho+R.rho) * ((R.vx - s.sR.vx) / (aL + aR));
    
    return s;
}
RiemannSolution Riemann::PVRS(){
    double aL = sqrt(gamma * L.p/L.rho), aR = sqrt(gamma * R.p/L.rho); // Sound Speeds
    return PVRS(aL, aR, 0);
}

//MARK: Two-Rarefaction
RiemannSolution Riemann::TRRS(double aL, double aR){
    RiemannSolution s = RiemannSolution(*this);

    double _LR = pow(L.p / R.p, _Ginv) * aR/aL;
    //Velocity
    s.sL.vx = (_LR*L.vx + R.vx + _2_Gm1*(_LR*aL-aR)) / (_LR + 1);
        s.sR.vx = s.sL.vx;
    //Pressure
    s.sL.p = L.p * pow(1 - _Gm1_2*(s.sL.vx - L.vx)/aL, _2G_Gm1);
        s.sR.p = s.sL.p;
    //Density
    s.sL.rho = L.rho * pow(s.sL.p / L.p , _Ginv);
    s.sR.rho = R.rho * pow(s.sR.p / R.p , _Ginv);
    
    return s;
}
RiemannSolution Riemann::TRRS(){
    double aL = sqrt(gamma * L.p/L.rho), aR = sqrt(gamma * R.p/L.rho); // Sound Speeds
    return TRRS(aL, aR);
}

//MARK: Two-Shock
RiemannSolution Riemann::TSRS(double aL, double aR, double pGuess){
    RiemannSolution s = RiemannSolution(*this);

    double gL = sqrt((_2_Gp1 / L.rho) / (pGuess + _Gm1_Gp1 * L.p) );
    double gR = sqrt((_2_Gp1 / R.rho) / (pGuess + _Gm1_Gp1 * R.p) );
    //Pressure
    s.sL.p = (gL*L.p + gR*R.p + (L.vx - R.vx)) / (gL + gR);
        s.sR.p = s.sL.p;
    //Velocities
    s.sL.vx = (L.vx + R.vx)/2 + (gL*(L.p - s.sL.p) - gR*(R.p - s.sR.p))/2;
        s.sR.vx = s.sL.vx;
    //Density
    s.sL.rho = L.rho * (s.sL.p + _Gm1_Gp1*L.p)/(_Gm1_Gp1*s.sL.p + L.p);
    s.sR.rho = R.rho * (s.sR.p + _Gm1_Gp1*R.p)/(_Gm1_Gp1*s.sR.p + R.p);

    return s;
}
RiemannSolution Riemann::TSRS(){
    double aL = sqrt(gamma * L.p/L.rho), aR = sqrt(gamma * R.p/L.rho); // Sound Speeds
    double p_pvrs = (L.p + R.p)/2 + (L.vx - R.vx) * (L.rho + R.rho) * (aL + aR) / 8;
    return TSRS(aL, aR, fmax(p_pvrs,0));
}


//MARK: Adaptive Solvers
//PVRS, pivot to iterative exact if conditions not met
RiemannSolution Riemann::PVRS_Iter(){
    double aL = sqrt(gamma * L.p/L.rho), aR = sqrt(gamma * R.p/L.rho); // Sound Speeds
    double p_pvrs = (L.p + R.p)/2 + (L.vx - R.vx) * (L.rho + R.rho) * (aL + aR) / 8;
    double pMin = fmin(L.p, R.p), pMax = fmax(L.p, R.p);
    
    if (pMin < p_pvrs && p_pvrs < pMax && pMax / pMin < Adaptive_PVRS_Ratio) { //PVRS
        return PVRS(aL, aR, p_pvrs);
    } else { //Pivot to Exact
        return exact(p_pvrs);
    }
}
//PVRS, pivot to TRRS/TSRS if conditions not met
RiemannSolution Riemann::PVRS_TXRS(){
    double aL = sqrt(gamma * L.p/L.rho), aR = sqrt(gamma * R.p/L.rho); // Sound Speeds
    double rhoBar = (L.rho + R.rho)/2, aBar = (aL + aR)/2;
    double p_pvrs = (L.p + R.p)/2 + (L.vx - R.vx)/2 * rhoBar*aBar;
    double pMin = fmin(L.p, R.p), pMax = fmax(L.p, R.p);
    
    if (pMin < p_pvrs && p_pvrs < pMax && pMax / pMin < Adaptive_PVRS_Ratio) { //PVRS
        return PVRS(aL, aR, p_pvrs);
    } else if (p_pvrs < pMin) { //Two Rarefactions
        return TRRS(aL, aR);
    } else { //Two Shocks
        return TSRS(aL, aR, p_pvrs);
    }
}
    

//TRRS, pivot to iterative exact if conditions not met
RiemannSolution Riemann::TRRS_Iter(){
    return TRRS_Iter( (L.p + R.p)/2 );
}
RiemannSolution Riemann::TRRS_Iter(double pGuess){
    double p = L.p < R.p ? f(L.p, L) : f(R.p, R);
    if(p > 0 ){
        return TRRS();
    } else {
        return exact(pGuess);
    }
}
