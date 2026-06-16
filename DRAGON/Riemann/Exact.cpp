//
//  Riemann/Exact.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 28/04/2026.
//  Implementation based mostly on Toro (2009). https://doi.org/10.1007/b79761
//

#include "Riemann.hpp"
#include <math.h>
#include "Constants.h"




//MARK: Setup
Riemann::Riemann(PrimitiveState _L, PrimitiveState _R){ L = _L; R = _R; }
RiemannSolution::RiemannSolution(Riemann problem){
    wL = problem.L;
    wR = problem.R;
    //Cop.y Transverse Values
    sL.v = wL.v;
    sR.v = wR.v;
}


//MARK: Velocity Jump Function
double Riemann::f(double p, PrimitiveState w){
    if (p > w.p) { //Shock
        double A= _2_Gp1/w.rho, B= _Gm1_Gp1 * w.p;
        return (p - w.p) * sqrt(A/(B+p));
    } else {//Rarefaction
        double a = sqrt(_gamma * w.p / w.rho);
        return _2_Gm1 * a * (pow(p/w.p, _Gm1_2G)-1);
    }
}
double df(double p, PrimitiveState w){
    if (p > w.p) { //Shock
        double A = _2_Gp1/w.rho, B= _Gm1_Gp1 * w.p;
        return sqrt(A/(p+B)) * (1-(p-w.p)/(2*(B+p)));
    }
    else {//Rarefaction
        double a = sqrt(_gamma * w.p / w.rho);
        return  pow(p/w.p, -_Gp1_2G) / (w.rho * a);
    }
}

//MARK: Exact Riemann Solver
RiemannSolution Riemann::exact(){
    return exact( (L.p + R.p)/2.0 );
}
RiemannSolution Riemann::exact(double pGuess){
    #ifdef Exact_Rarefactions_Check
    double p_min = fmin(L.p, R.p);
    if(f(p_min,L) + f(p_min,R) + R.v.x - L.v.x >= 0) return TRRS();
    #endif
    
    RiemannSolution s = RiemannSolution(*this);
    //Pressure
    s.sL.p = exact_StarP(pGuess);
        s.sR.p = s.sL.p;
    //velocity
    s.sL.v.x = exact_StarV(s.sL.p);
        s.sR.v.x = s.sL.v.x;
    //Density
    s.sL.rho = exact_StarRho(L, s.sL.p);
    s.sR.rho = exact_StarRho(R, s.sR.p);
    
    return s;
}



double Riemann::exact_StarP(double pGuess){
    double pStar = pGuess, CHA = 1; int iters = ExactRiemann_MaxIters;
    do{
        double fp = f(pStar, L) + f(pStar, R) + R.v.x - L.v.x;
        double dfdp = df(pStar, L) + df(pStar, R);
        double delta = fmin(fp/dfdp, 0.8*pStar);
        pStar -= delta;
        CHA = fabs(delta/(pStar+delta/2));
    } while(CHA > ExactRiemann_Tolerance &&  --iters != 0 );
    return pStar;
}
//v* and rho* given p*
double Riemann::exact_StarV(double pStar){
    return (L.v.x + R.v.x + f(pStar,R)-f(pStar,L))/2.0;
}
double Riemann::exact_StarRho(PrimitiveState w, double p){
    if(p>w.p) return w.rho * (p+_Gm1_Gp1*w.p)/(_Gm1_Gp1*p+w.p);//Shock
    else return w.rho * pow( p/w.p, _Ginv); //Rarefaction
}


//MARK: Two-Rarefaction
RiemannSolution Riemann::TRRS(){
    double aL = sqrt(_gamma * L.p/L.rho), aR = sqrt(_gamma * R.p/R.rho); // Sound Speeds
    return TRRS(aL, aR);
}
RiemannSolution Riemann::TRRS(double aL, double aR){
    RiemannSolution s = RiemannSolution(*this);

    double _LR = pow(L.p / R.p, _Ginv) * aR/aL;
    //Velocity
    s.sL.v.x = (_LR*L.v.x + R.v.x + _2_Gm1*(_LR*aL-aR)) / (_LR + 1);
        s.sR.v.x = s.sL.v.x;
    //Pressure
    s.sL.p = L.p * pow(1 - _Gm1_2*(s.sL.v.x - L.v.x)/aL, _2G_Gm1);
        s.sR.p = s.sL.p;
    //Density
    s.sL.rho = L.rho * pow(s.sL.p / L.p , _Ginv);
    s.sR.rho = R.rho * pow(s.sR.p / R.p , _Ginv);
    
    return s;
}
