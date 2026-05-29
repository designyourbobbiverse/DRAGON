//
//  Riemann/Exact.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#include "Riemann.hpp"
#include <math.h>
#include "Constants.h"




//MARK: Setup
Riemann::Riemann(PrimitiveState _L, PrimitiveState _R){ L = _L; R = _R; }
RiemannSolution::RiemannSolution(Riemann problem){
    wL = problem.L;
    wR = problem.R;
    //Copy Transverse Values
    sL.vy = wL.vy;
    sR.vy = wR.vy;
    sL.vz = wL.vz;
    sR.vz = wR.vz;
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
    RiemannSolution s = RiemannSolution(*this);
    //Pressure
    s.sL.p = exact_StarP(pGuess);
        s.sR.p = s.sL.p;
    //velocity
    s.sL.vx = exact_StarV(s.sL.p);
        s.sR.vx = s.sL.vx;
    //Density
    s.sL.rho = exact_StarRho(L, s.sL.p);
    s.sR.rho = exact_StarRho(R, s.sR.p);
    
    return s;
}



double Riemann::exact_StarP(double pGuess){
    double pStar = pGuess, CHA = 1; int iters = ExactSolver_MaxIters;
    do{
        double fp = f(pStar, L) + f(pStar, R) + R.vx - L.vx;
        double dfdp = df(pStar, L) + df(pStar, R);
        double delta = fmin(fp/dfdp, 0.8*pStar);
        pStar -= delta;
        CHA = fabs(delta/(pStar+delta/2));
    } while(CHA > ExactSolver_Tolerance &&  --iters != 0 );
    return pStar;
}
//v* and rho* given p*
double Riemann::exact_StarV(double pStar){
    return (L.vx + R.vx + f(pStar,R)-f(pStar,L))/2.0;
}
double Riemann::exact_StarRho(PrimitiveState w, double p){
    if(p>w.p) return w.rho * (p+_Gm1_Gp1*w.p)/(_Gm1_Gp1*p+w.p);//Shock
    else return w.rho * pow( p/w.p, _Ginv); //Rarefaction
}



