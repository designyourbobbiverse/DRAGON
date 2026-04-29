//
//  Riemann.cpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#include "Riemann.hpp"
#include <math.h>
#include "Constants.h"




//MARK: Setup
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
double Riemann::f(double p, PrimativeState w){
    if (p > w.p) { //Shock
        double A= _2_Gp1/w.rho, B= _Gm1_Gp1 * w.p;
        return (p-w.p) * sqrt(A/(B+p));
    } else {//Rarefaction
        double a = sqrt(gamma * w.p / w.rho);
        return _2_Gm1 * a * (pow(p/w.p, _Gm1_2G)-1);
    }
}
double df(double p, PrimativeState w){
    if (p > w.p) { //Shock
        double A = _2_Gp1/w.rho, B= _Gm1_Gp1 * w.p;
        return sqrt(A/(p+B)) * (1-(p-w.p)/(2*(B+p)));
    }
    else {//Rarefaction
        double a = sqrt(gamma * w.p / w.rho);
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
double Riemann::exact_StarRho(PrimativeState w, double p){
    if(p>w.p) return w.rho * (p+_Gm1_Gp1*w.p)/(_Gm1_Gp1*p+w.p);//Shock
    else return w.rho * pow( p/w.p, _Ginv); //Rarefaction
}



//MARK: Solution Sampling
PrimativeState RiemannSolution::sample(double x_t){
    PrimativeState state;
    //Handle left vs Right side
    bool isLeft = x_t < sR.vx;
    if(isLeft){ mirror(); x_t=-x_t; }
    
    double a = sqrt(gamma * wR.p / wR.rho);
    
    //Determine Zone
    int zone = 0; // 1 = outside, 2 = fan, 3 = star
    if (sR.p > wR.p){ //shock
        double scale = sqrt(_Gp1_2G*sR.p/wR.p + _Gm1_2G );
        zone = ( (x_t-wR.vx)  > scale * a ) ? 1 : 3;
    } else{ //Rarefraction
        if ( (x_t-wR.vx) > a ) zone = 1;
        else if ( (x_t-sR.vx) > a ) zone = 2;
        else if ( (x_t-sR.vx) > a *pow(sR.p/wR.p, _Gm1_2G) ) zone = 2;
        else zone = 3;
    }
    //Calculate the State at x/t
    switch(zone){
        case 1: state=wR; break; //Outer Region
        case 3: state=sR;break; //Star Region
        case 2://Fan
            state = wR;
            double scale = _2_Gp1 - _Gm1_Gp1 * (wR.vx-x_t)/a;
            state.rho *= pow(scale, _2_Gm1);
            state.p *= pow(scale,_2G_Gm1);
            state.vx = _2_Gp1 * (x_t-a + wR.vx * _Gm1_2);
            break;
    }
    //Be a good citizen, restore original state if we mirrored
    if(isLeft) mirror();

    return state;
}
ConservativeState RiemannSolution::flux(){ return flux(0); }
ConservativeState RiemannSolution::flux(double x_t){
    PrimativeState w = sample(x_t);
    return ConservativeState(w).flux(w.vx);
}

void RiemannSolution::mirror(){
    PrimativeState temp = wL;
    wL = wR;
    wR = temp;
    
    temp = sL;
    sL = sR;
    sR = temp;
    
    wL.vx *= -1;
    wR.vx *= -1;
    sL.vx *= -1;
    sR.vx *= -1;
}
