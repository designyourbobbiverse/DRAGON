//
//  Exact.cpp
//  DRAGON/Hydro/Riemann
//
//  Created by Bobbie Markwick on 28/04/2026.
//  Implementation based mostly on Toro (2009). https://doi.org/10.1007/b79761
//

#include "Hydro/Riemann/Riemann.hpp"

#include "Constants.h" //For _gamma and related coefficients
#include <cmath> //For std::sqrt, pow, abs
#include <algorithm> //For std::min
#include <utility> //For std::swap
using namespace DRAGON;


//MARK: Setup
Riemann::Riemann(PrimitiveState _L, PrimitiveState _R){ L = _L; R = _R; }
RiemannSolution::RiemannSolution(Riemann problem){
    wL = problem.L;
    wR = problem.R;
    //Copy Transverse Values
    sL.v = wL.v;
    sR.v = wR.v;
}


//MARK: Velocity Jump Function

// Guide to _gamma coefficients
// Rules: _X_Y = X/Y, G = _gamma, m = -, p = +
// Example: _Gp1_2G = (_gamma+1)/(2*_gamma)
// Example: _2_Gm1 = 2/(_gamma - 1)

double Riemann::f(double p, PrimitiveState w){
    if (p > w.p) { //Shock
        double A = _2_Gp1/w.rho, B = _Gm1_Gp1 * w.p;
        return (p - w.p) * std::sqrt(A/(B+p));
    } else {//Rarefaction
        return _2_Gm1 * w.cs() * (std::pow(p/w.p, _Gm1_2G)-1);
    }
}
double df(double p, PrimitiveState w){
    if (p > w.p) { //Shock
        double A = _2_Gp1/w.rho, B = _Gm1_Gp1 * w.p;
        return std::sqrt(A/(p+B)) * (1-(p-w.p)/(2*(B+p)));
    } else {//Rarefaction
        return  std::pow(p/w.p, -_Gp1_2G) / (w.rho * w.cs());
    }
}

//MARK: Exact Riemann Solver
RiemannSolution Riemann::exact(){
    return exact( (L.p + R.p)/2.0 ); //Default Initial Guess = average of sides
}
RiemannSolution Riemann::exact(double pGuess){
    #ifdef Exact_Rarefactions_Check
    //Check to see if the wave will be two rarefactions. If so, a closed solution exists.
    double p_min = std::min(L.p, R.p);
    if(f(p_min,L) + f(p_min,R) + R.v.x - L.v.x >= 0) return TRRS();
    #endif
    
    RiemannSolution s = RiemannSolution(*this);
    //Pressure
    s.sL.p = exact_StarP(pGuess);
        s.sR.p = s.sL.p; //Pressure same across contact
    //velocity
    s.sL.v.x = exact_StarV(s.sL.p);
        s.sR.v.x = s.sL.v.x; //Velocity same across contact
    //Density
    s.sL.rho = exact_StarRho(L, s.sL.p);
    s.sR.rho = exact_StarRho(R, s.sR.p);
    
    return s;
}



double Riemann::exact_StarP(double pGuess){
    double pStar = pGuess, CHA = 1;
    int iters = CONFIG::ExactRiemann_MaxIters; //Keeps track of hwo many iterations we have left
    do{ //Newton's Method on the Velocity jump equation
        double fp = f(pStar, L) + f(pStar, R) + R.v.x - L.v.x;
        double dfdp = df(pStar, L) + df(pStar, R);
        double delta = std::min(fp/dfdp, 0.8*pStar); //0.8 helps keep from going negative/diverging
        pStar -= delta;
        CHA = std::abs(delta/(pStar+delta/2));
    } while(CHA > CONFIG::ExactRiemann_Tolerance &&  --iters != 0 );
    return pStar;
}
//v* and rho* given p*
double Riemann::exact_StarV(double pStar){
    return (L.v.x + R.v.x + f(pStar,R)-f(pStar,L))/2.0;
}
double Riemann::exact_StarRho(PrimitiveState w, double p){
    if(p>w.p) return w.rho * (p+_Gm1_Gp1*w.p)/(_Gm1_Gp1*p+w.p);//Shock
    else return w.rho * std::pow( p/w.p, _Ginv); //Rarefaction
}


//MARK: Two-Rarefaction
RiemannSolution Riemann::TRRS(){
    double aL = L.cs(), aR = R.cs();//Compute the Sound Speeds
    RiemannSolution s = RiemannSolution(*this);

    double _LR = std::pow(L.p / R.p, _Gm1_2G) * aR/aL;
    
    //Velocity
    s.sL.v.x = (_LR*L.v.x + R.v.x + _2_Gm1*(_LR*aL-aR)) / (_LR + 1);
        s.sR.v.x = s.sL.v.x;//Velocity same across contact
    //Pressure
    s.sL.p = L.p * std::pow(1 - _Gm1_2*(s.sL.v.x - L.v.x)/aL, _2G_Gm1);
        s.sR.p = s.sL.p;//Pressure same across contact
    //Density
    s.sL.rho = L.rho * std::pow(s.sL.p / L.p , _Ginv);
    s.sR.rho = R.rho * std::pow(s.sR.p / R.p , _Ginv);
    
    return s;
}



//MARK: Solution Sampling
//Convenience methods to sample the solution then get the flux
ConservativeState RiemannSolution::flux(){ return flux(0); }
ConservativeState RiemannSolution::flux(double x_t){
    PrimitiveState w = sample(x_t);
    return ConservativeState(w).flux(w.v);
}

//Get the state along any given x/t line
PrimitiveState RiemannSolution::sample(double x_t){
    PrimitiveState state;
    //Edge case: Right on the contact wave, do this to help ensure symmetry
    if(std::abs(sR.v.x - x_t) < 1e-12){
        double sql = std::sqrt(sL.rho), sqr = std::sqrt(sR.rho);
        return (sql*sL + sqr*sR)/(sql+sqr);
    }
    
    //Handle left vs Right side
    bool isLeft = x_t < sR.v.x;
    if(isLeft){ mirror(); x_t=-x_t; }
    
    //Calculate Hydro Sound Speed
    double a = wR.cs();
    //Determine Zone
    int zone = 0; // 1 = outside, 2 = fan, 3 = star
    if (sR.p > wR.p){ //shock
        double scale = std::sqrt(_Gp1_2G*sR.p/wR.p + _Gm1_2G );
        zone = ( (x_t-wR.v.x)  > scale * a ) ? 1 : 3; //Shocks have no fan, just jump form outside to star
    } else{ //Rarefraction
        if ( (x_t-wR.v.x) > a ) zone = 1; //Outside region
        else if ( (x_t-sR.v.x) > a ) zone = 2; //one extra compare that saves a call to pow
        else if ( (x_t-sR.v.x) > a * std::pow(sR.p/wR.p, _Gm1_2G) ) zone = 2; //Fan
        else zone = 3; //Inside the star region
    }
    //Calculate the State at x/t
    switch(zone){
        case 1: state=wR; break; //Outer Region
        case 3: state=sR;break; //Star Region
        case 2://Fan
            state = wR;
            double scale = _2_Gp1 - _Gm1_Gp1 * (wR.v.x-x_t)/a;
            state.rho *= std::pow(scale, _2_Gm1);
            state.p *= std::pow(scale,_2G_Gm1);
            state.v.x = _2_Gp1 * (x_t-a + wR.v.x * _Gm1_2);
            break;
    }
    //Be a good citizen, restore original state if we mirrored
    if(isLeft){ mirror(); state.v.x *= -1;}

    return state;
}

void RiemannSolution::mirror(){
    //Swap corresponding Left and Right States
    std::swap(wL,wR);
    std::swap(sL,sR);
    //Multiply all x components by -1
    wL.v.x *= -1;
    wR.v.x *= -1;
    sL.v.x *= -1;
    sR.v.x *= -1;
    #ifdef MHD
    wL.B.x *= -1;
    wR.B.x *= -1;
    sL.B.x *= -1;
    sR.B.x *= -1;
    #endif
}
