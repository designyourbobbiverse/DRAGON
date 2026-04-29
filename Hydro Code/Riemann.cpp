//
//  Riemann.cpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 28/04/2026.
//

#include "Riemann.hpp"
#include <math.h>
#include "Constants.h"


//MARK: Gamma-dependent constants
const double  _Gm1_2G = (gamma - 1.0)/(2.0 * gamma);
const double _Gp1_2G = (gamma + 1.0)/(2.0*gamma);
const double _2G_Gm1 = 2.0*gamma/(gamma - 1.0);
const double _Gm1_2 = (gamma - 1.0)/2.0;
const double _2_Gm1 = 2.0/(gamma - 1.0);
const double _2_Gp1 = 2.0/(gamma + 1.0);
const double _Gm1_Gp1 = (gamma - 1.0)/(gamma + 1.0);
const double _Ginv = 1.0/gamma;

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

//Velocity Jump Function
double f(double p, PrimativeState w){
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


//Components: Exact Iterative Soler
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


//MARK: HLL Solvers
ConservativeState Riemann::HLL(double sl, double sr){
    //Outside region
    if(sl > 0) return ConservativeState(L).flux(L.vx);
    if(sr < 0) return ConservativeState(R).flux(R.vx);
    
    ConservativeState UL = ConservativeState(L);
    ConservativeState UR = ConservativeState(R);
    ConservativeState FL = UL.flux(L.vx);
    ConservativeState FR = UR.flux(R.vx);

    return (sr*FL - sl*FR + sl*sr*(UR-UL))/(sr-sl);
}


ConservativeState Riemann::HLLC(double sl, double sr){
    //Outside region
    if(sl >= 0) return ConservativeState(L).flux(L.vx);
    if(sr <= 0) return ConservativeState(R).flux(R.vx);
    //Calculate the contact wave
    double _pr = R.p + R.rho*R.vx*(R.vx - sr), _pl = L.p + L.rho*L.vx*(L.vx - sl);
    double sc = (_pr - _pl) / (R.rho*(R.vx-sr) - L.rho*(L.vx-sl));
    //Left or Right
    PrimativeState X = (sc > 0 ? L : R);
    double sx = sc > 0 ? sl : sr;
    ConservativeState UX = ConservativeState(X);
    //Compute Star Region
    ConservativeState U = UX * (sx - X.vx) / (sx - sc);
    U.px = U.rho * sc;
    U.E +=  U.rho*(sc - X.vx)*(sc + X.p/(X.rho*(sx-X.vx)) );
    //Compute Flux
    return UX.flux(X.vx) + (U - UX) * sx;
}


ConservativeState Riemann::HLLC(){
    //Estimate the pressure
    double aL = sqrt(gamma * L.p/L.rho), aR = sqrt(gamma * R.p/L.rho); // Sound Speeds
    double p_pvrs = (L.p + R.p)/2 + (L.rho + R.rho)*(L.vx - R.vx)*(aL + aR)/8;
    //Left Speed
    double SL = aL;
    if(p_pvrs > L.p) SL *= sqrt(1 + _Gp1_2G*(p_pvrs/L.p - 1));
    SL = L.vx - SL;
    //Right Speed
    double SR = aR;
    if(p_pvrs > R.p) SL *= sqrt(1 + _Gp1_2G*(p_pvrs/R.p - 1));
    SR = R.vx + SR;
    
    return HLLC(SL, SR);
}

ConservativeState Riemann::HLL(){
    //Estimate the pressure
    double aL = sqrt(gamma * L.p/L.rho), aR = sqrt(gamma * R.p/L.rho); // Sound Speeds
    double p_pvrs = (L.p + R.p)/2 + (L.rho + R.rho)*(L.vx - R.vx)*(aL + aR)/8;
    //Left Speed
    double SL = aL;
    if(p_pvrs > L.p) SL *= sqrt(1 + _Gp1_2G*(p_pvrs/L.p - 1));
    SL = L.vx - SL;
    //Right Speed
    double SR = aR;
    if(p_pvrs > R.p) SL *= sqrt(1 + _Gp1_2G*(p_pvrs/R.p - 1));
    SR = R.vx + SR;
    
    return HLL(SL, SR);
}


//MARK: Roe Solver
//Roe Eigenvectors
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
ConservativeState Riemann::Roe(){
    double sqL = sqrt(L.rho), sqR = sqrt(R.rho);
    //weighted averages
    double _rho = sqL * sqR;
    double _vx = (sqL * L.vx + sqR * R.vx)/(sqL+sqR);
    double _vy = (sqL * L.vy + sqR * R.vy)/(sqL+sqR);
    double _vz = (sqL * L.vz + sqR * R.vz)/(sqL+sqR);
    double _H = (sqL*L.enthalpy() + sqR*R.enthalpy()) / (sqL + sqR);
    double _V2 = _vx*_vx + _vy*_vy + _vz*_vz;
    double _a = sqrt((gamma-1) * (_H - _V2/2));
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
    double aL = sqrt(gamma * L.p/L.rho), aSL = sqrt(gamma*SL.pressure()/SL.rho);
    double lambdaL = L.vx-aL, lambdaLS = SL.px/SL.rho - aSL;
    if(lambdaL < 0 &&  lambdaLS > 0 ) {//Left Rarefaction
        double _lambda = lambdaL * (lambda[0] - lambdaLS)/(lambdaL - lambdaLS);
        return UL.flux(L.vx) + _lambda*alpha[0]*K[0];
    }
    ConservativeState SR = UR - alpha[4]*K[4];
    double aR = sqrt(gamma * R.p/L.rho), aSR = sqrt(gamma*SR.pressure()/SR.rho);
    double lambdaR = R.vx+aR, lambdaRS = SR.px/SR.rho + aSR;
    if(lambdaR > 0 &&  lambdaRS  < 0) {//Right Rarefaction
        double _lambda = lambdaR * (lambda[4] - lambdaRS)/(lambdaR - lambdaRS);
        return UR.flux(R.vx) - _lambda*alpha[4]*K[4];
    }
#endif
    //Combine the waves
    ConservativeState F = (UL.flux(L.vx) + UR.flux(R.vx));
    for(int i=0;i<5;i++) F += alpha[i]*fabs(lambda[i])*K[i];
    return F/2;
}


//MARK: Other Approximate Solvers
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

//Two-Rarefaction Approximation
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

//Two-Shock Approximation
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










