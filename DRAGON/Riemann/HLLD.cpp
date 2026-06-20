//
//  Riemann/HLLD.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 15/06/2026.
//  Implementation based on https://robertcaddy.com/posts/HLLD-Algorithm/
//      Miyoshi and Kusano (2005). https://doi.org/10.1016/j.jcp.2005.02.017




#include "Riemann.hpp"
#include "Constants.h"
#include <cassert>

#ifdef MHD

//MARK: Setup
void compute_outer_star_transverse(PrimitiveState& wsK, const PrimitiveState& K, double _xK, double SK, double SM, double Bx){
    double Bx2_4pi = (Bx*Bx)/(4*M_PI); // B_x^2 / 4pi
    double denom = _xK*(SK - SM) - Bx2_4pi; //Denominator of B and v scaling factors
    //Set transverse components via vector arithmatic, then override x component
    //Magnetic Field
    wsK.B = (fabs(denom)<1e-12) ? K.B : K.B * (_xK*(SK-K.v.x)-Bx2_4pi)/denom;
    wsK.B.x = Bx;
    //Velocities
    wsK.v = K.v - (fabs(denom)<1e-12 ? vec3{0,0,0} :  K.B * Bx*(SM-K.v.x)/(4*M_PI*denom) );
    wsK.v.x = SM;
}




ConservativeState Riemann::HLLD(){
    //Set normal magnetic fields
    double Bx = (L.B.x+R.B.x)/2;
    L.B.x = Bx; R.B.x = Bx;
    //Calculate Outer Wave speeds
    double SL = fmin(L.v.x - L.c_fast(), R.v.x - R.c_fast());
    if(SL >= 0) return L.flux(); //Left intial region
    double SR = fmax(L.v.x + L.c_fast(), R.v.x + R.c_fast());
    if(SR <= 0) return R.flux(); //Right intial region
        
    //Calculate the Contact Wave Speed
    double pTL = L.p + (L.B*L.B)/(8*M_PI), pTR = R.p + (R.B*R.B)/(8*M_PI);
    double _xL = L.rho * (SL-L.v.x), _xR = R.rho * (SR-R.v.x);
    double SM  = ((_xR*R.v.x - _xL*L.v.x) - (pTR-pTL))  / (_xR - _xL);
    
    //Calculate the Alvfen Wave Speeds
    PrimitiveState wsL, wsR;
    wsL.rho = _xL / (SL - SM);
    double SsL = SM - fabs(Bx) / sqrt(4*M_PI * wsL.rho);
    wsR.rho = _xR / (SR - SM);
    double SsR = SM + fabs(Bx) / sqrt(4*M_PI * wsR.rho);

    
    //Calculate the regions between fast/alfven
    ConservativeState FsL, FsR;
    double pT = (_xR*pTL - _xL*pTR + _xL*_xR*(R.v.x - L.v.x)) / (_xR - _xL);
    
    if(SsR >= 0) {//If this isn't true, we are between right fast/alfven, so we won't need FsL
        compute_outer_star_transverse(wsL, L, _xL, SL, SM, Bx);
        wsL.p = pT - wsL.B*wsL.B / (8*M_PI);
        FsL =  L.flux() + (wsL - L)*SL; FsL.B.x = 0;
       
    }
    if(SsL <= 0) { //If this isn't true, we are between left fast/alfven, so we won't need FsR
        compute_outer_star_transverse(wsR, R, _xR, SR, SM, Bx);
        wsR.p = pT - wsR.B*wsR.B / (8*M_PI);
        FsR = R.flux() + (wsR - R)*SR; FsR.B.x = 0;
    }
            
    if(SsL >= 0) return FsL;//Check for Left fast-alfven region, if so can return the flux now
    if(SsR <= 0) return FsR;//Check for Right fast-alfven region, if so can return the flux now

    
    //Calculate the regions between Alfven/Contact
    double sqL = sqrt(wsL.rho), sqR = sqrt(wsR.rho);
    double sbx = Bx > 0 ? sqrt(4*M_PI) : (Bx < 0 ? -sqrt(4*M_PI) : 0);
    PrimitiveState wss;
    //Set transverse components via vector arithmatic, then override x component
    wss.v = (sqL*wsL.v + sqR*wsR.v + sbx*(wsR.B - wsL.B)/(4*M_PI))/(sqL+sqR);
    wss.v.x = SM;
    wss.B = (sqL*wsR.B + sqR*wsL.B + sqL*sqR*sbx*(wsR.v - wsL.v))/(sqL+sqR);
    wss.B.x = Bx;
    wss.p = pT - wss.B*wss.B / (8*M_PI);
    
    //Copy density and calculate final flux
    if(SM >= 0){//Left
        wss.rho = wsL.rho;
        return FsL + (wss-wsL)*SsL;
    } else {//Right
        wss.rho = wsR.rho;
        return FsR + (wss-wsR)*SsR;
    }
}
#elif RIEMANN_DEFAULT == RIEMANN_HLLD
#error HLLD Solver requires MHD
#endif
