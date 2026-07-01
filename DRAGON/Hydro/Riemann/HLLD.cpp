//
//  HLLD.cpp
//  DRAGON/Hydro/Riemann
//
//  Created by Bobbie Markwick on 15/06/2026.
//  Implementation based on https://robertcaddy.com/posts/HLLD-Algorithm/
//      Miyoshi and Kusano (2005). https://doi.org/10.1016/j.jcp.2005.02.017




#include "Riemann.hpp"
#include "Constants.h"
#include <cassert>

#ifdef MHD

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
    double cl = L.c_fast(), cr = R.c_fast();
    double SL = fmin(L.v.x - cl, R.v.x - cr);
    if(SL >= 0) return L.flux(); //Left intial region
    double SR = fmax(L.v.x + cl, R.v.x + cr);
    if(SR <= 0) return R.flux(); //Right intial region
    
    //Check for zero-normal case
    auto cutoff = fmax(L.B*L.B+ R.B*R.B,1e-100);
    if(Bx*Bx < 1e-12 * cutoff) return HLLD_zero_B(SL, SR);
            
    //Calculate the Contact Wave Speed
    double pTL = L.p + (L.B*L.B)/(8*M_PI), pTR = R.p + (R.B*R.B)/(8*M_PI);
    double _xL = L.rho * (SL-L.v.x), _xR = R.rho * (SR-R.v.x);
    double SM  = ((_xR*R.v.x - _xL*L.v.x) - (pTR-pTL))  / (_xR - _xL);
    
    //Calculate the Alvfen Wave Speeds
    PrimitiveState wsL, wsR;
    //ConservativeState usL, usR;
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
        
        #ifdef HLLD_PHYSICAL_SAFETY
        if(!wsL.isPhysical()) return HLLE();
        #endif
        
        FsL =  L.flux() + (wsL - L)*SL; FsL.B.x = 0;
    }
    if(SsL <= 0) { //If this isn't true, we are between left fast/alfven, so we won't need FsR
        compute_outer_star_transverse(wsR, R, _xR, SR, SM, Bx);
        wsR.p = pT - wsR.B*wsR.B / (8*M_PI);

        #ifdef HLLD_PHYSICAL_SAFETY
        if(!wsR.isPhysical()) return HLLE();
        #endif
        
        FsR = R.flux() + (wsR - R)*SR; FsR.B.x = 0;
    }
            
    if(SsL >= 0) return FsL;//Check for Left fast-alfven region, if so can return the flux now
    if(SsR <= 0) return FsR;//Check for Right fast-alfven region, if so can return the flux now

    auto ws = SM >= 0 ? wsL : wsR;

    //Calculate the regions between Alfven/Contact
    double sqL = sqrt(wsL.rho), sqR = sqrt(wsR.rho);
    double sbx = Bx > 0 ? sqrt(4*M_PI) : (Bx < 0 ? -sqrt(4*M_PI) : 0);
    PrimitiveState wss;
    // ConservativeState uss;
    //Set transverse components via vector arithmatic, then override x component
    wss.rho = ws.rho;
    wss.v = (sqL*wsL.v + sqR*wsR.v + sbx*(wsR.B - wsL.B)/(4*M_PI))/(sqL+sqR);
    wss.v.x = SM;
    wss.B = (sqL*wsR.B + sqR*wsL.B + sqL*sqR*sbx*(wsR.v - wsL.v))/(sqL+sqR);
    wss.B.x = Bx;
    wss.p = pT - wss.B*wss.B / (8*M_PI);
    
    #ifdef HLLD_PHYSICAL_SAFETY
    if(!wss.isPhysical()) return HLLE();
    #endif
    
    //Copy density and calculate final flux
    if(SM >= 0){//Left
        return FsL + (wss-wsL)*SsL;
    } else {//Right
        return FsR + (wss-wsR)*SsR;
    }
}


//MARK: Zero normal
void compute_outer_star_transverse_zero_B(PrimitiveState& wsK, const PrimitiveState& K, double _xK, double SK, double SM){
    //Set transverse components via vector arithmatic, then override x component
    //Magnetic Field
    wsK.B = (fabs(SK - SM)<1e-12) ? K.B : K.B * (SK-K.v.x)/(SK - SM);
    wsK.B.x = 0;
    //Velocities
    wsK.v = K.v;
    wsK.v.x = SM;
}

ConservativeState Riemann::HLLD_zero_B(double SL, double SR){
    
    //Calculate the Contact Wave Speed
    double pTL = L.p + (L.B*L.B)/(8*M_PI), pTR = R.p + (R.B*R.B)/(8*M_PI);
    double _xL = L.rho * (SL-L.v.x), _xR = R.rho * (SR-R.v.x);
    double SM  = ((_xR*R.v.x - _xL*L.v.x) - (pTR-pTL))  / (_xR - _xL);
    
    //Calculate the Alvfen Wave Speeds
    PrimitiveState wsL, wsR;
    wsL.rho = _xL / (SL - SM);
    wsR.rho = _xR / (SR - SM);
    
    //Calculate the regions between fast/alfven
    ConservativeState FsL, FsR;
    double pT = (_xR*pTL - _xL*pTR + _xL*_xR*(R.v.x - L.v.x)) / (_xR - _xL);
    
    if(SM >= 0) {//If this isn't true, we are between right fast/alfven, so we won't need FsL
        compute_outer_star_transverse_zero_B(wsL, L, _xL, SL, SM);
        wsL.p = pT - wsL.B*wsL.B / (8*M_PI);
        
        #ifdef HLLD_PHYSICAL_SAFETY
        if(!wsL.isPhysical()) return HLLE();
        #endif
        
        FsL =  L.flux() + (wsL - L)*SL; FsL.B.x = 0;
    }
    if(SM <= 0) { //If this isn't true, we are between left fast/alfven, so we won't need FsR
        compute_outer_star_transverse_zero_B(wsR, R, _xR, SR, SM);
        wsR.p = pT - wsR.B*wsR.B / (8*M_PI);

        #ifdef HLLD_PHYSICAL_SAFETY
        if(!wsR.isPhysical()) return HLLE();
        #endif
        
        FsR = R.flux() + (wsR - R)*SR; FsR.B.x = 0;
    }
            
    if(SM > 0) return FsL;//Check for Left fast-alfven region, if so can return the flux now
    if(SM < 0) return FsR;//Check for Right fast-alfven region, if so can return the flux now

    
    //Special case: right on the contact wave
    double sqL = sqrt(wsL.rho), sqR = sqrt(wsR.rho);
    PrimitiveState wss = (sqL*wsL + sqR*wsR)/(sqL+sqR);
    wss.rho = wsL.rho;
    wss.p = pT - wss.B*wss.B / (8*M_PI);

    #ifdef HLLD_PHYSICAL_SAFETY
    if(!wss.isPhysical()) return HLLE();
    #endif
    
    return FsL + (wss-wsL)*SM;
}

#elif RIEMANN_DEFAULT == RIEMANN_HLLD
#error HLLD Solver requires MHD
#endif
