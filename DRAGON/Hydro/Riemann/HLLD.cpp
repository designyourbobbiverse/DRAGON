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
namespace HLL{ PrimitiveState roeAvg(PrimitiveState L, PrimitiveState R); }

void compute_outer_star(ConservativeState& usK, const PrimitiveState& K, double _xK, double SK, double SM, double pT, double Bx){
    double Bx2_4pi = (Bx*Bx) * _1_4pi; // B_x^2 / 4pi
    double denom = _xK*(SK - SM) - Bx2_4pi; //Denominator of B and v scaling factors
    //Set transverse components via vector arithmatic, then override x component
    //Magnetic Field
    usK.B = (fabs(denom)<1e-12) ? K.B : K.B * (_xK*(SK-K.v.x)-Bx2_4pi)/denom;
    usK.B.x = Bx;
    //Velocities
    auto vsK = K.v - (fabs(denom)<1e-12 ? vec3{0,0,0} :  K.B * Bx*(SM-K.v.x)/denom * _1_4pi );
    vsK.x = SM;
    usK.p = usK.rho * vsK;
    //Energy
    usK.E = K.energy() * (SK - K.v.x);
    usK.E += pT*SM - (K.p + K.B*K.B*_1_8pi)*K.v.x;
    usK.E += Bx * (K.v*K.B - vsK*usK.B) * _1_4pi;
    usK.E /= (SK - SM);
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
    double pTL = L.p + (L.B*L.B)*_1_8pi, pTR = R.p + (R.B*R.B)*_1_8pi;
    double _xL = L.rho * (SL-L.v.x), _xR = R.rho * (SR-R.v.x);
    double SM  = ((_xR*R.v.x - _xL*L.v.x) - (pTR-pTL))  / (_xR - _xL);
    
    //Calculate the Alvfen Wave Speeds
    ConservativeState usL, usR;
    usL.rho = _xL / (SL - SM);
    double SsL = SM - fabs(Bx) / (sq4pi * sqrt(usL.rho));
    usR.rho = _xR / (SR - SM);
    double SsR = SM + fabs(Bx) / (sq4pi * sqrt(usR.rho));
   
    
    //Calculate the regions between fast/alfven
    ConservativeState FsL, FsR;
    PrimitiveState wsL, wsR;
    double pT = (_xR*pTL - _xL*pTR + _xL*_xR*(R.v.x - L.v.x)) / (_xR - _xL);
    
    if(SsR >= 0) {//If this isn't true, we are between right fast/alfven, so we won't need FsL
        compute_outer_star(usL, L, _xL, SL, SM, pT, Bx);
        wsL.rho = usL.rho;
        wsL.v = usL.p / usL.rho;
        wsL.B = usL.B;

        #ifdef HLLD_PHYSICAL_SAFETY
        wsL.p = pT - wsL.B*wsL.B*_1_8pi;
        if(!wsL.isPhysical() || !usL.isPhysical()) return HLLE();
        #endif
        
        FsL =  L.flux() + (usL - L)*SL; FsL.B.x = 0;
    }
    if(SsL <= 0) { //If this isn't true, we are between left fast/alfven, so we won't need FsR
        compute_outer_star(usR, R, _xR, SR, SM, pT, Bx);
        wsR.rho = usR.rho;
        wsR.v = usR.p / usR.rho;
        wsR.B = usR.B;

        #ifdef HLLD_PHYSICAL_SAFETY
        wsR.p = pT - wsR.B*wsR.B*_1_8pi;
        if(!wsR.isPhysical() || !usR.isPhysical()) return HLLE();
        #endif
        
        FsR = R.flux() + (usR - R)*SR; FsR.B.x = 0;
    }
            
    if(SsL >= 0) return FsL;//Check for Left fast-alfven region, if so can return the flux now
    if(SsR <= 0) return FsR;//Check for Right fast-alfven region, if so can return the flux now

    auto ws = SM >= 0 ? wsL : wsR;

    //Calculate the regions between Alfven/Contact
    double sqL = sqrt(wsL.rho), sqR = sqrt(wsR.rho);
    double _sqL = sqL / (sqL + sqR), _sqR = sqR / (sqL + sqR);
    double sbx = Bx > 0 ? sq4pi : (Bx < 0 ? -sq4pi : 0);
    ConservativeState uss;
    PrimitiveState wss;
    wss.rho = ws.rho;
    uss.rho = wss.rho;
    //Set transverse components via vector arithmatic, then override x component
    wss.v = _sqL*wsL.v + _sqR*wsR.v + (sbx*(usR.B - usL.B)*_1_4pi)/(sqL+sqR);
    wss.v.x = SM;
    uss.p = uss.rho * wss.v;
    wss.B = _sqL*usR.B + _sqR*usL.B + sbx*(wsR.v - wsL.v)*sqL*sqR/(sqL+sqR);
    wss.B.x = Bx;
    uss.B = wss.B;
    //Energy
    auto dE = sbx * (ws.v*ws.B - wss.v*wss.B) * _1_4pi * sqrt(wss.rho);
    uss.E = SM >= 0 ? usL.E - dE : usR.E + dE;
    
    #ifdef HLLD_PHYSICAL_SAFETY
    wss.p = pT - wss.B*wss.B*_1_8pi;
    if(!wss.isPhysical() || !uss.isFinite()) return HLLE();
    #endif
    
    //Copy density and calculate final flux
    if(SM >= 0){//Left
        return FsL + (uss-usL)*SsL;
    } else {//Right
        return FsR + (uss-usR)*SsR;
    }
}


//MARK: Zero normal
void compute_outer_star_zero_B(ConservativeState& usK, const PrimitiveState& K, double _xK, double SK, double SM, double pT){
    //Set transverse components via vector arithmatic, then override x component
    //Magnetic Field
    usK.B = (fabs(SK - SM)<1e-12) ? K.B : K.B * (SK-K.v.x)/(SK - SM);
    usK.B.x = 0;
    //Velocities
    auto vsK = K.v;
    vsK.x = SM;
    usK.p = usK.rho * vsK;
    //Energy
    usK.E = K.energy() * (SK - K.v.x);
    usK.E += pT*SM - (K.p + K.B*K.B*_1_8pi)*K.v.x;
    usK.E /= (SK - SM);
}

ConservativeState Riemann::HLLD_zero_B(double SL, double SR){
    
    //Calculate the Contact Wave Speed
    double pTL = L.p + (L.B*L.B)*_1_8pi, pTR = R.p + (R.B*R.B)*_1_8pi;
    double _xL = L.rho * (SL-L.v.x), _xR = R.rho * (SR-R.v.x);
    double SM  = ((_xR*R.v.x - _xL*L.v.x) - (pTR-pTL))  / (_xR - _xL);
    
    //Calculate the Alvfen Wave Speeds
    ConservativeState usL, usR;
    usL.rho = _xL / (SL - SM);
    usR.rho = _xR / (SR - SM);
    
    //Calculate the regions between fast/alfven
    ConservativeState FsL, FsR;
    PrimitiveState wsL, wsR;
    double pT = (_xR*pTL - _xL*pTR + _xL*_xR*(R.v.x - L.v.x)) / (_xR - _xL);
    
    if(SM >= 0) {//If this isn't true, we are between right fast/alfven, so we won't need FsL
        compute_outer_star_zero_B(usL, L, _xL, SL, SM, pT);
        wsL.rho = usL.rho;
        wsL.v = usL.p / usL.rho;
        wsL.B = usL.B;

        #ifdef HLLD_PHYSICAL_SAFETY
        wsL.p = pT - wsL.B*wsL.B*_1_8pi;
        if(!wsL.isPhysical() || !usL.isFinite()) return HLLE();
        #endif
        
        FsL =  L.flux() + (usL - L)*SL; FsL.B.x = 0;
    }
    if(SM <= 0) { //If this isn't true, we are between left fast/alfven, so we won't need FsR
        compute_outer_star_zero_B(usR, R, _xR, SR, SM, pT);
        wsR.rho = usR.rho;
        wsR.v = usR.p / usR.rho;
        wsR.B = usR.B;

        #ifdef HLLD_PHYSICAL_SAFETY
        wsR.p = pT - wsR.B*wsR.B*_1_8pi;
        if(!wsR.isPhysical() || !usR.isFinite()) return HLLE();
        #endif
        
        FsR = R.flux() + (usR - R)*SR; FsR.B.x = 0;
    }
            
    if(SM > 0) return FsL;//Check for Left fast-alfven region, if so can return the flux now
    if(SM < 0) return FsR;//Check for Right fast-alfven region, if so can return the flux now

    
    //Special case: right on the contact wave (down to the bit)
    //Average the two outputs to help preserve symmetry
    double sql = sqrt(wsL.rho), sqr = sqrt(wsR.rho);
    double _sql = sql / (sql + sqr), _sqr = sqr / (sql + sqr);
    return _sql * FsL + _sqr * FsR;
}

#elif RIEMANN_DEFAULT == RIEMANN_HLLD
#error HLLD Solver requires MHD
#endif
