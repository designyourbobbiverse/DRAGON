//
//  Godunov.cpp
//  Hydro Code
//
//  Created by Bobbie Markwick on 12/05/2026.
//

#include "Godunov.hpp"
#include "Constants.h"
#include "Riemann.hpp"


void advance1D(PrimitiveState* w, int size, double time){
    ConservativeState flux = Riemann(w[0],w[0]).flux();
    
    ConservativeState fL = flux, fR;
    for(int i=0; i<size; i++){
        PrimitiveState wR = i==size-1 ? w[i] : w[i+1];
        ConservativeState u = ConservativeState(w[i]);
      /*  fR = tvdFlux(w, i, size, time);
       */
        u += (fL - fR) * (time/stepSize);
        w[i] = PrimitiveState(u);
        fL = fR;
    }
}


/*
//WAF
double rk(double ck, state*w, int i, int size){
    double r = 0;
    if (ck > 0 && i > 0)  r = (w[i].rho - w[i-1].rho)/(w[i+1].rho - w[i].rho);
    else if (ck < 0 && i < size-2)  r = (w[i+2].rho - w[i+1].rho)/(w[i+1].rho - w[i].rho);
    if ( r != r ) return -1;
    
    double phi ;
    if (r <= 0) phi = 1;
    else if (r > 1) phi = fabs(ck);
    else  phi =  1-(1-fabs(ck))*r;

    double sign = fabs(ck) < 0.001 ? 0 : ck/fabs(ck);
 //   if (sign != sign) abort();
    return phi * sign / 2;
}


stateConsv tvdFlux(state* w, int i, int size, double time){
    state starL, starR;
    state wL = w[i], wR = i==size-1 ? w[i] : w[i+1];
    ExactRiemann(wL, wR, &starL, &starR);
    //Handle Rarefractions
    starL = riemannSample( wL, starL, 0);
    starR = riemannSample( wR, starR, 0);
    
    double  u =  starL.u, a = sqrt(gamma * fabs(u < 0 ? starR.p/starR.rho : starL.p/starL.rho));
    double c1 = (u-a) * time / stepSize;
    double c2 = u * time / stepSize;
    double c3 = (u+a) * time / stepSize;
    double r1 = rk(c1, w,i,size), r2 = rk(c2, w,i,size), r3 = rk(c3, w,i,size);
    
    
    stateConsv F1 = flux(conservativeForm(wL), wL.u);
    stateConsv F2 = flux(conservativeForm(starL), starL.u);
    stateConsv F3 = flux(conservativeForm(starR), starR.u);
    stateConsv F4 = flux(conservativeForm(wR), wR.u);
    
    stateConsv waf;
    waf.rho =  (r1+0.5)*F1.rho + (r2-r1)*F2.rho + (r3-r2)*F3.rho + (0.5-r3)*F4.rho;
    waf.px =  (r1+0.5)*F1.px + (r2-r1)*F2.px + (r3-r2)*F3.px + (0.5-r3)*F4.px;
    waf.py =  (r1+0.5)*F1.py + (r2-r1)*F2.py + (r3-r2)*F3.py + (0.5-r3)*F4.py;
    waf.E =  (r1+0.5)*F1.E + (r2-r1)*F2.E + (r3-r2)*F3.E + (0.5-r3)*F4.E;
    return waf;
}
*/
