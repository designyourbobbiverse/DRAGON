//  fluid solver
//  Created by Bobbie Markwick on 21/07/21.
//

#include <iostream>
#include <math.h>

const double gamma = 5.0/3.0;
const double  G1 = (gamma - 1.0)/(2.0 * gamma);
const double G2 = (gamma + 1.0)/(2.0*gamma);
const double G3 = 2.0*gamma/(gamma - 1.0);
const double G4 = 2.0/(gamma - 1.0);
const double G5 = 2.0/(gamma + 1.0);
const double G6 = (gamma - 1.0)/(gamma + 1.0);
const double G7 = (gamma - 1.0)/2.0;
const double G8 = 1.0/gamma;
const double G9 = gamma / (gamma - 1.0);

const double TOL = 1E-13;
const double stepSize = 0.01;
const double timeStep = 0.01;


            //State Components
            struct state { double rho,  u, v,  p; };

            struct stateConsv { double rho,  px, py, E; };


            stateConsv conservativeForm(state w){
                stateConsv u;
                u.rho = w.rho;
                u.px = w.rho * w.u;
                u.py = w.rho * w.v;
                u.E = w.p/(gamma-1) + w.rho * (pow(w.u,2) + pow(w.v,2))/2;
                return u;
            }
            double pressure(stateConsv u){
                return (gamma-1.0) * (u.E - (pow(u.px,2)+ pow(u.py,2))/(2 * u.rho));
            }
            state primativeForm(stateConsv u){
                state w;
                w.rho = u.rho;
                w.u = u.px/u.rho;
                w.v = u.py/u.rho;
                w.p = pressure(u);
                return w;
            }
            double H(state w){
                return G9 * (w.p/w.rho) + (pow(w.u,2)+pow(w.v,2))/2;
            }
            stateConsv flux(stateConsv u, double v){
                stateConsv F;
                double p = pressure(u);
                F.rho = u.px;
                F.px = u.px * v + p;
                F.py = u.py * v;
                F.E = v * (u.E + p);
                return F;
            }
            double S(state w){
                return sqrt(w.u*w.u + w.v*w.v) + sqrt(gamma*w.p/w.rho);
            }
            state transform(state w){
                state wT;
                wT.rho = w.rho; wT.p = w.p;
                wT.u = w.v; wT.v = w.u;
                return wT;
            }


            //Exact Riemann Solver Components
            double f(double p, state w){
                if (p > w.p) { //Shock
                    double A= G5/w.rho, B= G6 * w.p;
                    return (p-w.p) * sqrt(A/(B+p));
                } else {//Rarefaction
                    double a = sqrt(gamma * w.p / w.rho);
                    return G4 * a * (pow(p/w.p, G1)-1);
                }
            }
            double df(double p, state w){
                if (p > w.p) { //Shock
                    double A= G5/w.rho, B= G6 * w.p;
                    return sqrt(A/(p+B)) * (1-(p-w.p)/(2*(B+p)));
                }
                else {//Rarefaction
                    double a = sqrt(gamma * w.p / w.rho);
                    return  pow(p/w.p, -G2) / (w.rho * a);
                }
            }
            double solvePStar(state wL, state wR, double pGuess){
                double pStar = pGuess, CHA = 1;
                do{
                    double fp = f(pStar, wL) + f(pStar, wR) + wR.u-wL.u;
                    double dfdp = df(pStar, wL) + df(pStar, wR);
                    double delta = fmin(fp/dfdp, 0.8*pStar);
                    pStar -= delta;
                    CHA = fabs(delta/(pStar+delta/2));
                } while(CHA > TOL);
                return pStar;
            }
            double solvePStar(state wL, state wR){
                return solvePStar(wL, wR, (wL.p+wR.p)/2.0 );
            }
            double solveUStar(state wL, state wR, double pStar){
                return (wL.u+wR.u + f(pStar,wR)-f(pStar,wL))/2.0;
            }
            double solveRhoStar(state w, double p){
                if(p>w.p) return w.rho * (p+G6*w.p)/(G6*p+w.p);
                else return w.rho * pow( p/w.p, 1.0/gamma);
            }
            void ExactRiemann(state wL, state wR, state* starL, state* starR){
                starL->p=solvePStar(wL, wR); starR->p = starL->p;
                starL->u=solveUStar(wL, wR, starL->p); starR->u = starL->u;
                starL->rho = solveRhoStar(wL, starL->p);
                starR->rho = solveRhoStar(wR, starR->p);
                starL->v = wL.v; starR->v = wR.v;
            }
            //Riemann Sampler
            state riemannSample(state w, state s, double x_t){
                double a = sqrt(gamma * w.p / w.rho);
                int zone = 0; // 1 = outside, 2 = star, 3 = fan
                if (s.p > w.p){ //shock
                    double scale = sqrt( G2*s.p/w.p + G1 ), sign = (x_t> s.u)?1:-1;
                    std::cout<<sign*(x_t-w.u)  - scale * a<<"\n";
                    zone = ( sign*(x_t-w.u)  > scale * a ) ? 1 : 2;
                } else{ //Rarefraction
                    double scale = pow(s.p/w.p, G1), sign = (x_t> s.u)?1:-1;
                    if ( (x_t-w.u)*sign > a ) zone = 1;
                    else if ( fabs(x_t-s.u)  < scale * a ) zone = 2;
                    else zone = 3;
                }
                state wRet;
                switch(zone){
                    case 1: wRet=w; break;
                    case 2: wRet=s;break;
                    case 3:
                        double as= (x_t > s.u ? a : -a);
                        double scale = G5 - G6 * (w.u-x_t)/as;
                        wRet.rho = w.rho * pow(scale, G4);
                        wRet.p = w.p * pow(scale,G3);
                        wRet.u = G5 * (x_t-as + w.u * G7);
                        break;
                }
                //Debug
               //wRet.rho = zone; wRet.u = (x_t > u ? 1 : -1); wRet.p = (p > w.p ? 1:-1);
                wRet.v = w.v;
                return wRet;
            }

        //Godunov
        stateConsv flux(state wL, state wR){
            //Using Roe return roe(wL, wR);
            //Using Exact Riemann
            state starL, starR;
            ExactRiemann(wL, wR, &starL, &starR);
           // ANRS(wL, wR, &starL, &starR);
            state w = riemannSample( (0 > starL.u ? wR : wL), (0 > starR.u ? starR : starL), 0);
            return flux(conservativeForm(w), w.u);
        }


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

void advance1D(state* w, int size, double time){
    stateConsv fL = flux(w[0], w[0]), fR;
    for(int i=0; i<size; i++){
        state wR = i==size-1 ? w[i] : w[i+1];
        stateConsv u = conservativeForm(w[i]);
        fR = tvdFlux(w, i, size, time);
        u.rho += (fL.rho-fR.rho)*(time / stepSize);
        u.px += (fL.px - fR.px)*(time / stepSize);
        u.py += (fL.py - fR.py)*(time / stepSize);
        u.E += (fL.E - fR.E)*(time / stepSize);
        w[i] = primativeForm(u);
        fL = fR;
    }
}

void advanceX(state** w, int sizex, int sizey, double time){
    for(int i=0; i<sizey; i++) advance1D(w[i], sizex, time);
}

void advanceY(state** w, int sizex, int sizey, double time){
    state  wTemp[sizey];
    for(int i=0; i<sizex; i++) {
        for(int j=0; j<sizey; j++) wTemp[j] = transform(w[j][i]);
        advance1D(wTemp, sizey, time);
        for(int j=0; j<sizey; j++) w[j][i] = transform(wTemp[j]);
    }
}



void advance2D(state** w, int sizex, int sizey, double time){
    if (time <= 0)  return;
    double timeStep = time;
    for(int i=0; i<sizey; i++){for(int j=0; j<sizex; j++){
        timeStep = fmin(timeStep, 0.4*stepSize/S(w[i][j]));
    }}
    if (timeStep != timeStep || fabs(timeStep) < 1E-6) timeStep = 1E-6;
    
    advanceX(w, sizex, sizey, timeStep/2);
    advanceY(w, sizex, sizey, timeStep);
    advanceX(w, sizex, sizey, timeStep/2);

    advance2D(w, sizex,sizey, time-timeStep);
}



const double x0 = 0.5;
const double tFinal = 0.15;

int main(int argc, const char * argv[]) {
    state wL; wL.rho=2.0; wL.u = -1E5; wL.p = 2*2.475E10;
    state wR; wR.rho=1.0; wR.u= 2E5; wR.p = 2.475E10;
    state starLE, starRE, wLE = wL, wRE = wR;
    ExactRiemann(wLE, wRE, &starLE, &starRE);
    std::cout <<  ( wL.u + wR.u + f(starLE.p,wL) - f(starRE.p,wR))/2.0 <<"\n";
    std::cout <<  ( starLE.p) <<"\n";

    state ans = riemannSample(wR, starRE, 980.0*60/2.0);
    std::cout<< starLE.rho << ","<< starLE.u << ","<< starLE.p<<"\n";
    std::cout<< starRE.rho << ","<< starRE.u << ","<< starRE.p<<"\n";
    
    std::cout<< ans.rho << ","<< ans.u << ","<< ans.p<<"\n";
    std::cout<< (ans.rho*(ans.u - 981*60))<<"\n";
    return 0;
    
    
  //  wL.u /= sqrt(2); wL.v = wL.u; wR.v = wR.u;
    /* Riemann Solver
    state starL, starR;
    ANRS(wL, wR, &starL, &starR);
    for(int i=0; i<101; i++){
        double x_t = (i*stepSize - x0)/tFinal;
        
        state w = riemannSample( (x_t > starL.u ? wR : wL), (x_t >  starL.u ? starR : starL),  x_t);
        state wE = riemannSample( (x_t > starLE.u ? wR : wL),  (x_t >  starLE.u ? starRE : starLE),  x_t);

        std::cout << i*stepSize<<","<<w.rho<<","<<wE.rho<<","<<w.u<<","<<wE.u<<","<<w.p<<","<<wE.p<<"\n";
    }
     //*/
    //* Godunov Method
    int sizex=100, sizey=100;
    state** w = (state**) malloc(sizeof(state*) * sizey);
    for(int j=0;j<sizey;j++){
        w[j] = (state*) malloc(sizeof(state) * sizex);
        for(int i=0;i<sizex;i++) w[j][i] = (i+j)/2 < x0/stepSize ? wL : wR;
    }
    
    for(double t = 0; t<tFinal;t+=timeStep) advance2D(w, sizex,sizey, timeStep);
    for(int i=0;i<100;i++){
        double x_t = sqrt(2)*(i*stepSize-x0) / tFinal;
        state exact = riemannSample( (x_t > starLE.u ? wRE : wLE),  (x_t >  starLE.u ? starRE : starLE),  x_t);
        //std::cout << x_t<<","<<w.rho<<","<<w.u<<","<<w.p<<"\n";

        std::cout << i*stepSize<<","<<w[i][i].rho<<","<<exact.rho<<","<<w[i][i].u<<","<<exact.u<<","<<w[i][i].p<<","<<exact.p<<"\n";
    }
    //*/
    return 0;
}
