//
//  CFL.cpp
//  DRAGON/Hydro/CFL
//
//  Created by Bobbie Markwick on 09/06/2026.
//  Implementation based in part on Toro (2009). https://doi.org/10.1007/b79761
//      Courant, Friedrichs, and Lewy (1928). https://doi.org/10.1007/BF01448839
//

#include "CFL.hpp"
#include "Config.h"
#include "Constants.h"
#include <math.h>
#include <iostream>
using namespace CONFIG;

//MARK: Individual Speeds
//CFL speed = highest (|v|+a)/dL of any one side
double CFL::cfl_max_speed(const PrimitiveState& W, double dx, double dy, double dz){
#ifdef MHD //Use Magnetosonic fast speed for MHD
    double a = W.c_fast_max();
#else //Use sound speed for Hydro
    double a = W.cs();
#endif
    double speed = 0;
    if (dx > 1e-14) speed = (fabs(W.v.x) + a)/dx;
    if (dy > 1e-14)  speed = std::max(speed, (fabs(W.v.y) + a)/dy);
    if (dz > 1e-14) speed = std::max(speed, (fabs(W.v.z) + a)/dz);
    return speed;
}

//CFL speed = sum of (|v|+a)/dL over all applicable sides
double CFL::cfl_add_speed(const PrimitiveState& W, double dx, double dy, double dz){
#ifdef MHD //Use Magnetosonic fast speed for MHD
    double a = W.c_fast_max();
#else //Use sound speed for Hydro
    double a = W.cs();
#endif
    double speed = 0;
    if (dx > 1e-14) speed = (fabs(W.v.x) + a)/dx;
    if (dy > 1e-14) speed +=  (fabs(W.v.y) + a)/dy;
    if (dz > 1e-14) speed += (fabs(W.v.z) + a)/dz;
    return speed;
}

//CFL speed = [sum of ((|v|+a)/dL)^p]^(1/p)
double CFL::cfl_pow_speed(const PrimitiveState& W, double p, double dx, double dy, double dz){
#ifdef MHD //Use Magnetosonic fast speed for MHD
    double a = W.c_fast_max();
#else //Use sound speed for Hydro
    double a = W.cs();
#endif
    double speed = 0;
    if (dx > 1e-14) speed = pow((fabs(W.v.x) + a)/dx, p);
    if (dy > 1e-14) speed +=  pow((fabs(W.v.y) + a)/dy,p);
    if (dz > 1e-14) speed += pow((fabs(W.v.z) + a)/dz,p);
    return pow(speed, 1.0/p);
}

//MARK: Configuration-Control
#if CFL_CALCULATION == CHOOSE_RUNTIME || defined(TESTMODE)
namespace CONFIG {
    int cfl_choice = CFL_ADD;
}
#endif

double CFL::cfl_speed(const PrimitiveState& W, double dx, double dy, double dz){
#if CFL_CALCULATION == CHOOSE_RUNTIME || defined(TESTMODE) //Choose at runtime (which is always the case when unit testing)
    switch(CONFIG::cfl_choice){
        case CFL_MAX: return cfl_max_speed(W, dx, dy, dz);
        case CFL_ADD: return cfl_add_speed(W, dx, dy, dz);
        default: return cfl_pow_speed(W, CONFIG::cfl_choice, dx, dy, dz);
    }
#elif CFL_CALCULATION == CFL_MAX //Use highest individual dimension
    return cfl_max_speed(W, dx, dy, dz);
#elif CFL_CALCULATION == CFL_ADD //Add the speeds together
    return cfl_add_speed(W, dx, dy, dz);
#elif CFL_CALCULATION > 0 //Use the  (sum of []^p)^(1/p) method
    return cfl_pow_speed(W, CFL_CALCULATION, dx, dy, dz);
#endif
}

//MARK: Grid Minimum
//Calculate CFL time limit over the entire grid
double CFL::cfl_time(const Grid1D& g){
    //Calculate the highest speed over the entire grid
    double max_speed = 0;
    int ng = std::min(1,g.getGhosts());
    for(int i = -ng; i<g.getSize()+ng; i++){
        const PrimitiveState& W = g[i];
        #ifdef MHD //Use Magnetosonic fast speed for MHD
        double a = W.c_fast_max();
        #else //Use sound speed for Hydro
        double a = W.cs();
        #endif
      double speed = fabs(W.v.x) + a;
      max_speed = fmax(max_speed, speed);
    }
    if(max_speed <= 0) throw std::runtime_error("CFL timestep failed: max signal speed is zero");
    //Then convert to time and apply the coefficient
    return CFL_coeff * g.dx / max_speed;
}

double CFL::cfl_time(const Grid2D& g){
    //Calculate the highest speed over the entire grid
    double max_speed = 0;
    int ng = std::min(1,g.getGhosts());
    for(int i = -ng; i<g.getSizeX()+ng; i++){
        for(int j = -ng; j<g.getSizeY()+ng; j++){
            double speed = cfl_speed(g[i,j], g.dx, g.dy);
            max_speed = fmax(max_speed, speed);
        }
    }
    if(max_speed <= 0) throw std::runtime_error("CFL timestep failed: max signal speed is zero");
    //Then convert to time and apply the coefficient
    return CFL_coeff / max_speed;
}

double CFL::cfl_time(const Grid3D& g){
    //Calculate the highest speed over the entire grid
    double max_speed = 0;
    int ng = std::min(1,g.getGhosts());
    for(int i = -ng; i<g.getSizeX()+ng; i++){
        for(int j = -ng; j<g.getSizeY()+ng; j++){
            for(int k = -ng; k<g.getSizeZ()+ng; k++){
                double speed = cfl_speed(g[i,j,k], g.dx, g.dy, g.dz);
                max_speed = fmax(max_speed, speed);
            }
        }
    }
    if(max_speed <= 0) throw std::runtime_error("CFL timestep failed: max signal speed is zero");
    //Then convert to time and apply the coefficient
    return CFL_coeff / max_speed;
}
