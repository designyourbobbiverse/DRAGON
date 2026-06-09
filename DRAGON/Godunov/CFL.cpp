//
//  CFL.cpp
//  DRAGON
//
//  Created by Bobbie Markwick on 09/06/2026.
//
#include "CFL.hpp"
#include "Config.h"
#include "Constants.h"
#include <math.h>

//MARK: Individual Speeds
double CFL::cfl_max_speed(const PrimitiveState& W, double dx, double dy, double dz){
    double a = sqrt(_gamma * W.p / W.rho);
    double speed = 0;
    if (dx > 1e-14) speed = (fabs(W.vx) + a)/dx;
    if (dy > 1e-14)  speed = std::max(speed, (fabs(W.vy) + a)/dy);
    if (dz > 1e-14) speed = std::max(speed, (fabs(W.vz) + a)/dz);
    return speed;
}
double CFL::cfl_add_speed(const PrimitiveState& W, double dx, double dy, double dz){
    double a = sqrt(_gamma * W.p / W.rho);
    double speed = 0;
    if (dx > 1e-14) speed = (fabs(W.vx) + a)/dx;
    if (dy > 1e-14) speed +=  (fabs(W.vy) + a)/dy;
    if (dz > 1e-14) speed += (fabs(W.vz) + a)/dz;
    return speed;
}
double CFL::cfl_pow_speed(const PrimitiveState& W, double p, double dx, double dy, double dz){
    double a = sqrt(_gamma * W.p / W.rho);
    double speed = 0;
    if (dx > 1e-14) speed = pow((fabs(W.vx) + a)/dx, 2);
    if (dy > 1e-14) speed +=  pow((fabs(W.vy) + a)/dy,2);
    if (dz > 1e-14) speed += pow((fabs(W.vz) + a)/dz,2);
    return pow(speed, 1.0/p);
}

//MARK: Configuration-Control
#if CFL_CALCULATION == CHOOSE_RUNTIME || defined TESTMODE
namespace CONFIG {
    int cflChoice = CFL_ADD;
}
#endif

double CFL::cfl_speed(const PrimitiveState& W, double dx, double dy, double dz){
#if CFL_CALCULATION == CHOOSE_RUNTIME || defined TESTMODE
    switch(CONFIG::cflChoice){
        case CFL_MAX: return cfl_max_speed(W, dx, dy);
        case CFL_ADD: return cfl_add_speed(W, dx, dy);
        case CFL_V2: return cfl_pow_speed(W, 2, dx, dy);
        default: return cfl_pow_speed(W, CONFIG::cflChoice, dx, dy);
    }
#elif CFL_CALCULATION == CFL_MAX
    return cfl_max_speed(W, dx, dy, dz);
#elif CFL_CALCULATION == CFL_ADD
    return cfl_add_speed(W, dx, dy, dz);
#elif CFL_CALCULATION == CFL_V2
    return cfl_pow_speed(W, 2, dx, dy, dz);
#endif
}

//MARK: Grid Minimum
double CFL::cfl_time(const Grid1D& g){
    double max_speed = 0;
    for(int i = 0; i<g.getSize(); i++){
      const PrimitiveState& W = g[i];
      double speed = fabs(W.vx) + sqrt(_gamma * W.p / W.rho);
      max_speed = fmax(max_speed, speed);
    }
    if(max_speed <= 0) throw std::runtime_error("CFL timestep failed: max signal speed is zero");
    return CFL_coeff * g.dx / max_speed;
}

double CFL::cfl_time(const Grid2D& g){
    double max_speed = 0;
    for(int i = 0; i<g.getSizeX(); i++){
        for(int j = 0; j<g.getSizeY(); j++){
            double speed = cfl_speed(g[i,j], g.dx, g.dy);
            max_speed = fmax(max_speed, speed);
        }
    }
    if(max_speed <= 0) throw std::runtime_error("CFL timestep failed: max signal speed is zero");
    return CFL_coeff / max_speed;
}

double CFL::cfl_time(const Grid3D& g){
    double max_speed = 0;
    for(int i = 0; i<g.getSizeX(); i++){
        for(int j = 0; j<g.getSizeY(); j++){
            for(int k = 0; k<g.getSizeZ(); k++){
                double speed = cfl_speed(g[i,j,k], g.dx, g.dy, g.dz);
                max_speed = fmax(max_speed, speed);
            }
        }
    }
    if(max_speed <= 0) throw std::runtime_error("CFL timestep failed: max signal speed is zero");
    return CFL_coeff / max_speed;
}
