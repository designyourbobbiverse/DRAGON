//
//  Problem2D.cpp
//  DRAGON/Examples/MHD/MHDRotor
//
//  Created by Bobbie Markwick on 8/07/2026.
//  Balsara and Spicer (1999). https://doi.org/10.1006/jcph.1998.6153
//  Toth (2000). https://doi.org/10.1006/jcph.2000.6519
//

#include "main/Problem.hpp"
#include "Refinement/DistGrid.hpp"
using namespace DRAGON;

#include <cmath>        //For std::sqrt

typedef DistGrid2D MyGrid;//Choose the dimension of your grid here

constexpr double rho0 = 10.0;
constexpr double rho_amb = 1.0;
constexpr double p_amb = 1;
constexpr double omega = 20.0;
const double B0 = 5;

constexpr double r0 = 0.1;
constexpr double r1 = 0.115;

constexpr int n = 256;



Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(n,n, 1.0/n,1.0/n);
    grid->boundary = Boundary::Outflow();
    
    return *grid;
}


PrimitiveState Problem::initialFluidState(double x, double y, double z){
    //Initialize the fluid state w at point (x,y,z).
        //(dx/2,dy/2) corresponds to the [0,0] cell. As such, we need to convert
    x -= 0.5; y-=0.5;
    double r = std::sqrt(x*x + y*y);
    
    PrimitiveState w;
    w.p = p_amb;
    
    if(r<r0){
        w.rho = rho0;
        w.v = omega * vec3{-y,x,0};
    } else if (r < r1){
        double f = (r1 - r) / (r1 - r0);
        w.rho = rho_amb + (rho0-rho_amb) *  f;
        w.v = (omega * f * r0/r) * vec3{-y,x,0};
    } else {
        w.rho = rho_amb;
        w.v = {0,0,0};
    }
    
    return w;
}
vec3 Problem::initialMagneticPotential(double x, double y, double z){
    //Initialize the vector potential at point (x,y,z)
        // (0,0) corresponds to the [0,0] corner. As such, we need to convert
        x -= 0.5; y-=0.5;
    //Magnetic Fields will be initialized from this potential to ensure div B = 0
    return  vec3{0, 0, B0 * y};
}


void Problem::completeProblemInit(Grid& problem){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    //Here you can do any initialization not covered by initialFluidState and initialMagneticPotential
    
}

void Problem::beforeCycle(Grid &problem, int cycle, double t){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    
    //Do any processing you need to do between cycles
    //This is called before the Godunov scheme. It isn't called for the initial frame
    
}


void Problem::afterCycle(Grid &problem, int cycle, double t){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    
    //Do any processing you need to do between cycles
    //This is called after the Godunov scheme but before the file write. It isn't called for the initial frame
    
}

void Problem::problemComplete(Grid& problem, double t){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);

    //This is called only after the final time is reached.
}
