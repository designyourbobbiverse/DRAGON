//
//  Problem3D.cpp
//  DRAGON/Examples/MHD/Blast
//
//  Created by Bobbie Markwick on 8/07/2026.
//  Gardiner and Stone (2008). https://arxiv.org/abs/0712.2634
//

#include "main/Problem.hpp"
#include "Refinement/DistGrid.hpp"
using namespace DRAGON;

#include <cmath>        //For std::sqrt
#include "Constants.h"  //For sq4pi

typedef DistGrid3D MyGrid;//Choose the dimension of your grid here

constexpr double rho = 1.0;

constexpr double p_amb = 1;
constexpr double p_blast = 100;
constexpr double r0 = 0.125;

const double B0 = 10 * sq4pi;

constexpr bool diagonal = true;

constexpr int n = 256;



Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(n,n,n, 1.0/n,1.0/n,1.0/n);
    grid->boundary = Boundary::Outflow();
    
    return *grid;
}


PrimitiveState Problem::initialFluidState(double x, double y, double z){
    //Initialize the fluid state w at point (x,y,z).
        //(dx/2,dy/2,dz/2) corresponds to the [0,0,0] cell. As such, we need to convert
    x -= 0.5; y-=0.5; z-=0.5;
    double r = std::sqrt(x*x + y*y + z*z);

    PrimitiveState w;
    w.rho = rho;
    w.v = {0,0,0};
    w.p = r < r0 ? p_blast : p_amb;
    return w;
}
vec3 Problem::initialMagneticPotential(double x, double y, double z){
    //Initialize the vector potential at point (x,y,z)
        // (0,0,0)corresponds to the [0,0,0] cell. As such, we need to convert
        x -= 0.5; y-=0.5; z-=0.5;
    //Magnetic Fields will be initialized from this potential to ensure div B = 0
    return  diagonal ? vec3{0, B0 * (x-z)/std::sqrt(2) ,0} : vec3{0, B0 * x, 0 };
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
