//
//  Problem2D.cpp
//  DRAGON/Examples/Hydro/Hydrostatic-Equilibrium
//
//  Created by Bobbie Markwick on 22/09/2026.
//

#include "main/Problem.hpp"
#include "Refinement/DistGrid.hpp"
#include "Source/Sources.hpp"
using namespace DRAGON;


typedef DistGrid2D MyGrid;//Choose the dimension of your grid here

constexpr double rho0 = 1;
constexpr double p0 = 1;
constexpr double g = 10;

constexpr int n = 256;



Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(n/2,n, 1.0/n,1.0/n);
    grid->boundary = Boundary::Periodic("X") + Boundary::Reflective("Y");
    
    grid->sources = Source::UniformGravity(0,-g,0);
    
    return *grid;
}


PrimitiveState Problem::initialFluidState(double x, double y, double z){
    //Initialize the fluid state w at point (x,y,z).
        //(dx/2,dy/2) corresponds to the [0,0] cell. As such, we need to convert
    PrimitiveState w;
    w.rho = rho0;
    w.p = p0 + (1-y)*rho0*g;
    
    return w;
}
vec3 Problem::initialMagneticPotential(double x, double y, double z){
    //This funciton isn't called in Pure Hydro
    return {0,0,0};
}


void Problem::completeProblemInit(Grid& problem){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    //Here you can do any initialization not covered by initialFluidState and initialMagneticPotential
    //For example, you can initialise passive scalars here
    
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
