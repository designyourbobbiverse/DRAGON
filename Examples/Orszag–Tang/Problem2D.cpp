//
//  Problem2D.cpp
//  DRAGON/Examples/Orszag–Tang
//
//  Created by Bobbie Markwick on 22/07/2026.
//

#include "Problem.hpp"
#include "DistGrid.hpp"

#include <cmath> //For std::sin, cos

typedef DistGrid2D MyGrid;//Choose the dimension of your grid here

constexpr double rho0 = 25.0 / (36*M_PI);
constexpr double p_amb = 5.0 / (12*M_PI);
constexpr double B0 = 1.0;

constexpr int n = 512;
constexpr double dx = 2*M_PI/n;


Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(n,n, dx, dx);
    grid->boundary = Boundary::Periodic();
    return *grid;
}

PrimitiveState Problem::initialFluidState(double x, double y, double z){
    //Initialize the fluid state w at point (x,y,z).
        //dx/2, dy/2 corresponds to the [0,0] cell. z will always be zero.
    PrimitiveState w;
    w.rho =  rho0;
    w.p = p_amb;
    w.v = {-std::sin(y), std::sin(x),0};
    return w;
}
vec3 Problem::initialMagneticPotential(double x, double y, double z){
    //Initialize the vector potential at point (x,y,z)
        // (0,0,0) corresponds to the [0,0,0] cell. z will always be zero in 2D
    //Magnetic Fields will be initialized from this potential to ensure div B = 0
    double Az = B0 * (std::cos(y) + 0.25*std::cos(2*x));
    return {0,0,Az};
}


void Problem::completeProblemInit(Grid& problem){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    //Here you can do any initialization not covered by initialFluidState and initialMagneticPotential
    
}


void Problem::beforeCycle(Grid &problem, int cycle, double t){
    MyGrid& grid = *dynamic_cast<MyGrid*>(&problem);
    
    //Do any processing you need to do between cycles
    //This is called before the Godunov scheme. It isn't called for the initial frame
    problemComplete(problem, t);
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
