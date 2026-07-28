//
//  Problem2D.cpp
//  DRAGON/Examples/LoopAdvection
//
//  Created by Bobbie Markwick on 20/07/2026.
//

#include "Problem.hpp"
#include "DistGrid.hpp"
#include <cmath>
#include <iostream>

typedef DistGrid2D MyGrid;//Choose the dimension of your grid here

constexpr double rho0 = 1.0;
constexpr double p_amb = 5.0;
constexpr double vx = 2.0;
constexpr double vy = 0.5;

constexpr bool gaussian = false;
constexpr double r0 = gaussian ? 0.1 : 0.2;
constexpr double B0 = 1e-6;

constexpr int n = 256;
constexpr double dx = 2.0/n;


Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(n,n/2, dx, dx);
    grid->boundary = Boundary::Periodic();
    return *grid;
}

PrimitiveState Problem::initialFluidState(double x, double y, double z){
    //Initialize the fluid state w at point (x,y,z).
        //dx/2, dy/2 corresponds to the [0,0] cell. z will always be zero.
    PrimitiveState w;
    w.rho =  rho0;
    w.p = p_amb;
    w.v = {vx,vy,0};
    return w;
}
vec3 Problem::initialMagneticPotential(double x, double y, double z){
    //Initialize the vector potential at point (x,y,z)
        // (0,0,0) corresponds to the [0,0,0] cell. z will always be zero in 2D
    //As such, we need to adjust them to be center-relative
    x -= 1.0;
    y -= 0.5;
    //Magnetic Fields will be initialized from this potential to ensure div B = 0
    double r = sqrt(x*x + y*y);
    double Az = gaussian ? B0*exp(-0.5*pow(r/r0,2)) : B0*fmax(r0-r,0);
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
    double L1 = 0;
    double L2 = 0;
    double Linf = 0;
    
    double cx = fmod(vx * t + 1.0, 2.0);
    double cy = fmod(vy * t + 0.5, 1.0);

    
    for(int i=0; i<n;i++){
        for(int j=0; j<n/2; j++){
            double x = (i+0.5) * dx - cx;
            double y = (j+0.5) * dx - cy;
            double r = sqrt(x*x + y*y);
            double B2_exact = r <= r0 ? B0*B0 : 0;

            double err = fabs(grid[i,j].B*grid[i,j].B - B2_exact) / (B0*B0);
            if(err > Linf) Linf = err;
            L1 += err ;
            L2 += err*err;
            
        }
    }
    L2 = sqrt(L2 / (n*n));
    
    std::cout<<"L1 error: "<<L1 / (n*n)<<"\n";
    std::cout<<"L2 error: "<<L2<<"\n";
    std::cout<<"L-infinity error: "<<Linf<<"\n";
}
