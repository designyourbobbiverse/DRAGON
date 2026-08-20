//
//  Problem2D.cpp
//  DRAGON/Examples/Hydro/Kelvin-Helmholtz
//
//  Created by Bobbie Markwick on 8/07/2026.
//

#include "Problem.hpp"
#include "DistGrid.hpp"
using namespace DRAGON;

#include <cmath> //For std::sin etc

typedef DistGrid2D MyGrid;//Choose the dimension of your grid here


constexpr double width   = 0.1;
constexpr double y_lower = 0.5 - width/2;
constexpr double y_upper = 0.5 + width/2;
constexpr double a       = 0.05;    // tanh smoothing width (interface)

constexpr double rho1   = 2.0;     // density, layer y_lower and y_upper
constexpr double rho2   = 1.0;     // density, outer layer
constexpr double V0     = 0.5;         // Veloicty of the outer layer (= - velocity of inner layer
constexpr double p0     = 2.5;     // uniform initial pressure


constexpr double sigma  = 0.2;     // Gaussian localization width for vy seed
constexpr double amp    = 0.01;    // vy perturbation amplitude
constexpr int    kmode  = 2;       // wavenumber (integer, box-periodic in x)

constexpr int n = 512;



Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(n,n, 1.0/n,1.0/n);
    grid->boundary = Boundary::Outflow();
    
    return *grid;
}


PrimitiveState Problem::initialFluidState(double x, double y, double z){
    //Initialize the fluid state w at point (x,y,z).
        //(dx/2,dy/2) corresponds to the [0,0] cell. As such, we need to convert
    PrimitiveState w;
    w.p = p0;
   
    
    
    const double s1 = 0.5 * (1.0 + std::tanh((y - y_lower) / a));
    const double s2 = 0.5 * (1.0 + std::tanh((y - y_upper) / a));
    const double layerFrac = s1 * (1.0 - s2);
     
    w.rho = rho2 + (rho1 - rho2) * layerFrac;
    const double vx = V0 - 2 * V0 * (1-layerFrac);
    
    //Perturbation
    const double env1 = std::exp(-(y - y_lower) * (y - y_lower) / (2.0 * sigma * sigma));
    const double env2 = std::exp(-(y - y_upper) * (y - y_upper) / (2.0 * sigma * sigma));
    const double vy = amp * std::sin(2.0 * M_PI * kmode * x) * (env1 + env2);

    w.v = {vx, vy, 0};
    
    return w;
}
vec3 Problem::initialMagneticPotential(double x, double y, double z){
    //This funciton isn't called in Pure Hydro
    return {0,0,0};
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
