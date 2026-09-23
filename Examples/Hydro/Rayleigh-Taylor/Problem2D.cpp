//
//  Problem2D.cpp
//  DRAGON/Examples/Hydro/Rayleigh-Taylor
//
//  Created by Bobbie Markwick on 23/09/2026.
//

#include "main/Problem.hpp"
#include "Refinement/DistGrid.hpp"
#include "Source/Sources.hpp"
using namespace DRAGON;

#include <cmath> //For std::sin etc
#include "Constants.h" //For _pi

typedef DistGrid2D MyGrid;//Choose the dimension of your grid here

constexpr double rho1   = 2.0;     // density, upper layer
constexpr double rho2   = 1.0;     // density, lower layer
constexpr double g      = 3.0;     // gravity
constexpr double p0     = rho1*g;     // pressure at the boundary


constexpr double sigma  = 0.2;     // Gaussian localization width for vy seed
constexpr double amp    = 0.025;    // vy perturbation amplitude
constexpr int    kmode  = 2;       // wavenumber (integer, box-periodic in x)

constexpr int n = 512;



Grid& Problem::makeProblem(){
    //Construct your grid object. Don't worry about initial setup, you'll do that later
    auto grid = new MyGrid(n,n, 1.0/n,1.0/n);
    grid->boundary = Boundary::Periodic("X");
    
    grid->passives().add("Lower");
    grid->passives().add("Upper");
    
    grid->sources = Source::UniformGravity(0,-g,0);

    return *grid;
}


PrimitiveState Problem::initialFluidState(double x, double y, double z){
    //Initialize the fluid state w at point (x,y,z).
        //(dx/2,dy/2) corresponds to the [0,0] cell. As such, we need to convert
    double y0 = 0.5 + amp * std::cos(2.0 * _pi * kmode * x);
    
    PrimitiveState w{};
    w.rho = y>y0 ? rho1 : rho2;
    if (y>0.5) {
        w.p = p0 - rho1 * (y-0.5) * g;
    } else {
        w.p = p0 - rho2 * (y-0.5) * g;
    }
    
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
    
    for(int i=0; i<grid.getSizeX(); i++){
        double y0 = 0.5 + amp * std::cos(2.0 * _pi * kmode * (i+0.5)*grid.dx);
        for(int j=0; j<grid.getSizeY(); j++){
            double y = (j+0.5)*grid.dy;
            grid.passives()[i,j,"Lower"] = (y <= y0 ? 1 : 0);
            grid.passives()[i,j,"Upper"] = (y >= y0 ? 1 : 0);
        }
    }
    
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
